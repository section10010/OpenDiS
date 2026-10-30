"""Native counterfactual force/velocity audit on actual evolved MT geometries."""
import argparse,json,sys,os
from pathlib import Path
import numpy as np
from run_shared import bootstrap,seed_from_npz
from material_identity import MaterialIdentity
from shared_contact import SharedContact


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--g0-root',type=Path,required=True)
    p.add_argument('--repo',required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--case',default='mesh2p5_timehalf_MT');a=p.parse_args();bootstrap(a.g0_root)
    from configuration import config,save_json
    from native import Simulation,imports
    if a.out.exists():raise RuntimeError('Refuse overwrite')
    aliases=json.loads((a.root/'case_aliases.json').read_text()).get('aliases',{}) if (a.root/'case_aliases.json').exists() else {}
    a.out.mkdir(parents=True);parent=a.root/'cases'/aliases.get(a.case,a.case);spec=json.loads((parent/'spec.json').read_text());c=config(**spec['common_native_updates'])
    sources=json.loads((parent/'sources.json').read_text());os.environ.update(OMP_NUM_THREADS='2',OMP_PROC_BIND='false');py,_=imports(a.repo);py.initialize(num_threads=2);checks=[]
    try:
        def centered(sim):
            if 'force_evaluation_origin_nm' in spec:
                from centered_force import install
                install(sim,spec['force_evaluation_origin_nm'])
            return sim
        paths=[parent/'initial.npz',*sorted(parent.glob('frame_*.npz')),parent/'final.npz']
        # All saved frames, not a favourable hand-selected instant.
        for path in paths:
            seed=seed_from_npz(path,c);sim=centered(Simulation(c,a.repo,initial=seed));sim.contact=SharedContact(MaterialIdentity(sources,spec['interface']))
            sim.contact.validate(sim,'force_audit');f,v=sim.evaluate();full_f=f.copy();full_v=v.copy();full_reaction=sim.velocity_diagnostic['reaction_N'].copy()
            aa=sim.arrays(False);meta=sim.contact.identity.inspect(aa);cross=[];single_force=np.zeros_like(f)
            for host in ('M','T'):
                node=meta['node_host']==host;edge=meta['segment_host']==host;indices=np.flatnonzero(node);remap=np.full(len(node),-1);remap[indices]=np.arange(len(indices))
                sub=dict(cell=seed['cell'],nodes={k:val[node] for k,val in seed['nodes'].items()},segs=dict(nodeids=remap[seed['segs']['nodeids'][edge]],burgers=seed['segs']['burgers'][edge],planes=seed['segs']['planes'][edge]))
                one=centered(Simulation(c,a.repo,initial=sub));one_f,_=one.evaluate();single_force[node]=one_f
                cross.append((full_f[node]-one_f)*sim.b**2)
            residual=float(np.linalg.norm(sum((d.sum(axis=0) for d in cross),np.zeros(3)))/sum(np.linalg.norm(d,axis=1).sum() for d in cross))
            # Direct cross-pair native reference avoids subtraction of two much
            # larger full fields. Do not enforce balance by projecting forces.
            pairs=[[int(i),int(j)] for i in np.flatnonzero(meta['segment_host']=='M') for j in np.flatnonzero(meta['segment_host']=='T')]
            reference=Simulation(c,a.repo,initial=seed)
            if 'force_evaluation_origin_nm' in spec:
                from centered_force import CenteredForce
                pair_net=CenteredForce(reference,reference.force,spec['force_evaluation_origin_nm']).network(reference.net).get_disnet(reference.exa.ExaDisNet)
            else:pair_net=reference.g
            direct=np.array(py.compute_force_segseglist(pair_net.net,mu=c['mu_pa'],nu=c['nu'],a=c['core_a_b'],segseglist=pairs))
            direct_residual=float(np.linalg.norm(direct.sum(axis=0))/np.linalg.norm(direct,axis=1).sum())
            pair_difference=float(np.linalg.norm((full_f-single_force)-direct)/np.linalg.norm(direct))
            # Deleting only cross-side elastic force on this identical native
            # geometry leaves constitutive law, graph, load and CTB unchanged.
            raw=np.array(py.compute_mobility(sim.g.net,mobility=sim.glide.mobility,nodeforces=single_force,nodetags=aa['tags']))
            removed=sim.contact.project_all(sim,sim.net,sim.state,raw)
            effect=float(np.linalg.norm((full_v-removed)*sim.b));assert effect>1e-12
            zero=centered(Simulation(config(**{**c,'stress_pa':np.zeros((3,3)).tolist()}),a.repo,initial=seed));zf,_=zero.evaluate()
            ij=aa['ids'];dr=(aa['positions_nm'][ij[:,1]]-aa['positions_nm'][ij[:,0]])*1e-9
            segment=.5*np.cross((np.array(c['stress_pa'])@(aa['burgers']*sim.b).T).T,dr);ext=np.zeros_like(full_f)
            for end in (0,1):np.add.at(ext,ij[:,end],segment)
            pk=float(np.linalg.norm((full_f-zf)*sim.b**2-ext)/np.linalg.norm(ext))
            np.savez_compressed(a.out/(path.stem+'_counterfactual.npz'),tags=aa['tags'],node_host=meta['node_host'],positions_nm=aa['positions_nm'],
                full_force_N=full_f*sim.b**2,without_cross_force_N=single_force*sim.b**2,direct_cross_pair_force_N=direct*sim.b**2,full_velocity_m_s=full_v*sim.b,
                without_cross_velocity_m_s=removed*sim.b,full_reaction_N=full_reaction,without_cross_reaction_N=sim.contact.reaction_native*sim.b**2)
            threshold=spec['criteria']['cross_force_balance_relative']
            checks.append(dict(snapshot=path.name,status='PASS' if max(residual,direct_residual)<threshold and pk<1e-10 else 'FAIL',
                cross_force_norm_N=float(np.linalg.norm((full_f-single_force)*sim.b**2)),cross_force_pair_residual=residual,
                direct_pair_residual=direct_residual,subtractive_vs_direct_relative=pair_difference,velocity_change_m_s=effect,PK_tensor_relative_error=pk))
        status='PASS' if all(c['status']=='PASS' for c in checks) else 'FAIL'
        save_json(a.out/'result.json',dict(status=status,checks=checks,frozen_balance_threshold=spec['criteria']['cross_force_balance_relative'],
            scope='All actual saved MT frames: N2 field subtraction AND direct native segment pairs; zero-stress PK; instantaneous cross-force deletion. No projection to enforce balance, no relaxed criterion, no long-time decoupled trajectory.'))
        print('G1_FORCE_AUDIT',status,len(checks))
    finally:py.finalize()


if __name__=='__main__':main()
