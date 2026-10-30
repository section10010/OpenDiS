"""Native positive, release, remesh, collision and rejection admission controls."""
import argparse,copy,json,sys,os,traceback
from pathlib import Path
import numpy as np
from run_shared import bootstrap,guarded_stage,snapshot
from shared_force_probe import initial_data
from material_identity import MaterialIdentity,network_arrays
from shared_contact import SharedContact


def main():
    p=argparse.ArgumentParser();p.add_argument('--g0-root',type=Path,required=True);p.add_argument('--repo',required=True)
    p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();bootstrap(a.g0_root)
    from configuration import config,save_json
    from geometry import build,residual
    from native import Simulation,imports
    if a.out.exists():raise RuntimeError('Refuse overwrite')
    a.out.mkdir(parents=True);spec=json.loads(a.spec.read_text());c=config(**spec['common_native_updates'])
    os.environ['OMP_NUM_THREADS']='2';py,_=imports(a.repo);py.initialize(num_threads=2);checks=[]
    def make(seed,sources,cc=c):
        sim=Simulation(cc,a.repo,initial=seed);sim.repo=a.repo
        if 'force_evaluation_origin_nm' in spec:
            from centered_force import install
            install(sim,spec['force_evaluation_origin_nm'])
        sim.contact=SharedContact(MaterialIdentity(sources,spec['interface']),cc['contact_tol_nm'])
        sim.contact.validate(sim,'initialization');return sim
    try:
        seed,_,sources=initial_data(c,'MT',build,spec);sim=make(seed,sources)
        # Move each complete source along its own affine slip plane to contact.
        touch=copy.deepcopy(seed);pos=touch['nodes']['positions']*sim.bnm
        for source in sources:
            tags={tuple(t) for t in source['initial_global_tags']};idx=np.array([i for i,t in enumerate(touch['nodes']['tags']) if tuple(t) in tags])
            sign=-1 if source['region']=='M' else 1;n=np.array(source['selected_slip']['n'])
            direction=np.array([0.,0.,1.])-n*n[2];direction/=np.linalg.norm(direction)
            gap=np.min(sign*(pos[idx,2]-spec['interface']['z_nm']))
            pos[idx]+=(-sign*gap/direction[2])*direction
        touch['nodes']['positions']=pos/sim.bnm;sim=make(touch,sources);f,v=sim.evaluate();d=sim.velocity_diagnostic
        raw=native=np.array(py.compute_mobility(sim.g.net,mobility=sim.glide.mobility,nodeforces=f,nodetags=sim.state['nodeveltags']))*sim.b
        np.testing.assert_allclose(raw,d['raw_velocity_m_s'],rtol=1e-12,atol=1e-8)
        reconstruction=np.array(py.compute_mobility(sim.g.net,mobility=sim.glide.mobility,nodeforces=f+d['reaction_N']/sim.b**2,nodetags=sim.state['nodeveltags']))*sim.b
        np.testing.assert_allclose(reconstruction,d['constrained_velocity_m_s'],rtol=1e-10,atol=1e-7)
        aa=sim.arrays(False);meta=sim.contact.identity.inspect(aa);gap=np.abs(aa['positions_nm'][:,2]-spec['interface']['z_nm'])
        for host in ('M','T'):
            i=int(np.flatnonzero((meta['node_host']==host)&(gap<1e-6))[0]);tag=aa['tags'][i];sign=-1 if host=='M' else 1
            assert sign*raw[i,2]<0 and np.linalg.norm(d['reaction_N'][i])>0
            assert abs(reconstruction[i,2])<1e-7
            inward=-f[i];vi=np.array(sim.glide.OneNodeMobility(sim.net,sim.state,tag,inward,False))
            vc=sim.contact.project_one(sim,sim.net,sim.state,tag,vi)
            assert sign*vi[2]>0;np.testing.assert_allclose(vc,vi,atol=1e-8)
            source=next(s for s in sources if s['region']==host);q=np.array(source['selected_slip']['q'])
            vt=np.array(sim.glide.OneNodeMobility(sim.net,sim.state,tag,q,False))
            np.testing.assert_allclose(sim.contact.project_one(sim,sim.net,sim.state,tag,vt),vt,rtol=1e-12,atol=1e-8)
            checks.append(dict(test=host+'_native_block_tangent_release_reaction',status='PASS'))
        snapshot(sim,sim.contact.identity,a.out/'contact_response.npz')
        # Native split must create real nodes; inherited source plane and host
        # are re-read rather than inferred from the newly allocated tag.
        fine=config(**{**c,'maxseg_nm':5.});remesh=make(seed,sources,fine);before=remesh.g.num_nodes()
        guarded_stage(remesh,'remesh',lambda:remesh.remesh.Remesh(remesh.net,remesh.state),a.out,0)
        after=remesh.g.num_nodes();assert after>before
        snapshot(remesh,remesh.contact.identity,a.out/'native_refined.npz')
        checks.append(dict(test='native_remesh_source_inheritance',status='PASS',nodes_before=before,nodes_after=after))
        coarse=config(**{**c,'maxseg_nm':20.,'minseg_nm':6.})
        coarsen=make(remesh.snapshot(),sources,coarse);coarse_before=coarsen.g.num_nodes()
        guarded_stage(coarsen,'remesh_coarsen',lambda:coarsen.remesh.Remesh(coarsen.net,coarsen.state),a.out,0)
        coarse_after=coarsen.g.num_nodes();assert coarse_after<coarse_before
        snapshot(coarsen,coarsen.contact.identity,a.out/'native_coarsened.npz')
        checks.append(dict(test='native_coarsen_source_inheritance',status='PASS',nodes_before=coarse_before,nodes_after=coarse_after))
        # Native tag reuse fixture: old M address now belongs to an actual T
        # segment. Origin is resolved from current edges, not stale tag storage.
        seedT,_,_=initial_data(c,'T',build,spec);seedT['nodes']['tags'][0]=seed['nodes']['tags'][0]
        reuse=make(seedT,sources);reuse.contact.identity.previous=copy.deepcopy(sim.contact.identity.previous)
        reuse.contact.validate(reuse,'native_address_reuse');assert reuse.contact.identity.inspect(reuse.arrays(False))['node_source_id'][0]==1
        checks.append(dict(test='native_recycled_address_material_identity',status='PASS'))
        # Real native Euler proposal intentionally violates the contact event
        # limit. Graph/time restoration and failed state are both retained.
        bad=make(seed,sources);bad.evaluate();caught=False
        try:guarded_stage(bad,'illegal_integration',lambda:bad.advance(20e-12),a.out,1)
        except RuntimeError:caught=True
        assert caught and (a.out/'rejected_0000001_illegal_integration/rejection.json').exists()
        checks.append(dict(test='actual_native_penetration_rejected_with_rollback',status='PASS'))
        # Put one vertex from each source at the common intersection of the two
        # affine slip planes with the CTB. ExaDiS may propose a cross-side merge;
        # I0 must reject it rather than silently create a material junction.
        collision_seed=copy.deepcopy(seed);A=[];rhs=[]
        for s in sources:
            n=np.array(s['selected_slip']['n']);A.append(n);rhs.append(n@np.array(s['selected_slip']['center_nm']))
        A.append([0,0,1]);rhs.append(spec['interface']['z_nm']);intersection=np.linalg.solve(np.array(A),np.array(rhs))
        indices=[int(np.flatnonzero(np.all(collision_seed['nodes']['tags']==s['initial_global_tags'][0],axis=1))[0]) for s in sources]
        collision_seed['nodes']['positions'][indices]=intersection/sim.bnm
        col=make(collision_seed,sources);snapshot(col,col.contact.identity,a.out/'cross_collision_initial.npz');n0=col.g.num_nodes();caught=False
        try:guarded_stage(col,'cross_material_collision',lambda:col.collision.HandleCol(col.net,col.state),a.out,2)
        except RuntimeError:caught=True
        n1=col.g.num_nodes()
        checks.append(dict(test='native_cross_material_collision',status='PASS' if caught else 'UNRESOLVED',
            initial_nodes=n0,proposal_nodes=n1,rejected=caught,note='No accepted cross-material reaction permitted in I0'))
        if not caught:raise RuntimeError('Native cross-material collision negative control did not exercise rejection')
        save_json(a.out/'result.json',dict(status='PASS',checks=checks,scope='Native contact, source identity and fail-closed I0 stage admission; no general CTB reactions'))
        print('G1_ADMISSION_PASS',len(checks))
    except Exception as exc:
        save_json(a.out/'result.json',dict(status='FAIL',checks=checks,error=repr(exc)));(a.out/'traceback.log').write_text(traceback.format_exc());raise
    finally:py.finalize()


if __name__=='__main__':main()
