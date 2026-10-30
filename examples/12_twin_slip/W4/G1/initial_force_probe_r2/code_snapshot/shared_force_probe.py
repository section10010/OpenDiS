"""G1 initial geometry/common-load/native elastic-coupling admission only.

One combined real ExaDiS network, not two separate-channel evolutions. This
probe intentionally refuses time integration until material-side contact and
native stage identity contracts have been implemented and qualified.
"""
import argparse,json,sys,os,shutil
from pathlib import Path
import numpy as np


def initial_data(c,population,build):
    rn=[];lk=[];region=[];sources=[];offset=0
    for label,phase,center in [('M','matrix',[1800.,2000.,1930.]),('T','twin',[2200.,2000.,2070.])]:
        local={**c,'phase':phase,'center_nm':center}
        points,links,src,_,_=build(local);links=links.copy();links[:,:2]+=offset
        source=src[0];selected={k:source[k] for k in ('phase','mode','slip_index','center_nm','b','n','q','p','rss_pa')}
        component_tags=[[0,offset+i] for i in range(len(points))]
        rn.extend(points);lk.extend(links);region.extend([label]*len(points));offset+=len(points)
        sources.append(dict(region=label,phase=phase,initial_component_id=0 if label=='M' else 1,
            initial_global_tags=component_tags,kind='propagation_loop_not_FR',selected_slip=selected,shared_interface_id='CTB0',
            allowed_normal=[0,0,-1 if label=='M' else 1]))
    rn=np.array(rn);lk=np.array(lk);region=np.array(region);mask=np.isin(region,list(population));old=np.flatnonzero(mask)
    remap=np.full(len(rn),-1);remap[old]=np.arange(len(old));keep=mask[lk[:,0].astype(int)]&mask[lk[:,1].astype(int)]
    seg=lk[keep];ids=remap[seg[:,:2].astype(int)];bnm=c['b_m']*1e9
    return dict(cell=dict(h=np.diag(np.array(c['box_nm'])/bnm),origin=np.zeros(3),is_periodic=c['periodic']),
        nodes=dict(tags=np.column_stack((np.zeros(len(old),int),old)),positions=rn[mask,:3]/bnm,constraints=rn[mask,3:4].astype(int)),
        segs=dict(nodeids=ids,burgers=seg[:,2:5],planes=seg[:,5:8])),region[mask],[x for x in sources if x['region'] in population]


def main():
    p=argparse.ArgumentParser();p.add_argument('--g0-root',type=Path,required=True);p.add_argument('--repo',required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--steps',type=int,default=0);a=p.parse_args()
    if a.steps!=0:raise NotImplementedError('G1 shared-side contact/identity stage integration NOT_RUN; no fallback to two W3 channels')
    gate=json.loads((a.g0_root/'analysis_r5/gate_review.json').read_text());assert gate['status']=='PASS' and gate['G1_admission']=='PASS'
    assert json.loads((a.g0_root/'A0/gate_result.json').read_text())['status']=='PASS'
    for filename in ('configuration.py','geometry.py','native.py'):
        assert (a.g0_root/'code'/filename).read_bytes()==(a.g0_root/'A0/production_snapshot'/filename).read_bytes()
    if a.out.exists():raise RuntimeError('Refuse overwrite')
    a.out.mkdir(parents=True);sys.path.insert(0,str(a.g0_root/'code'))
    (a.out/'code_snapshot').mkdir()
    for source in Path(__file__).resolve().parent.glob('*.py'):shutil.copy2(source,a.out/'code_snapshot'/source.name)
    from configuration import config,save_json
    from geometry import build,systems,Q,ROT,residual
    from native import Simulation,imports
    os.environ['OMP_NUM_THREADS']='2';py,_=imports(a.repo);py.initialize(num_threads=2);checks=[]
    try:
        sigma=np.diag([0.,0.,1100e6]);c=config(kind='propagation_loop',mode='HI',stress_pa=sigma.tolist(),
            integrator='EulerForward',vmax_m_s=None,Medge=32051.5,Mscrew=32051.5,walls=[],end_ns=.02)
        save_json(a.out/'common_native_config.json',c);networks={};forces={};fields={};metadata={}
        save_json(a.out/'rss_all24.json',{phase:[{**s,'n':s['n'].tolist(),'b':s['b'].tolist()} for s in systems(phase,sigma)] for phase in ('matrix','twin')})
        np.testing.assert_allclose((ROT@Q)@Q.T,ROT,atol=1e-14)
        for population in ('M','T','MT'):
            seed,regions,sources=initial_data(c,population,build);sim=Simulation(c,a.repo,initial=seed)
            f,v=sim.evaluate();aa=sim.arrays();forces[population]=f.copy()*sim.b**2;fields[population]=aa;networks[population]=sim
            np.savez_compressed(a.out/(population+'_initial.npz'),**aa,origin_region=regions)
            sim.g.write_data(str(a.out/(population+'_initial.data')))
            metadata[population]=dict(source_components=sources,origin_region=regions)
            assert sim.measures()['max_burgers_residual_b']<=1e-8
            for label in population:
                ix=np.flatnonzero(regions==label);sgn=-1 if label=='M' else 1
                assert np.min(sgn*(aa['positions_nm'][ix,2]-2000.))>0
                edge_mask=(regions[aa['ids'][:,0]]==label)
                rss=np.einsum('si,ij,sj->s',aa['burgers'][edge_mask],sigma,aa['planes'][edge_mask])
                expected=next(x['selected_slip']['rss_pa'] for x in sources if x['region']==label)
                np.testing.assert_allclose(rss,expected,rtol=1e-12,atol=1e-5)
            np.testing.assert_array_equal(sim.state['applied_stress'],sigma.ravel()[[0,4,8,5,2,1]])
            checks.append(dict(test=population+'_actual_native_geometry_load',status='PASS',nodes=len(aa['tags'])))
        combined=fields['MT'];tags=combined['tags'][:,1];cross=[]
        for population in ('M','T'):
            single=fields[population];ix=np.array([np.flatnonzero(tags==t)[0] for t in single['tags'][:,1]])
            for key in ('positions_nm','tags','constraints'):np.testing.assert_array_equal(combined[key][ix],single[key])
            delta=forces['MT'][ix]-forces[population];cross.append(delta)
            assert np.linalg.norm(delta)>1e-18
            checks.append(dict(test=population+'_coupling_from_other_side',status='PASS',cross_force_norm_N=float(np.linalg.norm(delta)),
                relative_to_single=float(np.linalg.norm(delta)/np.linalg.norm(forces[population]))))
        balance=float(np.linalg.norm(sum((x.sum(axis=0) for x in cross),np.zeros(3)))/sum(np.linalg.norm(x,axis=1).sum() for x in cross))
        assert balance<1e-8
        # Same actual MT configuration, zero applied stress: PK difference only.
        sim=networks['MT'];zero=Simulation(config(**{**c,'stress_pa':np.zeros((3,3)).tolist()}),a.repo,initial=sim.snapshot());fz,_=zero.evaluate()
        aa=fields['MT'];ij=aa['ids'];dr=(aa['positions_nm'][ij[:,1]]-aa['positions_nm'][ij[:,0]])*1e-9
        ext=np.zeros_like(forces['MT']);segment=.5*np.cross((sigma@(aa['burgers']*c['b_m']).T).T,dr)
        for end in (0,1):np.add.at(ext,ij[:,end],segment)
        pk_error=float(np.linalg.norm(forces['MT']-fz*sim.b**2-ext)/np.linalg.norm(ext));assert pk_error<1e-10
        checks.append(dict(test='common_macro_stress_native_PK_tensor',status='PASS',relative_error=pk_error))
        checks.append(dict(test='cross_network_pair_force_balance',status='PASS',normalized_residual=balance))
        save_json(a.out/'initial_material_metadata.json',metadata)
        save_json(a.out/'result.json',dict(status='PASS',scope='G1 initial geometry/common macro stress/one-network native elastic coupling only',
            checks=checks,common_interface=dict(id='CTB0',z_nm=2000.),orientation=dict(M=Q,T=ROT@Q),
            shared_contact_dynamics='NOT_RUN',host_identity_remesh_collision_topology='NOT_RUN',G1_gate='NOT_RUN',
            VCoNi='NOT_QUALIFIED',native_module=py.__file__,command=sys.argv))
        print('G1_INITIAL_FORCE_PROBE_PASS',len(checks),'G1_GATE_NOT_RUN')
    except Exception as exc:save_json(a.out/'result.json',dict(status='FAIL',checks=checks,error=repr(exc),G1_gate='NOT_RUN'));raise
    finally:py.finalize()


if __name__=='__main__':main()
