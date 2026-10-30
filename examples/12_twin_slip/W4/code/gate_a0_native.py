"""A0 native numerical response checks; this is not an A1 source qualification."""
import argparse,os
from pathlib import Path
import numpy as np
from configuration import config,save_json
from geometry import pure_shear,build
from native import imports,Simulation


def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True);os.environ['OMP_NUM_THREADS']='2';os.environ['OMP_PROC_BIND']='false'
    py,exa=imports(a.repo);py.initialize(num_threads=2);checks=[]
    try:
        c=config(gate='A0',stress_pa=pure_shear('matrix','HI',100e6).tolist(),vmax_m_s=None,integrator='EulerForward')
        s=Simulation(c,a.repo);f,v=s.evaluate();f=f.copy();v=v.copy()
        double=Simulation(config(**{**c,'Medge':2*c['Medge'],'Mscrew':2*c['Mscrew']}),a.repo)
        f2,v2=double.evaluate();np.testing.assert_allclose(f2,f,rtol=1e-12,atol=1e-3);np.testing.assert_allclose(v2,2*v,rtol=1e-11,atol=1e-2)
        checks.append(dict(test='native_drag_parameter_response',velocity_norm_ratio=float(np.linalg.norm(v2)/np.linalg.norm(v)),status='PASS'))
        for key in ('Medge','Mscrew'):
            changed=Simulation(config(**{**c,key:2*c[key]}),a.repo);_,vv=changed.evaluate()
            assert np.linalg.norm(vv-v)>1e-5*np.linalg.norm(v)
            checks.append(dict(test='native_individual_'+key+'_response',velocity_norm_ratio=float(np.linalg.norm(vv)/np.linalg.norm(v)),status='PASS'))
        for changes in ({'source_gauge_nm':360.},{'backbone_offset_nm':160.},{'source_count':2}):
            cc=config(**{**c,**changes});model=Simulation(cc,a.repo);aa=model.arrays();src=model.sources[0]
            pts=aa['positions_nm'][src['initial_node_indices']]-np.array(src['center_nm'])
            np.testing.assert_allclose(np.ptp(pts@src['q']),cc['source_gauge_nm'],rtol=1e-12)
            np.testing.assert_allclose(np.ptp(pts@src['n']),cc['backbone_offset_nm'],rtol=1e-12)
            assert len(model.sources)==cc['source_count']
            checks.append(dict(test='native_geometry_readback',changed=changes,actual_nodes=model.g.num_nodes(),status='PASS'))
        s.advance(1e-14);double.advance(1e-14)
        initial=build(c)[0][:,:3];dx=s.g.get_positions()*s.bnm-initial;dx2=double.g.get_positions()*double.bnm-initial
        np.testing.assert_allclose(dx2,2*dx,rtol=1e-8,atol=1e-10)
        checks.append(dict(test='native_euler_actual_displacement_response',max_motion_nm=float(np.linalg.norm(dx,axis=1).max()),status='PASS'))
        forces=[]
        seed=Simulation(c,a.repo).snapshot()
        for factor in (0,1,2):
            t=Simulation(config(**{**c,'stress_pa':(np.array(c['stress_pa'])*factor).tolist()}),a.repo,initial=seed)
            ff,_=t.evaluate();forces.append(ff.copy())
        np.testing.assert_allclose(forces[2]-forces[1],forces[1]-forces[0],rtol=1e-9,atol=1e-2)
        assert np.linalg.norm(forces[1]-forces[0])>0
        checks.append(dict(test='native_PK_stress_tensor_response',status='PASS'))
        bulk=Simulation(c,a.repo);bf,bv=bulk.evaluate();idx=[i for i,k in enumerate(bulk.g.get_nodes_data()['constraints'].ravel()) if k==0][:6]
        errors=[]
        for i in idx:
            one=bulk.force.OneNodeForce(bulk.net,bulk.state,bulk.g.get_tags()[i],False,match_global=True)
            errors.append(float(np.linalg.norm(one-bf[i])/max(np.linalg.norm(bf[i]),1.)))
        assert max(errors)<1e-8
        checks.append(dict(test='native_direct_one_node_vs_global_force',max_relative_error=max(errors),status='PASS'))
        from test_native_pins_w1 import PinMerge
        pin=PinMerge();pin.test_reject_move_fixed_pin();pin.test_pin_tag_survives_when_second_argument();pin.test_distinct_pins_cannot_merge()
        checks.append(dict(test='W1_native_pin_contract',cases=3,status='PASS'))
        save_json(a.out/'config.json',c);np.savez_compressed(a.out/'native_response.npz',initial_positions_nm=initial,displacement_nm=dx,double_drag_displacement_nm=dx2,force_N=f*c['b_m']**2,velocity_m_s=v*c['b_m'])
        save_json(a.out/'result.json',dict(status='PASS',gate='A0_native_response',checks=checks,scope='Real native force, mobility and one-step numerical checks; not source/CTB qualification'))
        print('A0_NATIVE_PASS',len(checks),'checks')
    except Exception as exc:
        save_json(a.out/'result.json',dict(status='FAIL',gate='A0_native_response',checks=checks,error=repr(exc)));raise
    finally:py.finalize()


if __name__=='__main__':main()
