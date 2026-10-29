"""Native GLIDE/contact constitutive tests; prescribed forces are unit tests only."""
import argparse,os
from pathlib import Path
import numpy as np
from configuration import config,save_json
from geometry import pure_shear
from native import Simulation,imports
from contact import Contact


def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True);os.environ['OMP_NUM_THREADS']='2';py,_=imports(a.repo);py.initialize(num_threads=2);checks=[]
    try:
        for mode in ('HI','HII'):
            c=config(kind='propagation_loop',mode=mode,integrator='EulerForward',vmax_m_s=None,stress_pa=pure_shear('matrix',mode,200e6).tolist())
            sim=Simulation(c,a.repo);raw=sim.arrays();n=np.array(sim.sources[0]['n']);q=np.array(sim.sources[0]['q']);pvec=np.array(sim.sources[0]['p'])
            for which,sgn,idx in [('lower',1,np.argmin(raw['positions_nm'][:,2])),('upper',-1,np.argmax(raw['positions_nm'][:,2]))]:
                wall=dict(interface_id=which,z_nm=float(raw['positions_nm'][idx,2]),inward=sgn,enabled=True)
                cc={**c,'walls':[wall]};sim.contact=Contact(cc);tag=raw['tags'][idx]
                into=-sgn*np.sign(pvec[2])*pvec+.3*q;away=-into+.6*q
                vin=sim.glide.OneNodeMobility(sim.net,sim.state,tag,into,False)
                vc=sim.contact.project_one(sim,sim.net,sim.state,tag,vin)
                vt=sim.glide.OneNodeMobility(sim.net,sim.state,tag,q,False)
                vtc=sim.contact.project_one(sim,sim.net,sim.state,tag,vt)
                vo=sim.glide.OneNodeMobility(sim.net,sim.state,tag,away,False)
                voc=sim.contact.project_one(sim,sim.net,sim.state,tag,vo)
                scale=max(np.linalg.norm(vin),1e-99)
                assert abs(vc[2])/scale<1e-12 and abs(vc@n)/scale<1e-12 and into@vc>=-1e-12*scale
                assert np.linalg.norm(vtc-vt)<1e-12*scale and np.linalg.norm(vtc)>0
                assert np.linalg.norm(voc-vo)<1e-12*scale and sgn*vo[2]>0
                dt=sim.contact.event_dt(sim,np.zeros_like(raw['positions_nm']));assert np.isinf(dt)
                all_v=np.array(sim.py.compute_mobility(sim.g.net,mobility=sim.glide.mobility,
                    nodeforces=np.tile(into,(len(raw['tags']),1)),nodetags=raw['tags']))
                indexed=sim.contact.project_all(sim,sim.net,sim.state,all_v)
                reference=np.array([sim.contact.project_one(sim,sim.net,sim.state,t,v) for t,v in zip(raw['tags'],all_v)])
                np.testing.assert_allclose(indexed,reference,rtol=1e-12,atol=1e-10)
                sim.contact.validate(sim,'unit_test')
                checks.append(dict(mode=mode,wall=which,blocked_normal_ratio=abs(vc[2])/scale,
                    retained_tangent_ratio=float(np.linalg.norm(vtc)/np.linalg.norm(vt)),release_ratio=float(np.linalg.norm(voc)/np.linalg.norm(vo)),status='PASS'))
            # Real force and native Euler hit the first contact without clipping.
            cc={**c,'walls':[dict(interface_id='lower',z_nm=1900.,inward=1,enabled=True),dict(interface_id='upper',z_nm=2100.,inward=-1,enabled=True)]}
            sim.contact=Contact(cc);_,v=sim.evaluate();event=sim.contact.event_dt(sim,v)
            assert np.isfinite(event) and event>0
            sim.advance(event);sim.contact.validate(sim,'native_event_step')
            assert min(abs(sim.arrays(False)['positions_nm'][:,2]-1900.).min(),abs(sim.arrays(False)['positions_nm'][:,2]-2100.).min())<1e-8
            sim.evaluate();sim.operations()
            checks.append(dict(mode=mode,test='native_event_hit_and_all_operations',event_dt_s=event,status='PASS'))
        save_json(a.out/'result.json',dict(status='PASS',checks=checks,scope='native unilateral glide-compatible mobility, tangent and release'))
        print('CONTACT_NATIVE_PASS',len(checks))
    except Exception as exc:save_json(a.out/'result.json',dict(status='FAIL',checks=checks,error=repr(exc)));raise
    finally:py.finalize()


if __name__=='__main__':main()
