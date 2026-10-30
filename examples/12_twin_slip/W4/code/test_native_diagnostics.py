"""Real native force/mobility/reaction validation before G0 dynamics."""
import argparse,os
from pathlib import Path
import numpy as np
from configuration import config,save_json
from geometry import pure_shear
from native import Simulation,imports
from contact import Contact


def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise RuntimeError('Refuse existing output')
    a.out.mkdir(parents=True);os.environ['OMP_NUM_THREADS']='2';py,_=imports(a.repo);py.initialize(num_threads=2);checks=[]
    try:
        for mode in ('HI','HII'):
            c=config(kind='propagation_loop',mode=mode,stress_pa=pure_shear('matrix',mode,300e6).tolist(),integrator='EulerForward',vmax_m_s=100.)
            sim=Simulation(c,a.repo);z=sim.arrays(False)['positions_nm'][:,2].min()
            sim.c['walls']=[dict(interface_id='lower',z_nm=float(z),inward=1,enabled=True)];sim.contact=Contact(sim.c)
            f,v=sim.evaluate();d={k:x.copy() for k,x in sim.velocity_diagnostic.items()};tags=sim.g.get_tags()
            native_raw=np.array(py.compute_mobility(sim.g.net,mobility=sim.glide.mobility,nodeforces=f,nodetags=tags))*sim.b
            np.testing.assert_allclose(native_raw,d['raw_velocity_m_s'],rtol=1e-12,atol=1e-9)
            reconstructed=np.array(py.compute_mobility(sim.g.net,mobility=sim.glide.mobility,nodeforces=f+d['reaction_N']/sim.b**2,nodetags=tags))*sim.b
            np.testing.assert_allclose(reconstructed,d['constrained_velocity_m_s'],rtol=1e-10,atol=1e-7)
            assert np.linalg.norm(d['reaction_N'])>0 and (d['multipliers_N']>=0).all()
            assert np.linalg.norm(d['capped_velocity_m_s'],axis=1).max()<=100*(1+1e-12)
            assert np.max(np.linalg.norm(d['raw_velocity_m_s'],axis=1))>100
            sim.c['vmax_m_s']=None;sim.evaluate();u=sim.velocity_diagnostic
            np.testing.assert_allclose(u['raw_velocity_m_s'],d['raw_velocity_m_s'],rtol=1e-12,atol=1e-9)
            np.testing.assert_allclose(u['capped_velocity_m_s'],u['constrained_velocity_m_s'],rtol=1e-12)
            checks.append(dict(mode=mode,status='PASS',reaction_reconstruction=True,raw_contact_capped_distinguished=True))
        save_json(a.out/'result.json',dict(status='PASS',checks=checks,native_engine=py.__file__))
        print('NATIVE_DIAGNOSTIC_PASS',len(checks))
    except Exception as exc:save_json(a.out/'result.json',dict(status='FAIL',checks=checks,error=repr(exc)));raise
    finally:py.finalize()


if __name__=='__main__':main()
