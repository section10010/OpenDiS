"""Bounded native pair-kernel diagnosis; no change to production dynamics."""
import argparse,json,os,copy
from pathlib import Path
import numpy as np
from run_shared import bootstrap,seed_from_npz
from material_identity import MaterialIdentity


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--g0-root',type=Path,required=True);p.add_argument('--repo',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();bootstrap(a.g0_root)
    from native import Simulation,imports
    from configuration import config,save_json
    if a.out.exists():raise RuntimeError('Refuse overwrite')
    a.out.mkdir(parents=True);parent=a.root/'cases/mesh2p5_MT_budget_r2';c=config(**json.loads((parent/'spec.json').read_text())['common_native_updates']);sources=json.loads((parent/'sources.json').read_text())
    os.environ.update(OMP_NUM_THREADS='2',OMP_PROC_BIND='false');py,_=imports(a.repo);py.initialize(num_threads=2)
    try:
        seed=seed_from_npz(parent/'final.npz',c);s=Simulation(c,a.repo,initial=seed);meta=MaterialIdentity(sources,dict(id='CTB0',z_nm=2000.)).inspect(s.arrays(False))
        pairs=[[int(i),int(j)] for i in np.flatnonzero(meta['segment_host']=='M') for j in np.flatnonzero(meta['segment_host']=='T')]
        checks=[];reference=None
        for label,shift,reverse in [('original',[0,0,0],False),('swapped_pair_orientation',[0,0,0],True),('translated_origin',[-2000,-2000,-2000],False),('translated_and_swapped',[-2000,-2000,-2000],True)]:
            data=copy.deepcopy(seed);data['nodes']['positions']+=np.array(shift)/s.bnm;net=Simulation(c,a.repo,initial=data)
            listing=[[j,i] for i,j in pairs] if reverse else pairs
            f=np.array(py.compute_force_segseglist(net.g.net,mu=c['mu_pa'],nu=c['nu'],a=c['core_a_b'],segseglist=listing))*s.b**2
            if reference is None:reference=f.copy()
            error=float(np.linalg.norm(f-reference)/np.linalg.norm(reference));res=float(np.linalg.norm(f.sum(axis=0))/np.linalg.norm(f,axis=1).sum())
            np.savez_compressed(a.out/(label+'.npz'),force_N=f,tags=data['nodes']['tags'],positions_nm=data['nodes']['positions']*s.bnm)
            checks.append(dict(test=label,translation_nm=shift,relative_field_difference=error,pair_balance_residual=res))
        save_json(a.out/'result.json',dict(status='FAIL' if any(c['pair_balance_residual']>1e-8 for c in checks) else 'PASS',checks=checks,
            interpretation='Diagnostic only. Establishes whether subtraction, origin or ordered-pair convention alone removes the failed balance. No balancing projection or production force modification. Independent high-precision reference NOT_RUN.'))
        print('NATIVE_FORCE_INVARIANCE_DIAGNOSIS',checks)
    finally:py.finalize()


if __name__=='__main__':main()
