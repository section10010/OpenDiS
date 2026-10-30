"""Native replay after adding capability refusals; no equation changes."""
import argparse,json,os,time
from pathlib import Path
from run_shared import bootstrap,run_case
import numpy as np


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--g0-root',type=Path,required=True);p.add_argument('--repo',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();bootstrap(a.g0_root)
    from native import imports
    from configuration import save_json
    if a.out.exists():raise RuntimeError('Refuse overwrite')
    a.out.mkdir(parents=True);os.environ.update(OMP_NUM_THREADS='2',OMP_PROC_BIND='false');py,_=imports(a.repo);py.initialize(num_threads=2);checks=[]
    try:
        for pop in ('M','T','MT'):
            spec=json.loads((a.root/'case_specs'/('base_'+pop+'.json')).read_text());r=run_case(spec,pop,a.repo,a.out/pop)
            original=json.loads((a.root/'cases'/('base_'+pop)/'result.json').read_text());err={}
            for host in pop:err[host]=abs(r['swept_area_nm2'][host]-original['swept_area_nm2'][host])/original['swept_area_nm2'][host]
            assert max(err.values())<spec['criteria']['primary_refinement_relative']
            with np.load(a.out/pop/'initial.npz') as x,np.load(a.root/'cases'/('base_'+pop)/'initial.npz') as y:
                for k in ('positions_nm','tags','ids','burgers','planes'):np.testing.assert_array_equal(x[k],y[k])
            checks.append(dict(population=pop,status='PASS',same_initial_graph='EXACT',unsigned_sweep_relative=err))
        save_json(a.out/'result.json',dict(status='PASS',checks=checks,change='Only refusal checks added to run_case before building. Shared mobility, identity, equations and stage loop unchanged.'))
        print('G1_GUARD_ADMISSION_PASS',len(checks))
    finally:py.finalize()


if __name__=='__main__':main()
