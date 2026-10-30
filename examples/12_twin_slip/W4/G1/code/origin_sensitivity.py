"""Frozen alternative numerical-origin trajectories, same physical graphs."""
import argparse,json,copy,subprocess,sys,time,os
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--g0-root',type=Path,required=True);p.add_argument('--repo',required=True);p.add_argument('--prepare-only',action='store_true');a=p.parse_args();root=a.root;manifest=root/'origin_sensitivity_manifest.json'
    if not manifest.exists():
        jobs=[]
        for pop in ('M','T','MT'):
            spec=json.loads((root/'case_specs'/('mesh5_'+pop+'.json')).read_text());assert spec['force_evaluation_origin_nm']==[2000.,2000.,2000.]
            spec['force_evaluation_origin_nm']=[1900.,2000.,2000.];name='origin_shift100_'+pop
            (root/'case_specs'/(name+'.json')).write_text(json.dumps(spec,indent=2));jobs.append(dict(case_id=name,population=pop))
        manifest.write_text(json.dumps(dict(status='FROZEN_BEFORE_ORIGIN_TESTS',jobs=jobs,max_total_wall_s=900,
            criterion='Each paired MT-minus-single signed area retains sign and differs by no more than 2 percent of the original paired difference. No physical geometry is moved.'),indent=2))
    if a.prepare_only:print('ORIGIN_SENSITIVITY_FROZEN');return
    start=time.monotonic();os.environ.update(OMP_NUM_THREADS='2',OMP_PROC_BIND='false',OPENBLAS_NUM_THREADS='1');results={}
    for j in json.loads(manifest.read_text())['jobs']:
        name=j['case_id'];out=root/'cases'/name;receipt=root/'logs'/(name+'_receipt.json')
        if not receipt.exists():
            if time.monotonic()-start>900:raise RuntimeError('Origin sensitivity budget reached')
            cmd=[sys.executable,str(Path(__file__).resolve().parent/'run_shared.py'),'--g0-root',str(a.g0_root),'--repo',a.repo,'--spec',str(root/'case_specs'/(name+'.json')),'--population',j['population'],'--out',str(out)];t=time.monotonic()
            with (root/'logs'/(name+'.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=810)
            receipt.write_text(json.dumps(dict(command=cmd,exit_code=r.returncode,wall_s=time.monotonic()-t),indent=2));assert r.returncode==0
        results[j['population']]=json.loads((out/'result.json').read_text());assert results[j['population']]['status']=='PASS'
    checks=[]
    for host in ('M','T'):
        old_mt=json.loads((root/'cases/mesh5_MT/result.json').read_text())['swept_area_nm2'][host]
        old_one=json.loads((root/'cases'/('mesh5_'+host)/'result.json').read_text())['swept_area_nm2'][host]
        base=old_mt-old_one;new=results['MT']['swept_area_nm2'][host]-results[host]['swept_area_nm2'][host];error=abs(new-base)/abs(base)
        checks.append(dict(host=host,base_difference_nm2=base,shifted_origin_difference_nm2=new,relative_difference=error,status='PASS' if error<=.02 and base*new>0 else 'FAIL'))
    (root/'origin_sensitivity_result.json').write_text(json.dumps(dict(status='PASS' if all(c['status']=='PASS' for c in checks) else 'FAIL',checks=checks,physical_geometry='UNCHANGED'),indent=2));print('ORIGIN_SENSITIVITY',checks)


if __name__=='__main__':main()
