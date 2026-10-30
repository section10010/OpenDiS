"""One bounded 1.25 nm refinement after preserved 5->2.5 nm tail failure."""
import argparse,json,sys,time,subprocess,os
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--g0-root',type=Path,required=True);p.add_argument('--repo',required=True);a=p.parse_args();root=a.root;manifest=root/'spatial_extension_manifest.json'
    if not manifest.exists():
        preflight=json.loads((root/'mesh1p25_force_preflight.json').read_text());assert preflight['status']=='PASS'
        jobs=[]
        for pop in ('M','T','MT'):
            base=json.loads((root/'case_specs'/('mesh2p5_timehalf_'+pop+'.json')).read_text())
            base['common_native_updates'].update(maxseg_nm=1.25,minseg_nm=.3125,max_wall_s=1500)
            name='mesh1p25_'+pop;(root/'case_specs'/(name+'.json')).write_text(json.dumps(base,indent=2));jobs.append(dict(case_id=name,population=pop))
        manifest.write_text(json.dumps(dict(status='FROZEN_BEFORE_NATIVE_EXTENSION',jobs=jobs,max_case_wall_s=1500,max_total_wall_s=3000,
            reason='Actual MT T-tail 5->2.5 nm change=0.06336789758 exceeds original 0.02. Temporal refinement passed. Preserve failure; independent 1.25 spatial refinement only.',
            unchanged='Original polygon seed, material, stress, contact, mobility, time step and node displacement limit',
            criterion='Original 2 percent primary / 5 percent event. Last spatial pair is 2.5->1.25; old 5->2.5 remains UNRESOLVED. No further automatic refinement.'),indent=2))
    start=time.monotonic();os.environ.update(OMP_NUM_THREADS='2',OMP_PROC_BIND='false',OPENBLAS_NUM_THREADS='1')
    for j in json.loads(manifest.read_text())['jobs']:
        name=j['case_id'];receipt=root/'logs'/(name+'_receipt.json')
        if receipt.exists():assert json.loads(receipt.read_text())['exit_code']==0;continue
        if j['population']=='MT':
            # Do not oversubscribe the most costly production release case.
            while not (root/'logs/remesh_off_MT_receipt.json').exists():
                if time.monotonic()-start>3000:raise RuntimeError('CENSORED: bounded resource wait')
                print('WAIT_MAIN_MATRIX_BEFORE_FINE_MT',flush=True);time.sleep(20)
        remaining=3000-(time.monotonic()-start)
        if remaining<1500:raise RuntimeError('CENSORED: insufficient frozen remaining resource allowance')
        cmd=[sys.executable,str(Path(__file__).resolve().parent/'run_shared.py'),'--g0-root',str(a.g0_root),'--repo',a.repo,'--spec',str(root/'case_specs'/(name+'.json')),'--population',j['population'],'--out',str(root/'cases'/name)];t=time.monotonic()
        print('RUN',name,flush=True)
        with (root/'logs'/(name+'.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=1560)
        receipt.write_text(json.dumps(dict(command=cmd,exit_code=r.returncode,wall_s=time.monotonic()-t),indent=2));print('DONE',name,r.returncode,round(time.monotonic()-t,2),flush=True)
        if r.returncode:raise RuntimeError('Fine native case stopped: '+name)
    (root/'spatial_extension_result.json').write_text(json.dumps(dict(status='PASS',scope='Three native numerical runs completed; physical/convergence decision belongs to independent analyze_shared',jobs=json.loads(manifest.read_text())['jobs']),indent=2));print('SPATIAL_EXTENSION_NATIVE_COMPLETE',flush=True)


if __name__=='__main__':main()
