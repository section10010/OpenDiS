"""Fresh unchanged-core A0 and four actual G0 native trajectory replays."""
import argparse,json,sys,time,subprocess,os
from pathlib import Path
import numpy as np


def main():
    p=argparse.ArgumentParser();p.add_argument('--g0-root',type=Path,required=True);p.add_argument('--repo',required=True)
    p.add_argument('--build',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();g0=a.g0_root;out=a.out
    if out.exists():raise RuntimeError('Refuse overwrite')
    os.environ.update(OMP_NUM_THREADS='2',OMP_PROC_BIND='false',OPENBLAS_NUM_THREADS='1')
    cmd=[sys.executable,str(g0/'code/qualify_a0.py'),'--repo',a.repo,'--build',a.build,'--root',str(out)]
    # qualify_a0 requires a new empty root and creates it itself.
    parentlog=out.parent/'g0_regression_A0.log';t=time.monotonic()
    with parentlog.open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=900)
    (out/'fresh_a0_receipt.json').write_text(json.dumps(dict(command=cmd,exit_code=r.returncode,wall_s=time.monotonic()-t),indent=2))
    assert r.returncode==0
    checks=[]
    for name in ('configuration.py','geometry.py','native.py','run.py','contact.py','diagnostics.py'):
        assert (out/'code'/name).read_bytes()==(g0/'code'/name).read_bytes()
    def execute(label,cmd,cwd=None):
        t=time.monotonic()
        with (out/(label+'.log')).open('w') as f:r=subprocess.run(cmd,cwd=cwd,stdout=f,stderr=subprocess.STDOUT,timeout=240)
        (out/(label+'_receipt.json')).write_text(json.dumps(dict(command=cmd,exit_code=r.returncode,wall_s=time.monotonic()-t),indent=2));assert r.returncode==0,label
    for label,script in [('contact_admission','test_contact_native.py'),('diagnostic_admission','test_native_diagnostics.py')]:
        execute(label,[sys.executable,str(out/'code'/script),'--repo',a.repo,'--out',str(out/label)])
    execute('python_regressions',[sys.executable,'-m','unittest','test_configuration','test_lineage','test_diagnostics','-v'],out/'code')
    for family in ('loop_HI','threader_HII'):
        for wall in ('none','both'):
            name=f'G0_{family}_{wall}_uncapped_diagnostic';dest=out/'cases'/name;dest.parent.mkdir(exist_ok=True)
            execute(name,[sys.executable,str(out/'code/run.py'),'--repo',a.repo,'--config',str(g0/'configs'/(name+'.json')),'--out',str(dest),'--a0',str(out/'A0/gate_result.json')])
            old=json.loads((g0/'cases'/name/'result.json').read_text());new=json.loads((dest/'result.json').read_text());assert new['status']=='PASS'
            error=0.
            with np.load(g0/'cases'/name/'final.npz') as x,np.load(dest/'final.npz') as y:
                for key in ('positions_nm','tags','ids','burgers','planes','force_N','velocity_m_s'):
                    assert x[key].shape==y[key].shape,(name,key)
                    scale=max(float(np.linalg.norm(x[key])),1e-25);error=max(error,float(np.linalg.norm(x[key]-y[key]))/scale)
            assert error<1e-10,(name,error)
            checks.append(dict(case_id=name,status='PASS',final_array_relative_max=error,steps=new.get('steps'),native_result='PASS'))
    (out/'result.json').write_text(json.dumps(dict(status='PASS',A0='PASS',pin_native='PASS',contact_native='PASS',diagnostic_native='PASS',
        unchanged_production_source='BYTE_IDENTICAL',four_native_trajectory_replays=checks,scope='Preserves authoritative G0 80-case R5. Replays representative uncapped none/both HI/HII; not a duplicate 80-case campaign.'),indent=2))
    print('G0_REGRESSION_PASS',len(checks))


if __name__=='__main__':main()
