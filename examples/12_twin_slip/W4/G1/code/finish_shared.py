"""Bounded final review after all frozen native jobs, not a PASS shortcut."""
import argparse,json,subprocess,sys,time,shutil,os,datetime
from pathlib import Path
import numpy as np


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--historical-root',type=Path,required=True);p.add_argument('--g0-root',type=Path,required=True);p.add_argument('--repo',required=True);a=p.parse_args();root=a.root;code=root/'code';start=time.monotonic()
    jobs=json.loads((root/'campaign_manifest.json').read_text())['jobs']
    while True:
        complete=0
        for j in jobs:
            path=root/'logs'/(j['case_id']+'_receipt.json')
            if path.exists():
                r=json.loads(path.read_text())
                if r['exit_code']!=0:raise RuntimeError('Native matrix stopped; dependent review refused: '+j['case_id'])
                complete+=1
        if complete==len(jobs):break
        if time.monotonic()-start>3700:raise RuntimeError('CENSORED: bounded review wait; no gate inferred')
        print('REVIEW_WAIT_NATIVE',complete,'OF',len(jobs),flush=True);time.sleep(20)
    def execute(label,cmd,timeout=180,env=None):
        t=time.monotonic()
        with (root/(label+'.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=timeout,env=env)
        (root/(label+'_receipt.json')).write_text(json.dumps(dict(command=cmd,exit_code=r.returncode,wall_s=time.monotonic()-t),indent=2))
        print(label,'EXIT',r.returncode,flush=True)
        if r.returncode:raise RuntimeError('Review command failed; evidence retained: '+label)
    env=dict(os.environ,PYTHONPATH=str(code)+os.pathsep+str(a.g0_root/'code'),OMP_NUM_THREADS='2',OMP_PROC_BIND='false',OPENBLAS_NUM_THREADS='1')
    execute('identity_unit_tests',[sys.executable,'-m','unittest','discover','-s',str(code),'-p','test_material_identity.py','-v'],env=env)
    # Small byte-for-byte copy of real G0 evidence. No regenerated PASS JSON.
    old=a.historical_root/'g0_regression';dest=root/'g0_regression';dest.mkdir(exist_ok=False)
    for path in old.rglob('*'):
        if not path.is_file() or '__pycache__' in path.parts:continue
        if path.suffix in ('.json','.log','.txt','.py') or path.name in ('initial.npz','final.npz','initial.data','final.data'):
            target=dest/path.relative_to(old);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target)
            assert path.read_bytes()==target.read_bytes()
    (dest/'REFERENCE.json').write_text(json.dumps(dict(actual_run=str(old),same_native_build=True,source='Fresh actual A0, 4 trajectories and raw logs copied exactly. Full velocity/frames remain in historical root and its full archive.',auxiliary_strict_replay='FAIL retained'),indent=2))
    identity=[]
    for pop in ('M','T','MT'):
        oldpath=a.historical_root/'cases'/('base_'+pop)/'initial.npz';newpath=root/'cases'/('base_'+pop)/'initial.npz'
        with np.load(oldpath) as x,np.load(newpath) as y:
            for key in ('positions_nm','tags','ids','burgers','planes','cell_nm'):np.testing.assert_array_equal(x[key],y[key])
        identity.append(dict(population=pop,status='PASS',comparison='Exact physical initial graph versus old uncentered model',old=str(oldpath),new=str(newpath)))
    (root/'physical_initial_identity.json').write_text(json.dumps(identity,indent=2))
    repo=Path(a.repo);git=lambda p:subprocess.check_output(['git','-C',str(p),'rev-parse','HEAD'],text=True).strip()
    (root/'execution_environment.json').write_text(json.dumps(dict(client_date='2026-10-30',observed_vm_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        python=sys.executable,python_version=sys.version,numpy=np.__version__,vm_repository_revision=git(repo),native_revision=git(repo/'core/exadis'),
        native_semantic_patch=subprocess.check_output(['git','-C',str(repo/'core/exadis'),'diff','--ignore-cr-at-eol','--numstat'],text=True),
        Windows_source_branch='codex/w4-g1-shared-ctb',source_milestone='4048896 plus per-case exact code_snapshot and review tool snapshots',threads=2,
        force_evaluation='Native direct force on uniformly translated graph; physical graph unchanged; no force balancing correction'),indent=2))
    execute('force_audit',[sys.executable,str(code/'force_audit.py'),'--root',str(root),'--g0-root',str(a.g0_root),'--repo',a.repo,'--out',str(root/'force_audit')])
    execute('analysis_r1',[sys.executable,str(code/'analyze_shared.py'),'--root',str(root),'--out',str(root/'analysis_r1')],timeout=300)
    execute('report_r1',[sys.executable,str(code/'report_shared.py'),'--root',str(root),'--analysis',str(root/'analysis_r1'),'--out',str(root/'review_r1')],timeout=180)
    print('FINAL_G1_GATE',json.loads((root/'analysis_r1/gate_review.json').read_text())['status'],flush=True)


if __name__=='__main__':main()
