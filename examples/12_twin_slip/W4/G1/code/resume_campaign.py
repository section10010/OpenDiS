"""One frozen resource-only extension after measured G1 CENSORED evidence.

Never resumes the unqualified cumulative ledger. Each aliased case starts again
at the exact original initial graph, under a new immutable case ID.
"""
import argparse,json,sys,time,subprocess,os,copy
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--g0-root',type=Path,required=True);p.add_argument('--repo',required=True)
    a=p.parse_args();root=a.root;manifest=json.loads((root/'campaign_manifest.json').read_text());extension=root/'resource_extension_r2.json'
    if not extension.exists():
        evidence=json.loads((root/'cases/mesh2p5_MT/result.json').read_text());assert evidence['status']=='CENSORED'
        pending=[];aliases={}
        for j in manifest['jobs']:
            path=root/'cases'/j['case_id']/'result.json'
            if path.exists() and json.loads(path.read_text())['status']=='PASS':continue
            if path.exists():assert json.loads(path.read_text())['status']=='CENSORED','Physical/numerical FAIL requires repair, not a budget extension'
            record=copy.deepcopy(j);record['actual_case_id']=j['case_id']+'_budget_r2';pending.append(record);aliases[j['case_id']]=record['actual_case_id']
            spec=json.loads((root/'case_specs'/(j['case_id']+'.json')).read_text());updated=copy.deepcopy(spec);updated['common_native_updates']['max_wall_s']=750
            difference=[k for k in spec['common_native_updates'] if spec['common_native_updates'][k]!=updated['common_native_updates'][k]];assert difference==['max_wall_s']
            (root/'case_specs'/(record['actual_case_id']+'.json')).write_text(json.dumps(updated,indent=2))
        extension.write_text(json.dumps(dict(status='FROZEN_BEFORE_RETRY',maximum_case_wall_s=750,total_extension_wall_s=3600,jobs=pending,
            basis=dict(censored_case='mesh2p5_MT',wall_s=evidence['wall_s'],physical_time_ns=evidence['time_ns'],target_ns=.025,
                linear_estimate_s=evidence['wall_s']*.025/evidence['time_ns'],extra_time_refinement_factor=2),
            allowed_change='max_wall_s only; all science criteria and physical/numerical parameters unchanged',
            restart='New complete runs, never concatenate CENSORED trajectories'),indent=2))
        (root/'case_aliases.json').write_text(json.dumps(dict(aliases=aliases,historical_censored='cases/mesh2p5_MT/result.json',
            original_orchestrator_stop='campaign_stop.json classified nonzero process exit as FAIL; actual native result is CENSORED. Fixed orchestration preserves native status.'),indent=2))
    ext=json.loads(extension.read_text());aliases=json.loads((root/'case_aliases.json').read_text())['aliases'];start=time.monotonic()
    os.environ.update(OMP_NUM_THREADS='2',OMP_PROC_BIND='false',OPENBLAS_NUM_THREADS='1')
    for j in ext['jobs']:
        name=j['actual_case_id'];receipt=root/'logs'/(name+'_receipt.json')
        if receipt.exists():assert json.loads(receipt.read_text())['exit_code']==0;continue
        if time.monotonic()-start>ext['total_extension_wall_s']:
            (root/'resource_extension_stop.json').write_text(json.dumps(dict(status='CENSORED',next_case=name)));return
        cmd=[sys.executable,'-u',str(Path(__file__).resolve().parent/'run_shared.py'),'--g0-root',str(a.g0_root),'--repo',a.repo,
            '--spec',str(root/'case_specs'/(name+'.json')),'--population',j['population'],'--out',str(root/'cases'/name)]
        if 'initial_case' in j:
            parent=root/'cases'/aliases.get(j['initial_case'],j['initial_case']);assert json.loads((parent/'result.json').read_text())['status']=='PASS'
            cmd+=['--initial',str(parent/'final.npz'),'--sources',str(parent/'sources.json')]
        print('RUN',name,flush=True);t=time.monotonic();timed_out=False
        with (root/'logs'/(name+'.log')).open('w') as f:
            try:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=810);code=r.returncode
            except subprocess.TimeoutExpired:code=None;timed_out=True
        receipt.write_text(json.dumps(dict(command=cmd,exit_code=code,timed_out=timed_out,wall_s=time.monotonic()-t),indent=2))
        if code!=0:
            result=root/'cases'/name/'result.json';status=json.loads(result.read_text()).get('status','FAIL') if result.exists() else 'CENSORED' if timed_out else 'FAIL'
            (root/'resource_extension_stop.json').write_text(json.dumps(dict(status=status,case=name)));raise RuntimeError('Extended case stopped: '+name)
        print('DONE',name,round(time.monotonic()-t,2),flush=True)
    print('G1_RESOURCE_EXTENSION_COMPLETE',len(ext['jobs']),flush=True)


if __name__=='__main__':main()
