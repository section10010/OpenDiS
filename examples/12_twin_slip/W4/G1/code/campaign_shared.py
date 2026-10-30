"""Frozen bounded G1 matched-population campaign; no thickness ensemble."""
import argparse,copy,json,sys,time,os,subprocess
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--g0-root',type=Path,required=True)
    p.add_argument('--repo',required=True);p.add_argument('--spec',type=Path,required=True);p.add_argument('--prepare-only',action='store_true');a=p.parse_args()
    root=a.root;root.mkdir(exist_ok=True,parents=True);manifest=root/'campaign_manifest.json';code=Path(__file__).resolve().parent
    if not manifest.exists():
        spec=json.loads(a.spec.read_text());assert spec['schema']=='W4.G1.SharedEvolution.v1'
        assert json.loads((root/'admission_r1/result.json').read_text())['status']=='PASS'
        assert spec['initial_geometry_updates']=={'maxseg_nm':10.}
        jobs=[];(root/'case_specs').mkdir()
        levels=[('base',10.,.125),('timehalf',10.,.0625),('timequarter',10.,.03125),
                ('mesh5',5.,.03125),('mesh2p5',2.5,.03125),('mesh2p5_timehalf',2.5,.015625)]
        for level,mesh,step in levels:
            for population in ('M','T','MT'):
                c=copy.deepcopy(spec);c['common_native_updates'].update(maxseg_nm=mesh,minseg_nm=min(2.,mesh/4),maxdt_ps=step,max_node_step_nm=step)
                jobs.append(dict(case_id=level+'_'+population,population=population,group='forward',spec=c))
        for population in ('M','T','MT'):
            fine=copy.deepcopy(next(j['spec'] for j in jobs if j['case_id']=='mesh2p5_timehalf_'+population))
            no_wall=copy.deepcopy(fine);no_wall['contact_enabled']=False
            jobs.append(dict(case_id='no_ctb_'+population,population=population,group='no_ctb',spec=no_wall))
            reverse=copy.deepcopy(fine);reverse['common_native_updates']['stress_pa']=[[-v for v in row] for row in fine['common_native_updates']['stress_pa']]
            reverse['common_native_updates']['end_ns']=.01
            jobs.append(dict(case_id='reverse_'+population,population=population,group='reverse',spec=reverse,initial_case='mesh2p5_timehalf_'+population))
        # Operator ablations change no loading, source or contact parameter.
        for operator in ('topology','remesh'):
            c=copy.deepcopy(next(j['spec'] for j in jobs if j['case_id']=='timequarter_MT'))
            c['common_native_updates'][operator]=False
            jobs.append(dict(case_id=operator+'_off_MT',population='MT',group='operator_ablation',spec=c))
        records=[]
        for j in jobs:
            name=j['case_id'];(root/'case_specs'/(name+'.json')).write_text(json.dumps(j.pop('spec'),indent=2));records.append(j)
        manifest.write_text(json.dumps(dict(status='FROZEN_BEFORE_CAMPAIGN',jobs=records,max_wall_s=spec['criteria']['campaign_max_wall_s'],
            initial_geometry='Exactly common 10 nm polygon seed for every resolution, phase labels preserved across ablations',
            causal_admission='Paired signed area difference agrees at final levels and exceeds twice its last refinement change',
            observables='per-origin unsigned integration sweep, signed plastic tensor, stored/contact line, first contact and release',
            no_seed_CI=True),indent=2))
    if a.prepare_only:print('G1_PREPARED');return
    jobs=json.loads(manifest.read_text())['jobs'];(root/'cases').mkdir(exist_ok=True);(root/'logs').mkdir(exist_ok=True)
    start=time.monotonic();os.environ.update(OMP_NUM_THREADS='2',OMP_PROC_BIND='false',OPENBLAS_NUM_THREADS='1')
    for j in jobs:
        name=j['case_id'];receipt=root/'logs'/(name+'_receipt.json')
        if receipt.exists():
            prior=json.loads(receipt.read_text());assert prior['exit_code']==0,'Prior failed case retained; use a new case ID after repair'
            continue
        if time.monotonic()-start>json.loads(manifest.read_text())['max_wall_s']:
            (root/'campaign_stop.json').write_text(json.dumps(dict(status='CENSORED',next_case=name)));break
        cmd=[sys.executable,'-u',str(code/'run_shared.py'),'--g0-root',str(a.g0_root),'--repo',a.repo,'--spec',str(root/'case_specs'/(name+'.json')),
             '--population',j['population'],'--out',str(root/'cases'/name)]
        if 'initial_case' in j:
            parent=root/'cases'/j['initial_case'];assert json.loads((parent/'result.json').read_text())['status']=='PASS'
            cmd+=['--initial',str(parent/'final.npz'),'--sources',str(parent/'sources.json')]
        print('RUN',name,flush=True);t=time.monotonic();timed_out=False
        with (root/'logs'/(name+'.log')).open('w') as f:
            try:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=json.loads((root/'case_specs'/(name+'.json')).read_text())['common_native_updates']['max_wall_s']+60);exit_code=r.returncode
            except subprocess.TimeoutExpired:exit_code=None;timed_out=True
        receipt.write_text(json.dumps(dict(command=cmd,exit_code=exit_code,timed_out=timed_out,wall_s=time.monotonic()-t),indent=2))
        if exit_code!=0:
            result_path=root/'cases'/name/'result.json'
            status=json.loads(result_path.read_text()).get('status','FAIL') if result_path.exists() else ('CENSORED' if timed_out else 'FAIL')
            (root/'campaign_stop.json').write_text(json.dumps(dict(status=status,case=name)));raise RuntimeError('Case stopped: '+name+' '+status)
        print('DONE',name,round(time.monotonic()-t,2),flush=True)
    print('G1_CAMPAIGN_RECORDED',len(list((root/'cases').glob('*/result.json'))),'OF',len(jobs),flush=True)


if __name__=='__main__':main()
