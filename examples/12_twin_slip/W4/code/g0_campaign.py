"""Freeze and sequentially execute the actual 24-input G0 screen and refinements.

All physical settings come from preserved plan inputs. Extra cases are explicit
one-factor diagnostics, not an ensemble or four-thickness production campaign.
"""
import argparse,json,sys,time,shutil,subprocess,os
from pathlib import Path
import numpy as np
from configuration import config,save_json,load
from geometry import build
from record import run_logged

MASKS=('none','lower','upper','both')
FAMILIES=(('loop','HI'),('threader','HII'))


def prepare(plan,root):
    folder=root/'frozen_inputs';folder.mkdir();jobs=[]
    manifest=json.loads((plan/'configs/W4_G0_generated/manifest.json').read_text())
    assert len(manifest['jobs'])==24
    for entry in manifest['jobs']:
        src=plan/'configs/W4_G0_generated'/entry['config'];c=load(src)
        shutil.copy2(src,folder/entry['config']);jobs.append(dict(case_id=entry['case_id'],config=c,group='screen',**{k:entry[k] for k in ('kind','mode','wall')}))
    byname={j['case_id']:j for j in jobs}
    for kind,mode in FAMILIES:
        for mask in MASKS:
            seed=byname[f'G0_{kind}_{mode}_{mask}_uncapped_diagnostic']['config']
            extras=[('mesh40',dict(maxseg_nm=40.)),('mesh10',dict(maxseg_nm=10.)),
                ('mesh10_timehalf',dict(maxseg_nm=10.,maxdt_ps=.25,max_node_step_nm=.125)),
                ('mobility_half_rescaled',dict(Medge=seed['Medge']/2,Mscrew=seed['Mscrew']/2,end_ns=.2,maxdt_ps=1.))]
            if mask in ('lower','both'):
                extras += [('tol_tenth',dict(contact_tol_nm=1e-7)),('rann_half',dict(rann_nm=.3)),
                           ('topology_off',dict(topology=False)),('remesh_off',dict(remesh=False))]
            for label,changes in extras:
                name=f'G0_{kind}_{mode}_{mask}_{label}';c=config(**{**seed,**changes,'name':name})
                jobs.append(dict(case_id=name,config=c,group='refinement' if label not in ('topology_off','remesh_off') else 'operator_ablation',kind=kind,mode=mode,wall=mask))
    # Paired arrays are exact within each resolution/diagnostic, names are provenance.
    groups={}
    for j in jobs:
        tail=j['case_id'].split(j['wall']+'_',1)[1];key=(j['kind'],j['mode'],tail)
        arrays=build(j['config'])[:2]
        if key in groups:
            assert all(np.array_equal(a,b) for a,b in zip(groups[key],arrays)),key
        else:groups[key]=arrays
    save_json(root/'frozen_criteria.json',dict(status='FROZEN_BEFORE_NATIVE_SCREEN',source='configs/gate_criteria.json',
        primary_refinement_relative=.02,event_time_refinement_relative=.05,burgers_residual_b=1e-8,
        contact_power_normalized=1e-8,contact_normal_velocity_normalized=1e-8,penetration_nm=1e-6,
        plastic_closure='norm <= 1e-12 + 1e-8 * native_strain_norm',
        contact_normalization='force norm times max(raw, constrained velocity norm); near-zero outputs retain input scale',
        nonadditivity='UNRESOLVED unless greater than sum of all four absolute refinement differences; no statistical CI from frames',
        cap_sensitivity='Compare actual-time responses; a >2 percent change invalidates cap-insensitive timing. Second-wall causal sign must survive each resolved diagnostic; uncapped is mathematical only.',
        resource=dict(max_campaign_wall_s=3600,max_case_wall_s='preserve each input',max_nodes=5000,max_segments=10000),
        planned_cases=len(jobs),screen_cases=24,scope='effective FCC I0; not VCoNi calibration or collective thickness comparison'))
    (root/'configs').mkdir()
    for j in jobs:save_json(root/'configs'/(j['case_id']+'.json'),j['config'])
    save_json(root/'campaign_manifest.json',dict(status='NOT_RUN',jobs=[{k:v for k,v in j.items() if k!='config'} for j in jobs]))
    return jobs


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--plan',type=Path,required=True);p.add_argument('--repo',required=True);p.add_argument('--prepare-only',action='store_true');a=p.parse_args()
    root=a.root.resolve();code=Path(__file__).resolve().parent
    if not (root/'campaign_manifest.json').exists():jobs=prepare(a.plan,root)
    else:jobs=[dict(j,config=load(root/'configs'/(j['case_id']+'.json'))) for j in json.loads((root/'campaign_manifest.json').read_text())['jobs']]
    if a.prepare_only:print('PREPARED',len(jobs));return
    assert json.loads((root/'A0/gate_result.json').read_text())['status']=='PASS'
    for name in ('contact_admission','diagnostic_admission'):
        assert json.loads((root/name/'result.json').read_text())['status']=='PASS',name
    started=time.monotonic();records=[];os.environ['OMP_NUM_THREADS']='2'
    (root/'cases').mkdir(exist_ok=True);(root/'logs').mkdir(exist_ok=True)
    for j in jobs:
        name=j['case_id'];out=root/'cases'/name;receipt=root/'logs'/(name+'_receipt.json')
        if receipt.exists():
            records.append(json.loads(receipt.read_text()));continue
        if time.monotonic()-started>3600:
            save_json(root/'campaign_stop.json',dict(status='CENSORED',reason='frozen campaign budget',next_case=name));break
        cmd=[sys.executable,'-u',str(code/'run.py'),'--repo',a.repo,'--config',str(root/'configs'/(name+'.json')),'--out',str(out),'--a0',str(root/'A0/gate_result.json')]
        print('RUN',name,flush=True)
        r=run_logged(cmd,code,receipt,root/'logs'/(name+'.log'),j['config']['max_wall_s']+60);records.append(r)
        result=json.loads((out/'result.json').read_text()) if (out/'result.json').exists() else dict(status='FAIL')
        save_json(root/'execution_progress.json',dict(completed=len(records),total=len(jobs),last_case=name,last_status=result['status']))
        if result['status']=='FAIL' or r['timed_out']:
            save_json(root/'campaign_stop.json',dict(status='FAIL',reason='native/numerical failure; dependent calculations stopped',case=name));break
    print('CAMPAIGN_RECORDED',len(records),'OF',len(jobs),flush=True)


if __name__=='__main__':main()
