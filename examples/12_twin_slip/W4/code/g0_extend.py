"""Bounded resolution repair after an actual UNRESOLVED coarse-grid review.

Does not relax criteria, alter original cases, or change material/source physics.
Selects only families with failed primary refinement and preregisters finer runs.
"""
import argparse,time,sys
import numpy as np
from pathlib import Path
from configuration import config,save_json
from g0_analyze import read
from g0_campaign import MASKS
from record import run_logged
from geometry import build


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--review',type=Path,required=True);p.add_argument('--repo',required=True);a=p.parse_args();r=a.root;code=Path(__file__).resolve().parent
    families=sorted({(x['kind'],x['mode']) for x in read(a.review/'primary_refinement.json') if x['status']=='UNRESOLVED'})
    assert read(a.review/'gate_review.json')['screen_completed']==24
    manifest=r/'extension_manifest.json'
    if manifest.exists():raise RuntimeError('Refuse overwriting extension; create a separately reviewed amendment')
    jobs=[]
    for kind,mode in families:
        for wall in MASKS:
            base=read(r/'configs'/f'G0_{kind}_{mode}_{wall}_mesh10_timehalf.json')
            for label,time_changes in [('mesh5_timehalf',{}),('mesh5_timequarter',dict(maxdt_ps=.125,max_node_step_nm=.0625))]:
                name=f'G0_{kind}_{mode}_{wall}_{label}';c=config(**{**base,'name':name,'maxseg_nm':5.,'minseg_nm':2.,**time_changes})
                save_json(r/'configs'/(name+'.json'),c);jobs.append(dict(case_id=name,kind=kind,mode=mode,wall=wall,group='resolution_repair'))
    seen={}
    for job in jobs:
        label=job['case_id'].split(job['wall']+'_',1)[1];key=(job['kind'],job['mode'],label)
        arrays=build(read(r/'configs'/(job['case_id']+'.json')))[:2]
        if key in seen:assert all(np.array_equal(x,y) for x,y in zip(seen[key],arrays))
        else:seen[key]=arrays
    save_json(manifest,dict(status='FROZEN_BEFORE_EXTENSION',reason='Actual tail observable exceeds unchanged 2 percent criterion',
        prior_review=str(a.review),families=families,jobs=jobs,max_extension_wall_s=1800,
        changes='Coupled spatial resolution 10/4 to 5/2 nm and independent time halving; no material/interface/load/source changes',
        stopping='No further automatic resolution expansion in this amendment; report UNRESOLVED/CENSORED if still insufficient'))
    started=time.monotonic()
    for job in jobs:
        name=job['case_id']
        if time.monotonic()-started>1800:
            save_json(r/'extension_stop.json',dict(status='CENSORED',next_case=name));break
        c=read(r/'configs'/(name+'.json'));cmd=[sys.executable,'-u',str(code/'run.py'),'--repo',a.repo,'--a0',str(r/'A0/gate_result.json'),
             '--config',str(r/'configs'/(name+'.json')),'--out',str(r/'cases'/name)]
        print('EXTENSION',name,flush=True)
        rec=run_logged(cmd,code,r/'logs'/(name+'_receipt.json'),r/'logs'/(name+'.log'),c['max_wall_s']+60)
        result=read(r/'cases'/name/'result.json') if (r/'cases'/name/'result.json').exists() else dict(status='FAIL')
        if result['status']=='FAIL' or rec['timed_out']:break


if __name__=='__main__':main()
