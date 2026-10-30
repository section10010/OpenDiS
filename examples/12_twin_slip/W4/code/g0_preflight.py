"""Real native positive, illegal-contact negative, and resource-censor controls."""
import argparse,json,sys
from pathlib import Path
from configuration import config,save_json
from record import run_logged


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--repo',required=True);a=p.parse_args()
    r=a.root;code=Path(__file__).resolve().parent;c=json.loads((r/'configs/G0_loop_HI_both_cap4000.json').read_text());checks=[]
    for label,changes,expected in [('positive',dict(end_ns=.0002),'PASS'),('resource',dict(max_steps=1),'CENSORED'),
        ('illegal_contact',dict(walls=[{**w,'z_nm':1980. if w['interface_id']=='lower' else w['z_nm']} for w in c['walls']]),'FAIL')]:
        name='preflight_'+label;cfg=config(**{**c,**changes,'name':name});path=r/(name+'_config.json');save_json(path,cfg)
        cmd=[sys.executable,str(code/'run.py'),'--repo',a.repo,'--config',str(path),'--out',str(r/name),'--a0',str(r/'A0/gate_result.json')]
        receipt=run_logged(cmd,code,r/(name+'_receipt.json'),r/(name+'.log'),60)
        result=json.loads((r/name/'result.json').read_text());assert result['status']==expected,(name,result)
        if expected=='CENSORED':assert (r/name/'checkpoint.npz').exists()
        if expected=='FAIL':assert receipt['exit_code']!=0 and (r/name/'failure.npz').exists()
        checks.append(dict(case=name,status='PASS',expected_native_status=expected,actual_native_status=result['status']))
    save_json(r/'preflight_controls.json',dict(status='PASS',checks=checks));print('REAL_PREFLIGHT_CONTROLS_PASS',len(checks))


if __name__=='__main__':main()
