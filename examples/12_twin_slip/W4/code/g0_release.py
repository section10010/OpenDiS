"""Forward-final network, opposite stress: real native release regressions."""
import argparse,sys
import csv
from pathlib import Path
import numpy as np
from g0_analyze import read
from configuration import config,save_json
from record import run_logged


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--repo',required=True);a=p.parse_args();r=a.root;code=Path(__file__).resolve().parent;checks=[]
    for kind,mode in [('loop','HI'),('threader','HII')]:
        forward=r/'cases'/f'G0_{kind}_{mode}_both_uncapped_diagnostic';c=read(forward/'config.json')
        name=f'control_release_{kind}_{mode}';c=config(**{**c,'name':name,'stress_pa':(-np.array(c['stress_pa'])).tolist(),'end_ns':.02,'frame_interval_ns':.002})
        cfg=r/(name+'_config.json')
        if cfg.exists():assert read(cfg)==c
        else:save_json(cfg,c)
        cmd=[sys.executable,str(code/'run.py'),'--repo',a.repo,'--a0',str(r/'A0/gate_result.json'),'--config',str(cfg),'--out',str(r/name),'--initial',str(forward/'final.npz')]
        if not (r/name/'result.json').exists():
            rec=run_logged(cmd,code,r/(name+'_receipt.json'),r/(name+'.log'),120);assert rec['exit_code']==0
        result=read(r/name/'result.json')
        with np.load(forward/'final.npz') as before,np.load(r/name/'initial.npz') as initial,np.load(r/name/'final.npz') as after:
            same=all(np.array_equal(before[k],initial[k]) for k in ('positions_nm','tags','constraints','ids','burgers','planes','cell_nm'))
            gaps=[float((w['inward']*(after['positions_nm'][:,2]-w['z_nm'])).min()) for w in c['walls']] if len(after['positions_nm']) else None
            with (r/name/'history.csv').open() as f:history=list(csv.DictReader(f))
            released=[x for x in history if float(x['lower_contact_nodes'])==0 and float(x['upper_contact_nodes'])==0 and min(float(x['lower_gap_nm']),float(x['upper_gap_nm']))>1e-6]
            from run import area_tensor
            geom=-(c['b_m']*1e9)*(area_tensor(after)-area_tensor(initial))/np.prod(c['box_nm']);geom=(geom+geom.T)/2
            closure=float(np.linalg.norm(np.array(result['plastic_strain'])-geom))
            assert same and released and result['status']=='PASS' and closure<=1e-12+1e-8*np.linalg.norm(geom)
            if gaps is not None:assert min(gaps)>1e-6
            else:assert result['stop_reason']=='network_annihilated'
            checks.append(dict(kind=kind,mode=mode,status='PASS',identical_actual_restart_network=same,final_wall_gaps_nm=gaps,
                terminal_state=result['stop_reason'],released_nonempty_states=len(released),first_released_time_ns=float(released[0]['time_ns']),
                final_native_plastic_closure=closure,empty_gap_convention='Undefined/null, never infinity or fabricated zero'))
    save_json(r/'dynamic_release_review_r2.json',dict(status='PASS',checks=checks));print('W4_DYNAMIC_RELEASE_PASS',len(checks))


if __name__=='__main__':main()
