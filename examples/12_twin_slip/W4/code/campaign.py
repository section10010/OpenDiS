"""Bounded native W3 runs. No four-thickness population campaign."""
import argparse,sys
from pathlib import Path
from configuration import config,save_json
from geometry import pure_shear
from record import run_logged
FAILURES=[]


def execute(root,name,c,timeout=None,initial=None):
    code=Path(__file__).resolve().parent;root=Path(root);root.mkdir(exist_ok=True,parents=True)
    (root/'configs').mkdir(exist_ok=True);cfg=root/'configs'/(name+'.json');save_json(cfg,c)
    cmd=[sys.executable,'-u',str(code/'run.py'),'--repo','/home/sakura/opendis/src','--a0',str(root/'A0/gate_result.json'),
         '--config',str(cfg),'--out',str(root/name)]
    if initial:cmd+=['--initial',str(initial)]
    receipt=run_logged(cmd,str(code),root/(name+'_receipt.json'),root/(name+'.log'),timeout or c['max_wall_s']+60)
    if receipt['exit_code']!=0 or receipt['timed_out']:FAILURES.append(name)
    return receipt


def a1(root,group):
    common=dict(gate='A1',stress_pa=pure_shear('matrix','HI',200e6).tolist(),source_gauge_nm=240.,backbone_offset_nm=80.,
                maxseg_nm=40.,minseg_nm=4.,frame_interval_ns=.025,max_wall_s=900.,end_ns=1.5,stop_emissions=3,threads=2)
    if group=='reference':execute(root,'A1_ref240_80_tau200_wide',config(**{**common,'box_nm':[16000.]*3}))
    elif group=='scan':
        for tau in (50,100,150,250,300):
            c=config(**{**common,'stress_pa':pure_shear('matrix','HI',tau*1e6).tolist(),'end_ns':.7,'stop_emissions':1,'max_wall_s':240.})
            execute(root,f'A1_scan_tau{tau}',c)
    elif group=='bracket':
        for tau in (50,75,100,125,150,250,300):
            c=config(**{**common,'stress_pa':pure_shear('matrix','HI',tau*1e6).tolist(),'end_ns':2.,'stop_emissions':1,'max_wall_s':240.})
            execute(root,f'A1_bracket_tau{tau}_unique',c)
    elif group=='geometry':
        for key,values in [('source_gauge_nm',(120.,480.)),('backbone_offset_nm',(40.,160.))]:
            for value in values:
                c=config(**{**common,key:value,'end_ns':.8,'stop_emissions':1,'max_wall_s':300.})
                execute(root,f'A1_sensitivity_{key}_{int(value)}',c)
    else:raise ValueError(group)


def b(root,group):
    minimal=group.endswith('_minimal');group=group.replace('_minimal','')
    for mode in ('HI','HII'):
        walls=[dict(interface_id='lower',z_nm=1900.,inward=1,enabled=False),dict(interface_id='upper',z_nm=2100.,inward=-1,enabled=False)]
        common=dict(gate='B1',kind='propagation_loop' if group=='loop' else 'prepared_threader',mode=mode,
          stress_pa=pure_shear('matrix',mode,300e6).tolist(),walls=walls,integrator='EulerForward',
          maxseg_nm=20.,minseg_nm=4.,maxdt_ps=.5,max_node_step_nm=.25,frame_interval_ns=.025,
          end_ns=.1 if minimal else .25,max_wall_s=300.,threads=2)
        for mask in ('none','lower','upper','both'):
            ww=[{**w,'enabled':mask=='both' or mask==w['interface_id']} for w in walls]
            c=config(**{**common,'walls':ww})
            execute(root,f'B_{group}_{mode}_{mask}_tau300'+('_t100' if minimal else ''),c)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root');p.add_argument('stage',choices=('A1','B'));p.add_argument('group');a=p.parse_args()
    (a1 if a.stage=='A1' else b)(a.root,a.group)
    sys.exit(1 if FAILURES else 0)
