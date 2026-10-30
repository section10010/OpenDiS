"""Bounded numerical refinement, independent radius checks and dynamic release."""
import sys,json
from pathlib import Path
import numpy as np
from campaign import execute,FAILURES
from configuration import config
from geometry import pure_shear


root=Path(sys.argv[1]);group=sys.argv[2]
if group=='fine_source':
    c=config(name='A1_mesh20_first',gate='A1',stress_pa=pure_shear('matrix','HI',200e6).tolist(),
             end_ns=.5,stop_emissions=1,frame_interval_ns=.025,max_wall_s=300.)
    execute(root,'A1_mesh20_first_junctionaware',c)
elif group=='wide_source':
    c=config(name='A1_wide_domain',gate='A1',box_nm=[16000.]*3,stress_pa=pure_shear('matrix','HI',200e6).tolist(),
        maxseg_nm=40.,minseg_nm=4.,end_ns=1.5,stop_emissions=3,frame_interval_ns=.025,max_wall_s=900.)
    execute(root,'A1_ref240_80_tau200_wide',c)
elif group=='B_refine':
    for kind in ('loop','threader'):
        for mode in ('HI','HII'):
            base=root/f'B_{kind}_{mode}_both_tau300_t100'
            c=json.loads((base/'config.json').read_text());c.update(max_node_step_nm=.125,maxdt_ps=.25)
            execute(root,f'B_{kind}_{mode}_both_refine_t100',config(**c))
elif group=='release':
    for mode in ('HI','HII'):
        base=root/f'B_loop_{mode}_both_tau300_t100';c=json.loads((base/'config.json').read_text())
        c.update(stress_pa=(-np.array(c['stress_pa'])).tolist(),end_ns=.02,frame_interval_ns=.002)
        execute(root,f'B_loop_{mode}_release',config(**c),initial=base/'final.npz')
elif group=='radius':
    base=root/'B_loop_HI_both_tau300_t100';c=json.loads((base/'config.json').read_text())
    for r in (40.,80.):
        for mask in ('none','both'):
            cc={**c,'loop_radius_nm':r,'walls':[{**w,'enabled':mask=='both'} for w in c['walls']]}
            execute(root,f'B_loop_HI_radius{int(r)}_{mask}',config(**cc))
else:raise ValueError(group)
sys.exit(1 if FAILURES else 0)
