"""Read actual native arrays/accepted-step histories; no prescribed thickness trends."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from configuration import save_json
from geometry import choose
from events import topology


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def arrays(path):return dict(np.load(path))
def history(folder):
    with (folder/'history.csv').open() as f:return [{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
def vector_area(a):
    if not len(a['positions_nm']):return np.zeros(3)
    x=a['positions_nm']-a['positions_nm'].mean(axis=0);ij=a['ids']
    return .5*np.cross(x[ij[:,0]],x[ij[:,1]]).sum(axis=0)
def measures(folder):
    cfg=read(folder/'config.json');r=read(folder/'result.json');a=arrays(folder/'final.npz');b=arrays(folder/'initial.npz')
    h=history(folder);ss=choose(cfg);n=ss['n'];area=float((vector_area(a)-vector_area(b))@n)
    common=np.intersect1d(b['tags'][:,1],a['tags'][:,1]);dx=[]
    for tag in common:
        i=np.flatnonzero(b['tags'][:,1]==tag)[0];j=np.flatnonzero(a['tags'][:,1]==tag)[0]
        if b['constraints'].ravel()[i]==0:dx.append(np.linalg.norm(a['positions_nm'][j]-b['positions_nm'][i]))
    d=dict(case=folder.name,status=r['status'],time_ns=r.get('time_ns'),steps=r['steps'],
        signed_area_gain_nm2=area,swept_area_nm2=r.get('swept_area_nm2'),line_nm=r.get('final',{}).get('line_nm'),
        initial_line_nm=r.get('initial',{}).get('line_nm'),sweep_equivalent_distance_nm=r.get('swept_area_nm2',0)/r.get('initial',{}).get('mobile_line_nm',1),
        retained_original_mobile_nodes=len(dx),retained_node_mean_displacement_nm=float(np.mean(dx)) if dx else None,
        plastic_area_closure=r.get('max_plastic_area_closure_absolute'),minimum_gap_nm=r.get('minimum_contact_gap_nm'))
    d['max_slip_plane_deviation_nm']=float(np.abs((a['positions_nm']-np.array(cfg['box_nm'])/2)@n).max(initial=0))
    for k,v in h[-1].items():
        if '_tail_nm' in k or '_band10_line_nm' in k or '_contact_nodes' in k or '_gap_nm' in k:d[k]=v
        if '_band10_line_nm' in k:d[k.replace('_line_nm','_rho_m2')]=v/(cfg['box_nm'][0]*cfg['box_nm'][1]*20)*1e18
    for wall in cfg['walls']:
        key=wall['interface_id']+'_contact_nodes';d[wall['interface_id']+'_first_contact_ns']=next((x['time_ns'] for x in h if x.get(key,0)>0),None)
    return d


def export_csv(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    r=a.root;o=a.out;o.mkdir(exist_ok=True,parents=True);rows=[];families=[];all_pair=True;all_contact=True
    fields=('positions_nm','tags','constraints','ids','burgers','planes','cell_nm')
    for kind in ('loop','threader'):
        for mode in ('HI','HII'):
            folders={m:r/f'B_{kind}_{mode}_{m}_tau300_t100' for m in ('none','lower','upper','both')}
            if not all((x/'result.json').exists() for x in folders.values()):continue
            initial={m:arrays(x/'initial.npz') for m,x in folders.items()};configs={m:read(x/'config.json') for m,x in folders.items()}
            equal=all(np.array_equal(initial['none'][k],b[k]) for b in initial.values() for k in fields)
            equal_cfg=all({k:v for k,v in configs['none'].items() if k!='walls'}=={k:v for k,v in c.items() if k!='walls'} for c in configs.values())
            numerical=all(read(x/'result.json')['status']=='PASS' and abs(read(x/'result.json')['time_ns']-.1)<1e-9 for x in folders.values())
            all_pair &= equal and equal_cfg and numerical
            d={m:measures(x) for m,x in folders.items()};rows.extend(d.values())
            family=dict(kind=kind,mode=mode,identical_initial_arrays=equal,identical_other_configuration=equal_cfg,numerical_all_PASS=numerical)
            for metric in ('swept_area_nm2','signed_area_gain_nm2','line_nm'):
                q={m:x[metric] for m,x in d.items()}
                family[metric+'_second_upper_increment']=q['both']-q['lower']
                family[metric+'_second_lower_increment']=q['both']-q['upper']
                family[metric+'_nonadditivity']=q['both']-q['lower']-q['upper']+q['none']
                if metric=='swept_area_nm2':family['second_upper_reduction_fraction']=1-q['both']/q['lower']
            for m in ('lower','upper','both'):
                result=read(folders[m]/'result.json');checks=result.get('contact_checks',{})
                valid=result['minimum_contact_gap_nm']>=-configs[m]['contact_tol_nm'] and all(checks.get(s,0)>=result['steps'] for s in ('integration_substep','collision','topology','remesh'))
                all_contact &= valid
            families.append(family)
    export_csv(o/'paired_metrics.csv',rows);save_json(o/'causal_effects.json',families)
    refinements=[]
    for kind in ('loop','threader'):
        for mode in ('HI','HII'):
            folder=r/f'B_{kind}_{mode}_both_refine_t100'
            if not (folder/'result.json').exists():continue
            d=measures(folder);base=next(x for x in rows if x['case']==f'B_{kind}_{mode}_both_tau300_t100')
            err={k:abs(d[k]/base[k]-1) for k in ('swept_area_nm2','signed_area_gain_nm2','line_nm')}
            refinements.append(dict(kind=kind,mode=mode,relative_errors=err,status='PASS' if max(err.values())<.05 and d['status']=='PASS' else 'FAIL'))
    save_json(o/'time_refinement.json',refinements)
    release=[]
    for mode in ('HI','HII'):
        folder=r/f'B_loop_{mode}_release'
        if not (folder/'result.json').exists():continue
        d=measures(folder);baseline=arrays(r/f'B_loop_{mode}_both_tau300_t100/final.npz');init=arrays(folder/'initial.npz')
        same=all(np.array_equal(baseline[k],init[k]) for k in fields)
        passed=d['status']=='PASS' and same and d['lower_contact_nodes']==0 and d['upper_contact_nodes']==0 and min(d['lower_gap_nm'],d['upper_gap_nm'])>1e-5
        release.append(dict(mode=mode,same_actual_restart_network=same,lower_gap_nm=d['lower_gap_nm'],upper_gap_nm=d['upper_gap_nm'],status='PASS' if passed else 'FAIL'))
    save_json(o/'dynamic_release.json',release)
    source=[]
    patterns=('A1_bracket_tau*_unique','A1_ref240_80_tau200_wide','A1_mesh20_first_junctionaware','A1_sensitivity*')
    for pattern in patterns:
        for folder in sorted(r.glob(pattern)):
            if not (folder/'result.json').exists():continue
            cfg=read(folder/'config.json');result=read(folder/'result.json');lineage=read(folder/'source_lineage.json')
            children=[x for x in lineage['children'] if x['confirmed']];ev=lineage['events'];h=history(folder)
            source.append(dict(case=folder.name,numerical_status=result['status'],tau_mpa=abs(choose(cfg)['rss_pa'])/1e6,
                gauge_nm=cfg['source_gauge_nm'],return_nm=cfg['backbone_offset_nm'],maxseg_nm=cfg['maxseg_nm'],
                time_ns=result.get('time_ns'),distinct_qualified_loops=result.get('confirmed_emissions',0),
                first_detachment_ns=min((x['birth_time_ns'] for x in children),default=None),
                confirmed_lineage_ids=[x['lineage_id'] for x in children],final_max_v_m_s=h[-1]['velocity_max_m_s']))
    export_csv(o/'source_summary.csv',source);save_json(o/'source_summary.json',source)
    source_case=r/'A1_ref240_80_tau200_wide';a1='NOT_RUN';source_validation={}
    if (source_case/'result.json').exists():
        cfg=read(source_case/'config.json');res=read(source_case/'result.json');lin=read(source_case/'source_lineage.json');aa=arrays(source_case/'final.npz')
        groups=topology(aa);radii=[float(np.sqrt(abs(vector_area(dict(aa,ids=aa['ids'][np.isin(aa['ids'][:,0],g['ix'])]))@choose(cfg)['n'])/np.pi)) for g in groups if g['free_cycle']]
        # Independent final geometry witness, separate from event counters.
        inside=bool(np.all(aa['positions_nm']>=0) and np.all(aa['positions_nm']<=np.array(cfg['box_nm'])))
        source_validation=dict(large_free_loop_count=sum(x>=cfg['source_gauge_nm']/2 for x in radii),free_loop_area_radii_nm=radii,
            mobile_anchored_parent=any(g['anchored'] and np.any(aa['constraints'].ravel()[g['ix']]!=7) for g in groups),
            inside_computational_box=inside,confirmed_event_count=lin['confirmed_detachments'])
        fine=next((x for x in source if x['case']=='A1_mesh20_first_junctionaware'),None)
        reference=next(x for x in source if x['case']=='A1_ref240_80_tau200_wide')
        source_validation['fine20_first_emission_reproduced']=bool(fine and fine['distinct_qualified_loops']>=1)
        source_validation['first_detachment_time_mesh_relative_difference']=abs(fine['first_detachment_ns']/reference['first_detachment_ns']-1) if fine and fine['first_detachment_ns'] else None
        sensitivity=[x for x in source if 'sensitivity' in x['case']]
        source_validation['four_geometric_sensitivities_completed']=len(sensitivity)==4 and all(x['numerical_status']=='PASS' and x['distinct_qualified_loops']>=1 for x in sensitivity)
        a1='PASS' if res['status']=='PASS' and res['confirmed_emissions']>=3 and source_validation['large_free_loop_count']>=3 and source_validation['mobile_anchored_parent'] and inside else 'FAIL'
        if not source_validation['fine20_first_emission_reproduced'] or source_validation['first_detachment_time_mesh_relative_difference']>.05:a1='FAIL'
        if not source_validation['four_geometric_sensitivities_completed']:a1='FAIL'
    save_json(o/'source_independent_geometry.json',source_validation)
    b0=('PASS' if all(x['status']=='PASS' and x['signed_area_gain_nm2']>0 for x in rows if '_loop_' in x['case'] and '_none_' in x['case']) else 'FAIL') if len(rows)==16 else 'NOT_RUN'
    b1=('PASS' if all_pair and all_contact and all(x['status']=='PASS' for x in refinements+release) and read(r/'B_contact_unit_final/result.json')['status']=='PASS' else 'FAIL') if len(families)==4 and len(refinements)==4 and len(release)==2 else 'NOT_RUN'
    save_json(o/'gate_report.json',dict(A0=read(r/'A0/gate_result.json')['status'],A1=a1,B0=b0,B1=b1,
        B1_scope='I0, nonperiodic single channel, existing-loop / prepared-threader propagation; not material strengthening',
        VCoNi_material_qualification='NOT_QUALIFIED',W4_four_thickness_campaign='NOT_RUN',same_TB_two_crystal_sides='NOT_RUN',
        source_critical_stress_scope='Independent constant-stress, two-ns activation window; no asymptotic critical-stress claim',
        second_wall_changes_existing_propagation=bool(len(families)==4 and all_pair and all(abs(x['swept_area_nm2_second_upper_increment'])>1 for x in families))))
    print(json.dumps(read(o/'gate_report.json'),ensure_ascii=False))


if __name__=='__main__':main()
