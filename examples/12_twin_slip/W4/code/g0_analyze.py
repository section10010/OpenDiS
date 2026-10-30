"""Evidence-based W4-G0 review, including small-contrast error propagation."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from configuration import save_json
from g0_campaign import FAMILIES,MASKS


def read(p):return json.loads(Path(p).read_text())
def csvwrite(p,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with Path(p).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def summarize(folder,audit=False):
    res=read(folder/'result.json');c=read(folder/'config.json')
    row=dict(case_id=folder.name,status=res['status'],time_ns=res.get('time_ns'),steps=res['steps'],target_end_ns=c['end_ns'],
        reached_configured_end=abs(res.get('time_ns',-1)-c['end_ns'])<=1e-9,stop_reason=res.get('stop_reason'))
    if res['status']!='PASS':return row
    d=read(folder/'velocity_diagnostics.json')
    with (folder/'history.csv').open() as stream:h=list(csv.DictReader(stream))
    last=h[-1];row.update(area_nm2=res['swept_area_nm2'],line_nm=res['final']['line_nm'],
        tail_nm=sum(float(last[w+'_tail_nm']) for w in ('lower','upper')),
        cap_line_time_fraction=d['line_time_cap_fraction'],cap_any_record_fraction=d['any_node_capped_record_fraction'],
        max_raw_velocity_m_s=d['maximum_raw_velocity_m_s'],max_constrained_velocity_m_s=d['maximum_constrained_velocity_m_s'],
        max_contact_power_residual=d['max_normalized_contact_power'],max_contact_normal_residual=d['max_normalized_contact_normal_violation'],
        max_burgers_residual=res['max_burgers_residual_b'],plastic_closure=res['max_plastic_area_closure_absolute'],
        external_work_J=d['external_work_from_native_plastic_J'],bulk_dissipation_left_J=d['left_rule_bulk_dissipation_J'],
        reaction_work_J=d['contact_reaction_work_J'])
    with np.load(folder/'initial.npz') as initial:
        b=initial['burgers'][0];n=initial['planes'][0]
        S=.5*(np.outer(b,n)+np.outer(n,b))
        gamma=float(np.einsum('ij,ij->',np.array(res['plastic_strain']),S)/np.sum(S*S))
        row['native_plastic_equivalent_area_nm2']=-gamma*np.prod(c['box_nm'])/(c['b_m']*1e9)
    for wall in ('lower','upper'):
        row[wall+'_first_contact_ns']=next((float(x['time_ns']) for x in h if float(x[wall+'_contact_nodes'])>0),None)
    row['sample_count_matches_steps']=d['accepted_samples']==res['steps']
    row['stage_checks_complete']=all(res['contact_checks'].get(k,0)>=res['steps'] for k in ('integration_substep','collision','topology','remesh'))
    row['contact_admissible']=res['minimum_contact_gap_nm'] is None or res['minimum_contact_gap_nm']>=-c['contact_tol_nm']
    if audit:
        total=0.;hit4000=0.;applied=0.;samples=0;anycap=0;max_cap_error=0.;max_stored_frame_weight_error=0.
        frame_weights={}
        from diagnostics import node_weights
        for fp in folder.glob('frame_*.npz'):
            with np.load(fp) as z:
                frame_weights[int(fp.stem.split('_')[1])+1]=(z['tags'].copy(),node_weights(z))
        for name in d['chunks']:
            with np.load(folder/'accepted_velocity'/name) as archive:
                z={key:archive[key] for key in archive.files}
                for k,dt in enumerate(z['dt_s']):
                    lo,hi=z['offsets'][k:k+2];w=z['node_line_weight_nm'][lo:hi];pre=z['constrained_velocity_m_s'][lo:hi];post=z['capped_velocity_m_s'][lo:hi]
                    norm=np.linalg.norm(pre,axis=1);hit=norm>4000.;total+=float(w.sum()*dt);hit4000+=float(w@hit*dt)
                    limit=c['vmax_m_s'];flag=np.zeros(len(w),bool) if limit is None else norm>limit
                    applied+=float(w@flag*dt);anycap+=bool(flag.any());samples+=1
                    expected=pre if limit is None else pre*np.minimum(1,limit/np.maximum(norm,1e-99))[:,None]
                    max_cap_error=max(max_cap_error,float(np.linalg.norm(expected-post,axis=1).max(initial=0)))
                    witness=frame_weights.get(int(z['step'][k]))
                    if witness is not None:
                        tags,ww=witness;assert np.array_equal(tags,z['tags'][lo:hi])
                        max_stored_frame_weight_error=max(max_stored_frame_weight_error,float(np.max(np.abs(ww-w))))
                assert np.min(z['multipliers_N'],initial=0)>=0
        assert samples==res['steps'] and max_cap_error<1e-7
        fraction=applied/total if total else None
        assert fraction is not None and abs(fraction-d['line_time_cap_fraction'])<1e-12
        row.update(uncapped_reference4000_line_time_fraction=hit4000/total,velocity_archive_audit='PASS',
            max_cap_reconstruction_error_m_s=max_cap_error,max_frame_weight_reconstruction_error_nm=max_stored_frame_weight_error)
    return row


def relative(coarse,fine):
    return abs(coarse-fine)/abs(fine) if fine else (0. if coarse==0 else None)


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--audit-velocities',action='store_true');p.add_argument('--prior-velocity-audit',type=Path);a=p.parse_args();r=a.root;o=a.out;o.mkdir(parents=True,exist_ok=True)
    prior={x['case_id']:x for x in read(a.prior_velocity_audit/'case_metrics.json')} if a.prior_velocity_audit else {}
    manifest=read(r/'campaign_manifest.json');rows=[];jobs=manifest['jobs']
    if (r/'extension_manifest.json').exists():jobs=jobs+read(r/'extension_manifest.json')['jobs']
    for job in jobs:
        folder=r/'cases'/job['case_id']
        if (folder/'result.json').exists():
            row={**job,**summarize(folder,a.audit_velocities)}
            if not a.audit_velocities and row['case_id'] in prior:
                old=prior[row['case_id']]
                assert all(old[k]==row[k] for k in ('status','time_ns','steps','area_nm2'))
                assert old.get('velocity_archive_audit')=='PASS'
                for k in ('uncapped_reference4000_line_time_fraction','velocity_archive_audit','max_cap_reconstruction_error_m_s','max_frame_weight_reconstruction_error_nm'):row[k]=old[k]
                row['velocity_audit_source']=str(a.prior_velocity_audit)
            rows.append(row)
        else:rows.append({**job,'status':'NOT_RUN'})
    csvwrite(o/'case_metrics.csv',rows);save_json(o/'case_metrics.json',rows);by={x['case_id']:x for x in rows}
    def case(kind,mode,mask,label):return by[f'G0_{kind}_{mode}_{mask}_{label}']
    contrasts=[];refinements=[];pairchecks=[];scaling=[];capchecks=[]
    for kind,mode in FAMILIES:
        labels=sorted({x['case_id'].split(x['wall']+'_',1)[1] for x in rows if x['kind']==kind and x['mode']==mode})
        for label in labels:
            if not all(f'G0_{kind}_{mode}_{mask}_{label}' in by for mask in MASKS):continue
            group=[case(kind,mode,mask,label) for mask in MASKS]
            if not all(x['status']=='PASS' and x['reached_configured_end'] for x in group):continue
            initial=[np.load(r/'cases'/x['case_id']/'initial.npz') for x in group]
            fields=('positions_nm','tags','constraints','ids','burgers','planes','cell_nm')
            equal=all(np.array_equal(initial[0][k],z[k]) for z in initial for k in fields)
            for z in initial:z.close()
            configs=[read(r/'cases'/x['case_id']/'config.json') for x in group]
            params=all({k:v for k,v in configs[0].items() if k not in ('name','walls')}=={k:v for k,v in c.items() if k not in ('name','walls')} for c in configs)
            assert equal and params
            pairchecks.append(dict(kind=kind,mode=mode,label=label,identical_initial_arrays=equal,other_parameters_equal=params))
            q=[x['area_nm2'] for x in group]
            plastic_q=[x['native_plastic_equivalent_area_nm2'] for x in group]
            contrasts.append(dict(kind=kind,mode=mode,label=label,none_nm2=q[0],lower_nm2=q[1],upper_nm2=q[2],both_nm2=q[3],
                second_upper_nm2=q[3]-q[1],second_lower_nm2=q[3]-q[2],nonadditivity_nm2=q[3]-q[1]-q[2]+q[0],
                native_plastic_second_upper_nm2=plastic_q[3]-plastic_q[1],native_plastic_nonadditivity_nm2=plastic_q[3]-plastic_q[1]-plastic_q[2]+plastic_q[0],
                reduction_vs_lower=1-q[3]/q[1]))
        for mask in MASKS:
            comparisons=[('space_20_to_10','uncapped_diagnostic','mesh10'),('time_half_at_10','mesh10','mesh10_timehalf')]
            if f'G0_{kind}_{mode}_{mask}_mesh5_timehalf' in by:
                comparisons += [('space_10_to_5','mesh10_timehalf','mesh5_timehalf'),('time_half_at_5','mesh5_timehalf','mesh5_timequarter')]
            for comparison,left,right in comparisons:
                x=case(kind,mode,mask,left);y=case(kind,mode,mask,right)
                if x['status']!='PASS' or y['status']!='PASS':continue
                err={key:relative(x[key],y[key]) for key in ('area_nm2','line_nm','tail_nm','native_plastic_equivalent_area_nm2')}
                ev={wall:relative(x[wall+'_first_contact_ns'],y[wall+'_first_contact_ns']) for wall in ('lower','upper') if x[wall+'_first_contact_ns'] is not None and y[wall+'_first_contact_ns'] is not None}
                missing=[wall for wall in ('lower','upper') if (x[wall+'_first_contact_ns'] is None)!=(y[wall+'_first_contact_ns'] is None) or ((mask=='both' or wall==mask) and wall not in ev)]
                passed=not missing and all(v is not None and v<=.02 for v in err.values()) and all(v is not None and v<=.05 for v in ev.values())
                refinements.append(dict(kind=kind,mode=mode,wall=mask,comparison=comparison,relative_errors=err,event_relative_errors=ev,missing_contact_events=missing,status='PASS' if passed else 'UNRESOLVED'))
            x=case(kind,mode,mask,'uncapped_diagnostic');y=case(kind,mode,mask,'mobility_half_rescaled')
            if x['status']=='PASS' and y['status']=='PASS':
                error=relative(x['area_nm2'],y['area_nm2']);scaling.append(dict(kind=kind,mode=mode,wall=mask,
                    relative_area_error=error,status='PASS' if error is not None and error<=.02 else 'FAIL',
                    convention='M/2, t*2 and maxdt*2; geometric-clock test, not a second material mechanism'))
            x=case(kind,mode,mask,'cap4000');y=case(kind,mode,mask,'uncapped_diagnostic')
            if x['status']=='PASS' and y['status']=='PASS':capchecks.append(dict(kind=kind,mode=mode,wall=mask,
                relative_area_change=relative(x['area_nm2'],y['area_nm2']),
                timing_cap_insensitive='PASS' if relative(x['area_nm2'],y['area_nm2'])<=.02 else 'NOT_QUALIFIED'))
    csvwrite(o/'four_wall_contrasts.csv',contrasts);save_json(o/'four_wall_contrasts.json',contrasts)
    save_json(o/'initial_pair_checks.json',pairchecks);save_json(o/'primary_refinement.json',refinements);save_json(o/'mobility_clock.json',scaling);save_json(o/'cap_sensitivity.json',capchecks)
    resolution=[]
    for kind,mode in FAMILIES:
        extended=f'G0_{kind}_{mode}_none_mesh5_timequarter' in by
        fine_label='mesh5_timequarter' if extended else 'mesh10_timehalf'
        finest=next((x for x in contrasts if x['kind']==kind and x['mode']==mode and x['label']==fine_label),None)
        if finest is None:continue
        errors={}
        for mask in MASKS:
            x=case(kind,mode,mask,'mesh10_timehalf' if extended else 'uncapped_diagnostic');y=case(kind,mode,mask,'mesh5_timehalf' if extended else 'mesh10');z=case(kind,mode,mask,fine_label)
            errors[mask]=abs(x['area_nm2']-y['area_nm2'])+abs(y['area_nm2']-z['area_nm2'])
        bound=sum(errors.values());second_bound=errors['both']+errors['lower']
        resolution.append(dict(kind=kind,mode=mode,finest_label=fine_label,finest_nonadditivity_nm2=finest['nonadditivity_nm2'],observed_refinement_error_sum_nm2=bound,
            nonadditivity_numerical_sign='PASS' if abs(finest['nonadditivity_nm2'])>bound else 'UNRESOLVED',
            second_upper_nm2=finest['second_upper_nm2'],second_refinement_error_sum_nm2=second_bound,
            second_wall_detectable=abs(finest['second_upper_nm2'])>second_bound,
            caveat='Observed refinement differences are an empirical error estimate, not a rigorous theorem or seed confidence interval'))
    save_json(o/'contrast_resolution.json',resolution)
    counts={s:sum(x['status']==s for x in rows) for s in ('PASS','FAIL','CENSORED','NOT_RUN')}
    invariants=all(x.get('sample_count_matches_steps') and x.get('stage_checks_complete') and x.get('contact_admissible') and x['max_burgers_residual']<=1e-8 and x['max_contact_power_residual']<=1e-8 for x in rows if x['status']=='PASS')
    complete=counts['PASS']==len(jobs)
    endpoints=all(x.get('reached_configured_end',False) for x in rows if x['status']=='PASS')
    release_path=r/'dynamic_release_review_r2.json'
    release_status=read(release_path)['status'] if release_path.exists() else 'NOT_RUN'
    latest=[]
    for kind,mode in FAMILIES:
        extended=f'G0_{kind}_{mode}_none_mesh5_timequarter' in by
        comparisons=('space_10_to_5','time_half_at_5') if extended else ('space_20_to_10','time_half_at_10')
        latest += [x for x in refinements if x['kind']==kind and x['mode']==mode and x['comparison'] in comparisons]
    refinement_complete=len(latest)==16 and all(x['status']=='PASS' for x in latest)
    second_robust=bool(len(resolution)==2 and all(x['second_wall_detectable'] for x in resolution) and all(x['second_upper_nm2']<0 for x in contrasts if x['label'] in ('cap4000','uncapped_diagnostic','tau250_uncapped')))
    overall='FAIL' if counts['FAIL'] or not invariants or release_status=='FAIL' else 'CENSORED' if counts['CENSORED'] else 'NOT_RUN' if not complete else 'PASS' if refinement_complete and second_robust and endpoints and release_status=='PASS' else 'UNRESOLVED'
    save_json(o/'gate_review.json',dict(gate='W4-G0',status=overall,case_counts=counts,screen_completed=sum(x['status']=='PASS' and x['group']=='screen' for x in rows),
        native_numerical_invariants='PASS' if invariants else 'FAIL',primary_refinement='PASS' if refinement_complete else 'UNRESOLVED',
        configured_endpoints_matched='PASS' if endpoints else 'UNRESOLVED',
        dynamic_release_regression=release_status,
        second_wall_mechanism='PASS' if second_robust else 'UNRESOLVED',
        actual_time_cap_insensitivity='PASS' if capchecks and all(x['timing_cap_insensitive']=='PASS' for x in capchecks) else 'NOT_QUALIFIED',
        nonlinear_high_speed_material_law='NOT_QUALIFIED',VCoNi='NOT_QUALIFIED',
        absolute_elastic_potential_energy='NOT_RUN; force binding has no energy API; reaction-work/dissipation and native plastic closure audited',
        local_topological_plastic_field='NOT_RUN',G1_admission='PASS' if overall=='PASS' else 'BLOCKED',
        scope='Effective FCC I0 numerical/causal mechanism only; no cap-independent physical timing, thickness trend or material-strength qualification'))
    print(json.dumps(read(o/'gate_review.json'),ensure_ascii=False))


if __name__=='__main__':main()
