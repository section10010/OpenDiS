"""Independent G1 trajectory/paired-causality review, without seed/frame CIs."""
import argparse,json,csv,sys
from pathlib import Path
import numpy as np


def load(path):return json.loads(Path(path).read_text())


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--spatial-extension',action='store_true');a=p.parse_args()
    root=a.root;out=a.out
    if out.exists():raise RuntimeError('Refuse overwrite')
    out.mkdir(parents=True);results={};metrics={};audits=[];paired=[];refinements=[]
    def save(name,value):(out/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    manifest=load(root/'campaign_manifest.json');jobs=manifest['jobs']
    if a.spatial_extension:
        assert load(root/'spatial_extension_result.json')['status']=='PASS'
        jobs=jobs+load(root/'spatial_extension_manifest.json')['jobs']
    criteria=load(root/'case_specs/base_MT.json')['criteria']
    aliases=load(root/'case_aliases.json').get('aliases',{}) if (root/'case_aliases.json').exists() else {}
    def case_dir(name):return root/'cases'/aliases.get(name,name)
    regional=[];reaction_ledger=[]
    Q=np.array([[1/np.sqrt(2),-1/np.sqrt(2),0],[1/np.sqrt(6),1/np.sqrt(6),-2/np.sqrt(6)],[1/np.sqrt(3)]*3])
    rotation=np.array([[.5,-np.sqrt(3)/2,0],[np.sqrt(3)/2,.5,0],[0,0,1]])
    for job in jobs:
        name=job['case_id'];folder=case_dir(name);d=load(folder/'result.json');results[name]=d
        rows=list(csv.DictReader((folder/'history.csv').open()));final=rows[-1]
        metrics[name]={**{k:float(v) if v else None for k,v in final.items()},'status':d['status'],'steps':d['steps']}
        # Independent crystal/graph audit does not call the producer's identity
        # classifier. Every saved native state is checked, not a selected frame.
        for frame in [folder/'initial.npz',*sorted(folder.glob('frame_*.npz')),folder/'final.npz']:
            with np.load(frame) as z:
                assert np.all(z['node_host'][z['ids'][:,0]]==z['segment_host'])
                assert np.all(z['node_host'][z['ids'][:,1]]==z['segment_host'])
                np.testing.assert_allclose(np.einsum('ij,ij->i',z['burgers'],z['planes']),0,atol=1e-8)
                for h in ('M','T'):
                    edge=z['segment_host']==h
                    if not edge.any():continue
                    R=Q if h=='M' else rotation@Q
                    np.testing.assert_allclose(np.abs(z['planes'][edge]@R),1/np.sqrt(3),atol=1e-8)
                    np.testing.assert_allclose(np.sort(np.abs(z['burgers'][edge]@R),axis=1),np.tile([0,1/np.sqrt(2),1/np.sqrt(2)],(edge.sum(),1)),atol=1e-8)
        source={'M':0,'T':1};dt_sum=0.;samples=0;max_power=0.;max_normal=0.;min_gap=float('inf');frame_evaluation_error=0.;frame_checks=0
        for chunk in sorted((folder/'accepted_velocity').glob('chunk_*.npz')):
            with np.load(chunk) as z:
                assert np.all(z['dt_s']>0);dt_sum+=float(z['dt_s'].sum());samples+=len(z['dt_s'])
                sign=np.where(z['node_host']=='M',-1.,1.);gap=sign*(z['positions_nm'][:,2]-2000.)
                assert np.array_equal(z['node_source_id'],np.array([source[h] for h in z['node_host']]))
                np.testing.assert_array_equal(z['capped_velocity_m_s'],z['constrained_velocity_m_s'])
                reaction=z['reaction_N'];vc=z['constrained_velocity_m_s'];vr=z['raw_velocity_m_s']
                np.testing.assert_allclose(reaction[:,:2],0,atol=1e-30)
                np.testing.assert_allclose(reaction[:,2],sign*z['multipliers_N'][:,0],rtol=1e-12,atol=1e-25)
                assert np.min(z['multipliers_N'],initial=0)>=0
                scale=np.maximum(np.linalg.norm(reaction,axis=1),np.linalg.norm(z['force_N'],axis=1))*np.maximum(np.linalg.norm(vc,axis=1),np.linalg.norm(vr,axis=1))
                power=np.abs(np.einsum('ij,ij->i',reaction,vc))/np.maximum(scale,1e-30)
                max_power=max(max_power,float(power.max(initial=0)))
                if d['contact_enabled']:
                    min_gap=min(min_gap,float(gap.min(initial=float('inf'))));assert gap.min(initial=0)>=-1e-6
                    active=gap<=1e-6
                    normal=np.maximum(0,-sign*vc[:,2])/np.maximum(np.linalg.norm(vr,axis=1),1e-30)
                    max_normal=max(max_normal,float(normal[active].max(initial=0)))
                else:assert np.all(reaction==0)
                for s,dt in enumerate(z['dt_s']):
                    lo,hi=z['offsets'][s:s+2];r=reaction[lo:hi];v=vc[lo:hi];host=z['node_host'][lo:hi]
                    frame=folder/f'frame_{int(z["step"][s])-1:07d}.npz'
                    if frame.exists():
                        with np.load(frame) as actual:
                            np.testing.assert_array_equal(actual['tags'],z['tags'][lo:hi])
                            np.testing.assert_array_equal(actual['positions_nm'],z['positions_nm'][lo:hi])
                            for field in ('raw_velocity_m_s','constrained_velocity_m_s','capped_velocity_m_s'):
                                err=float(np.linalg.norm(actual[field]-z[field][lo:hi])/max(np.linalg.norm(actual[field]),1e-30))
                                frame_evaluation_error=max(frame_evaluation_error,err)
                        frame_checks+=1
                    rm=r[host=='M'].sum(axis=0);rt=r[host=='T'].sum(axis=0);carrier=-(rm+rt)
                    reaction_ledger.append(dict(case_id=name,actual_case=folder.name,step=int(z['step'][s]),time_ns=float(z['time_s'][s])*1e9,dt_s=float(dt),
                        M_reaction_z_N=float(rm[2]),T_reaction_z_N=float(rt[2]),CTB_carrier_z_N=float(carrier[2]),
                        virtual_work_J=float(np.einsum('ij,ij->',r,v)*dt)))
        assert samples==d['steps'];assert abs(dt_sum-d['time_ns']*1e-9)<1e-22
        assert max_power<=1e-8 and max_normal<=1e-8
        assert frame_checks>0 and frame_evaluation_error<=1e-8
        assert d['max_burgers_residual_b']<=1e-8 and d['max_plastic_closure']<1e-12
        assert d['status']=='PASS'
        with np.load(folder/'initial.npz') as ini,np.load(folder/'final.npz') as fin:
            def tensor(z,h):
                edge=z['segment_host']==h;ij=z['ids'][edge];x=z['positions_nm']-np.array([2000.,2000.,2000.])
                return np.einsum('si,sj->ij',z['burgers'][edge],.5*np.cross(x[ij[:,0]],x[ij[:,1]]))
            c=load(folder/'config.json');parts={}
            for h in ('M','T'):
                beta=-c['b_m']*1e9*(tensor(fin,h)-tensor(ini,h))/np.prod(c['box_nm']);parts[h]=.5*(beta+beta.T)
            native=np.array(load(folder/'plastic.json')['native_cumulative']);regional_error=float(np.linalg.norm(parts['M']+parts['T']-native))
            assert regional_error<1e-12
            regional.append(dict(case_id=name,normalization='Each origin contribution uses total simulation volume, so M+T must recover native global strain',
                M=parts['M'].tolist(),T=parts['T'].tolist(),native=native.tolist(),closure=regional_error))
        audits.append(dict(case_id=name,actual_case=folder.name,status='PASS',accepted_samples=samples,dt_sum_ns=dt_sum*1e9,
            max_normalized_power=max_power,max_normalized_normal_violation=max_normal,minimum_gap_nm=None if not np.isfinite(min_gap) else min_gap,
            frame_left_evaluation_relative=frame_evaluation_error,frame_evaluation_checks=frame_checks,
            tag_material_source_audit='PASS',plastic_closure=d['max_plastic_closure']))
    levels=['base','timehalf','timequarter','mesh5','mesh2p5','mesh2p5_timehalf']
    if a.spatial_extension:levels.append('mesh1p25')
    # Every forward refinement has exactly the same initial graph, not merely
    # the same JSON radius. Single populations equal the MT subnetworks.
    reference=None
    for level in levels:
        with np.load(case_dir(level+'_MT')/'initial.npz') as mt:
            if reference is None:reference={k:mt[k].copy() for k in ('positions_nm','tags','ids','burgers','planes')}
            else:
                for k,v in reference.items():np.testing.assert_array_equal(mt[k],v)
            for host in ('M','T'):
                with np.load(case_dir(level+'_'+host)/'initial.npz') as one:
                    mask=mt['node_host']==host;indices=np.flatnonzero(mask);remap=np.full(len(mask),-1);remap[indices]=np.arange(len(indices))
                    for k in ('positions_nm','tags','constraints','node_source_id','origin_region'):np.testing.assert_array_equal(mt[k][mask],one[k])
                    edge=mt['segment_host']==host
                    np.testing.assert_array_equal(remap[mt['ids'][edge]],one['ids'])
                    for k in ('burgers','planes'):np.testing.assert_array_equal(mt[k][edge],one[k])
                paired.append(dict(level=level,population=host,initial_subnetwork='EXACT'))
    pairs=[('timehalf','timequarter','temporal_at_10'),('mesh5','mesh2p5','spatial_at_fine_time'),('mesh2p5','mesh2p5_timehalf','temporal_at_2p5')]
    if a.spatial_extension:pairs.append(('mesh2p5_timehalf','mesh1p25','final_spatial_at_fine_time'))
    for coarse,fine,label in pairs:
        for population in ('M','T','MT'):
            c=metrics[coarse+'_'+population];f=metrics[fine+'_'+population];err={}
            for h in population:
                for suffix in ('line_nm','tail_nm','swept_area_nm2','area_tensor_norm_nm2'):
                    key=h+'_'+suffix;x,y=c[key],f[key];err[key]=abs(x-y)/max(abs(y),1e-20)
            events={h:abs(results[coarse+'_'+population]['first_contact_ns'][h]-results[fine+'_'+population]['first_contact_ns'][h])/results[fine+'_'+population]['first_contact_ns'][h] for h in population}
            refinements.append(dict(label=label,population=population,primary_relative=err,event_relative=events,
                qualification_role='historical_coarse_level' if a.spatial_extension and label=='spatial_at_fine_time' else 'required',
                status='PASS' if max(err.values())<=criteria['primary_refinement_relative'] and max(events.values())<=criteria['event_relative'] else 'UNRESOLVED'))
    contrasts={}
    for level in levels:
        contrasts[level]={}
        for h in ('M','T'):
            single=metrics[level+'_'+h];both=metrics[level+'_MT'];delta=both[h+'_swept_area_nm2']-single[h+'_swept_area_nm2']
            contrasts[level][h]=dict(area_difference_nm2=delta,relative_to_single=delta/single[h+'_swept_area_nm2'],
                tail_difference_nm=both[h+'_tail_nm']-single[h+'_tail_nm'],
                first_contact_difference_ns=results[level+'_MT']['first_contact_ns'][h]-results[level+'_'+h]['first_contact_ns'][h])
    causal=[]
    origin=load(root/'origin_sensitivity_result.json') if (root/'origin_sensitivity_result.json').exists() else dict(status='NOT_RUN',checks=[])
    finest='mesh1p25' if a.spatial_extension else 'mesh2p5_timehalf'
    for h in ('M','T'):
        final=contrasts[finest][h]['area_difference_nm2']
        spatial=abs(final-contrasts['mesh2p5_timehalf'][h]['area_difference_nm2']) if a.spatial_extension else abs(contrasts['mesh2p5'][h]['area_difference_nm2']-contrasts['mesh5'][h]['area_difference_nm2'])
        temporal=abs(contrasts['mesh2p5_timehalf'][h]['area_difference_nm2']-contrasts['mesh2p5'][h]['area_difference_nm2'])
        origin_item=next((v for v in origin['checks'] if v['host']==h),None)
        origin_error=abs(origin_item['base_difference_nm2']-origin_item['shifted_origin_difference_nm2']) if origin_item else None
        error=spatial+temporal+(origin_error if origin_error is not None else 0.)
        signs=[np.sign(contrasts[level][h]['area_difference_nm2']) for level in (('mesh2p5','mesh2p5_timehalf','mesh1p25') if a.spatial_extension else ('mesh5','mesh2p5','mesh2p5_timehalf'))]
        causal.append(dict(host=h,delta_area_nm2=final,observed_spatial_difference_nm2=spatial,observed_time_difference_nm2=temporal,
            empirical_error_nm2=error,status='PASS' if abs(final)>2*error and len(set(signs))==1 else 'UNRESOLVED',
            observed_origin_difference_nm2=origin_error,origin_sensitivity=origin['status'],
            note='Frozen deterministic resolution criterion; empirical differences are neither a rigorous error bound nor a seed CI'))
    releases=[];ctb=[]
    for population in ('M','T','MT'):
        reverse=results['reverse_'+population];forward=results['mesh2p5_timehalf_'+population]
        with np.load(case_dir('reverse_'+population)/'initial.npz') as initial,np.load(case_dir('mesh2p5_timehalf_'+population)/'final.npz') as final:
            for k in ('positions_nm','tags','ids','burgers','planes'):np.testing.assert_array_equal(initial[k],final[k])
        releases.append(dict(population=population,first_full_release_ns=reverse['first_full_release_ns'],
            status='PASS' if all(reverse['first_full_release_ns'][h] is not None for h in population) else 'UNRESOLVED',
            restart_actual_forward_final='EXACT'))
        with np.load(case_dir('no_ctb_'+population)/'initial.npz') as no,np.load(case_dir('mesh2p5_timehalf_'+population)/'initial.npz') as yes:
            for k in ('positions_nm','tags','ids','burgers','planes'):np.testing.assert_array_equal(no[k],yes[k])
        for h in population:
            wall=metrics['mesh2p5_timehalf_'+population][h+'_swept_area_nm2'];free=metrics['no_ctb_'+population][h+'_swept_area_nm2']
            ctb.append(dict(population=population,host=h,wall_area_nm2=wall,no_ctb_area_nm2=free,area_reduction=(free-wall)/free,
                contact_tail_nm=metrics['mesh2p5_timehalf_'+population][h+'_tail_nm'],status='PASS' if wall<free else 'FAIL'))
    force_path=root/('force_audit_mesh1p25' if a.spatial_extension else 'force_audit')/'result.json'
    force=load(force_path) if force_path.exists() else dict(status='NOT_RUN')
    g0=load(root/'g0_regression/result.json') if (root/'g0_regression/result.json').exists() else dict(status='NOT_RUN')
    centered=(root/'centered_force_admission/result.json').exists()
    required=[v for v in refinements if v['qualification_role']=='required']+causal+releases+ctb
    status='FAIL' if force['status']=='FAIL' or g0['status']=='FAIL' or origin['status']=='FAIL' or any(v['status']=='FAIL' for v in ctb) else 'PASS' if (len(audits)==(29 if a.spatial_extension else 26) and all(v['status']=='PASS' for v in required) and force['status']=='PASS' and g0['status']=='PASS' and (not centered or origin['status']=='PASS')) else 'UNRESOLVED'
    save('gate_review.json',dict(gate='W4-G1',status=status,numerical_cases=len(audits),initial_subnetwork_pairing='PASS',
        source_inheritance_contact_release_topology='PASS' if load(root/('admission_r2' if (root/'admission_r2/result.json').exists() else 'admission_r1')/'result.json')['status']=='PASS' else 'FAIL',
        native_force_causality=force['status'],G0_regression=g0['status'],
        finest_forward_level=finest,main_operational_reference='mesh2p5_timehalf',last_spatial_extension='PASS' if a.spatial_extension and all(v['status']=='PASS' for v in refinements if v['label']=='final_spatial_at_fine_time') else 'NOT_RUN' if not a.spatial_extension else 'UNRESOLVED',
        temporal_refinement_grid_nm=2.5,independent_temporal_refinement_at_1p25='NOT_RUN' if a.spatial_extension else 'NOT_APPLICABLE',
        numerical_origin_sensitivity=origin['status'],centered_force_admission=load(root/'centered_force_admission/result.json')['status'] if centered else 'NOT_RUN',
        scope='Nonperiodic direct elastic force, one shared CTB, two affine single-slip ordinary loops, I0 fail-closed reactions. Not general multisystem/junction evolution.',
        unsupported_cross_material_reactions='REJECT_AND_STOP; native negative collision preserved',
        exact_cumulative_restart='NOT_QUALIFIED',periodic_layers='NOT_RUN',VCoNi='NOT_QUALIFIED',strengthening_or_thickness_trend='NOT_QUALIFIED'))
    save('case_audits.json',audits);save('paired_initial.json',paired);save('refinement.json',refinements);save('paired_contrasts.json',contrasts)
    save('causal_resolution.json',causal);save('reverse_release.json',releases);save('ctb_ablation.json',ctb);save('metrics.json',metrics)
    save('regional_plastic.json',regional)
    save('resource_history.json',dict(aliases=aliases,historical_censored=load(root/'cases/mesh2p5_MT/result.json') if aliases else None,
        complete_logical_cases=len(audits),note='CENSORED remains historical; effective cases start from the original seed under separate IDs. No ledger continuation.'))
    with (out/'interface_reaction_ledger.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(reaction_ledger[0]));writer.writeheader();writer.writerows(reaction_ledger)
    with (out/'case_metrics.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['case_id',*next(iter(metrics.values()))]);writer.writeheader()
        for name,v in metrics.items():writer.writerow(dict(case_id=name,**v))
    print('G1_REVIEW',status,'cases',len(audits),'causal',causal)


if __name__=='__main__':main()
