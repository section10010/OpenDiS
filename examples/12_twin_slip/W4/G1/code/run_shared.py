"""Accepted-step real ExaDiS shared-CTB evolution; G0 source remains immutable."""
import argparse,copy,csv,json,os,sys,time,shutil,traceback
from pathlib import Path
import numpy as np
from material_identity import MaterialIdentity,network_arrays
from shared_contact import SharedContact
from shared_force_probe import initial_data
from specification import validate_spec


def bootstrap(g0):
    gate=json.loads((g0/'analysis_r5/gate_review.json').read_text())
    assert gate['status']=='PASS' and gate['G1_admission']=='PASS'
    assert json.loads((g0/'A0/gate_result.json').read_text())['status']=='PASS'
    for name in ('configuration.py','geometry.py','native.py'):
        assert (g0/'code'/name).read_bytes()==(g0/'A0/production_snapshot'/name).read_bytes()
    sys.path.insert(0,str(g0/'code'))


def seed_from_npz(path,c):
    with np.load(path) as a:
        bnm=c['b_m']*1e9
        return dict(cell=dict(h=a['cell_nm']/bnm,origin=np.zeros(3),is_periodic=c['periodic']),
            nodes=dict(tags=a['tags'],positions=a['positions_nm']/bnm,constraints=a['constraints']),
            segs=dict(nodeids=a['ids'],burgers=a['burgers'],planes=a['planes']))


def tensors(a,identity):
    meta=identity.inspect(a);ij=a['ids'];p=a['positions_nm']-np.array([2000.,2000.,2000.])
    area=.5*np.cross(p[ij[:,0]],p[ij[:,1]])
    return {h:np.einsum('si,sj->ij',a['burgers'][meta['segment_host']==h],area[meta['segment_host']==h]) for h in ('M','T')}


def snapshot(sim,identity,path):
    a=sim.arrays();np.savez_compressed(path,**a,**identity.fields(a));return a


def guarded_stage(sim,name,operation,out,step):
    """Reject invalid native proposals, retaining both sides of the transaction.

    No automatic retries or silent skip. A rejected production proposal ends
    that case. Rollback is demonstrable graph/state rollback for safe recovery,
    not an exact cumulative-ledger restart claim (which remains unqualified).
    """
    from native import Simulation
    from configuration import save_json
    before=copy.deepcopy(sim.snapshot());before_state=copy.deepcopy(sim.state)
    before_identity=copy.deepcopy(sim.contact.identity);before_contact=copy.deepcopy(sim.contact)
    try:
        operation();sim.contact.validate(sim,name)
    except Exception as exc:
        target=Path(out)/f'rejected_{step:07d}_{name}';target.mkdir(parents=True,exist_ok=False)
        a=network_arrays(sim);np.savez_compressed(target/'first_invalid.npz',**a)
        restored=Simulation(sim.c,sim.repo,initial=before);restored.repo=sim.repo
        restored.state.update(before_state);restored.contact=before_contact;restored.contact.identity=before_identity
        legal=network_arrays(restored);np.savez_compressed(target/'restored_legal.npz',**legal,**before_identity.fields(legal))
        restored.contact.validate(restored,'rollback_verify')
        expected=before['nodes'];np.testing.assert_array_equal(legal['tags'],expected['tags'])
        np.testing.assert_array_equal(legal['positions_nm'],expected['positions']*sim.bnm)
        save_json(target/'rejection.json',dict(status='FAIL',stage=name,step=step,error=repr(exc),
            rollback_graph='PASS',rollback_time='PASS' if restored.state.get('time',0)==before_state.get('time',0) else 'FAIL',
            production_action='STOP; no skipped physical operation',exact_cumulative_restart='NOT_QUALIFIED'))
        raise RuntimeError('Native stage rejected and graph restored: '+str(target)) from exc


def run_case(spec,population,repo,out,initial=None,source_metadata=None,negative_penetration=False):
    from configuration import config,save_json
    from geometry import build,residual
    from native import Simulation
    from diagnostics import Diagnostics,Censored
    out=Path(out)
    if out.exists():raise RuntimeError('Refuse overwrite: '+str(out))
    out.mkdir(parents=True)
    c=config(**spec['common_native_updates'])
    validate_spec(spec,c)
    geometry_config=config(**{**c,**spec.get('initial_geometry_updates',{})})
    seed,_,sources=initial_data(geometry_config,population,build,spec)
    if initial is not None:seed=seed_from_npz(initial,c);sources=source_metadata
    save_json(out/'config.json',c);save_json(out/'spec.json',spec);save_json(out/'sources.json',sources)
    save_json(out/'command.json',dict(argv=sys.argv,python=sys.executable,cwd=os.getcwd(),population=population,initial=None if initial is None else str(initial)))
    code=out/'code_snapshot';code.mkdir()
    for p in Path(__file__).resolve().parent.glob('*.py'):shutil.copy2(p,code/p.name)
    identity=MaterialIdentity(sources,spec['interface'],spec['criteria']['source_plane_tol_nm'])
    sim=Simulation(c,repo,initial=seed);sim.repo=repo
    if 'force_evaluation_origin_nm' in spec:
        from centered_force import install
        install(sim,spec['force_evaluation_origin_nm'])
    sim.contact=SharedContact(identity,c['contact_tol_nm'],spec.get('contact_enabled',True));sim.contact.validate(sim,'initialization')
    diag=Diagnostics(out,c);start=time.monotonic();step=0;plastic=np.zeros((3,3));sweep={h:0. for h in ('M','T')}
    initial_a=snapshot(sim,identity,out/'initial.npz');sim.g.write_data(str(out/'initial.data'))
    base=tensors(initial_a,identity);stage_tensors={};max_closure=0.;max_normal=0.;max_res=0.;max_glide=0.
    first_contact={h:None for h in population};released={h:None for h in population};initial_contact=sim.contact.metrics(sim)
    last_frame=-1.;status='FAIL';error=None;rows=[];max_position_step=0.
    try:
        while sim.g.num_nodes():
            f,v=sim.evaluate();a=sim.arrays(False);now=sim.state.get('time',0.);meta=identity.inspect(a)
            res=float(residual(len(a['tags']),a['ids'],a['burgers']).max(initial=0));max_res=max(max_res,res)
            if res>1e-8:raise RuntimeError('Burgers conservation failed')
            if np.any(a['positions_nm']<0) or np.any(a['positions_nm']>np.array(c['box_nm'])):raise RuntimeError('Escaped nonperiodic box')
            gap_sign=np.where(meta['node_host']=='M',-1.,1.);gap=gap_sign*(a['positions_nm'][:,2]-spec['interface']['z_nm'])
            fields={k:x.copy() for k,x in sim.velocity_diagnostic.items()}
            scale=np.maximum(np.linalg.norm(fields['raw_velocity_m_s'],axis=1),np.linalg.norm(fields['constrained_velocity_m_s'],axis=1))
            normal=np.maximum(0.,-gap_sign*fields['constrained_velocity_m_s'][:,2])/np.maximum(scale,1e-30)
            if spec.get('contact_enabled',True):max_normal=max(max_normal,float(normal[gap<=c['contact_tol_nm']].max(initial=0)))
            if max_normal>1e-8:raise RuntimeError('Shared CTB normal complementarity failed')
            for e,(i,j) in enumerate(a['ids']):
                max_glide=max(max_glide,abs(float(a['planes'][e]@fields['constrained_velocity_m_s'][i])),abs(float(a['planes'][e]@fields['constrained_velocity_m_s'][j])))
            t=tensors(a,identity);geom=-sim.bnm*sum((t[h]-base[h] for h in ('M','T')),np.zeros((3,3)))/np.prod(c['box_nm'])
            expected=.5*(geom+geom.T);closure=float(np.linalg.norm(plastic-expected));max_closure=max(max_closure,closure)
            if closure>1e-12+1e-8*np.linalg.norm(expected):raise RuntimeError('Native plastic/global affine-area closure failed')
            metrics=sim.contact.metrics(sim)
            for h in population:
                if metrics[h+'_contact_nodes'] and first_contact[h] is None:first_contact[h]=now*1e9
                if initial_contact[h+'_contact_nodes'] and released[h] is None and metrics[h+'_gap_nm']>10*c['contact_tol_nm']:released[h]=now*1e9
            row=dict(step=step,time_ns=now*1e9,nodes=len(a['tags']),segments=len(a['ids']),plastic_closure=closure,
                max_force_N=float(np.linalg.norm(f,axis=1).max(initial=0))*sim.b**2,max_velocity_m_s=float(np.linalg.norm(v,axis=1).max(initial=0))*sim.b,
                **metrics,**{h+'_swept_area_nm2':sweep[h] for h in ('M','T')})
            for h in ('M','T'):
                # Tensor-valued signed plastic remains available; never equate
                # stored line length with plastic strain or invent scalar M/T RSS.
                row[h+'_area_tensor_norm_nm2']=float(np.linalg.norm(t[h]-base[h]))
            rows.append(row)
            if now*1e9-last_frame>=c['frame_interval_ns']-1e-12 or last_frame<0:
                snapshot(sim,identity,out/f'frame_{step:07d}.npz');last_frame=now*1e9
                print('STEP',population,step,round(now*1e9,6),'N',len(a['tags']),'CONTACT',[metrics[h+'_contact_nodes'] for h in population],flush=True)
            if now>=c['end_ns']*1e-9-1e-24:break
            if step>=c['max_steps'] or time.monotonic()-start>c['max_wall_s']:raise Censored('Frozen G1 time/step budget')
            if len(a['tags'])>c['max_nodes'] or len(a['ids'])>c['max_segments']:raise Censored('Frozen G1 network budget')
            dt=min(c['maxdt_ps']*1e-12,c['end_ns']*1e-9-now,c['max_node_step_nm']/max(np.linalg.norm(v,axis=1).max()*sim.bnm,1e-99),sim.contact.event_dt(sim,v))
            if negative_penetration:dt=20e-12
            if dt<1e-20:raise RuntimeError('Timestep collapse')
            before=a;previous=t
            guarded_stage(sim,'integration',lambda:sim.advance(dt),out,step)
            moved=sim.arrays(False);delta=moved['positions_nm']-before['positions_nm'];max_position_step=max(max_position_step,float(np.linalg.norm(delta,axis=1).max(initial=0)))
            ij=before['ids'];dr=before['positions_nm'][ij[:,1]]-before['positions_nm'][ij[:,0]]
            area=.5*np.cross(dr,delta[ij[:,0]]+delta[ij[:,1]])+.5*np.cross(delta[ij[:,1]],delta[ij[:,0]])
            for h in ('M','T'):sweep[h]+=float(np.linalg.norm(area[meta['segment_host']==h],axis=1).sum())
            for name,operation,enabled in [('integration_account',lambda:None,True),('collision',lambda:sim.collision.HandleCol(sim.net,sim.state),c['collision']),
                    ('topology',lambda:(sim.force.PreCompute(sim.net,sim.state),sim.topology.Handle(sim.net,sim.state)),c['topology']),
                    ('remesh',lambda:sim.remesh.Remesh(sim.net,sim.state),c['remesh'])]:
                if name!='integration_account' and enabled:guarded_stage(sim,name,operation,out,step)
                current=tensors(sim.arrays(False),identity)
                for h in ('M','T'):stage_tensors[(name,h)]=stage_tensors.get((name,h),np.zeros((3,3)))+current[h]-previous[h]
                diag.area(name,sum((current[h]-previous[h] for h in ('M','T')),np.zeros((3,3))));previous=current
            dep,_,_=sim.g.net.get_plastic_strain();plastic+=np.asarray(dep).reshape(3,3)
            diag.external_work_J+=float(np.einsum('ij,ij->',np.array(c['stress_pa']),np.array(dep))*np.prod(c['box_nm'])*1e-27)
            fields.update(node_host=meta['node_host'],node_source_id=meta['node_source_id'])
            diag.append(step+1,now,sim.state['dt'],before,fields,f*sim.b**2);step+=1
        final=tensors(sim.arrays(False),identity)
        geom=-sim.bnm*sum((final[h]-base[h] for h in ('M','T')),np.zeros((3,3)))/np.prod(c['box_nm'])
        closure=float(np.linalg.norm(plastic-.5*(geom+geom.T)));max_closure=max(max_closure,closure)
        if closure>1e-12+1e-8*np.linalg.norm(geom):raise RuntimeError('Final/empty native plastic closure failed')
        status='PASS'
    except Censored as exc:status='CENSORED';error=repr(exc)
    except Exception as exc:error=repr(exc);(out/'traceback.log').write_text(traceback.format_exc(),encoding='utf-8')
    finally:
        if rows:
            with (out/'history.csv').open('w',newline='',encoding='utf-8') as f:
                writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        diagnostic=diag.save()
        save_json(out/'identity_events.json',identity.events);save_json(out/'contact_events.json',sim.contact.events)
        save_json(out/'stage_area_tensors.json',{name+'_'+h:v for (name,h),v in stage_tensors.items()})
        save_json(out/'plastic.json',dict(native_cumulative=plastic,initial_tensors=base,source='Actual ExaDiS get_plastic_strain after accepted native stages'))
        try:
            snapshot(sim,identity,out/'final.npz');sim.g.write_data(str(out/'final.data'))
        except Exception:np.savez_compressed(out/'final_unqualified.npz',**network_arrays(sim))
        result=dict(status=status,error=error,population=population,steps=step,time_ns=sim.state.get('time',0.)*1e9,wall_s=time.monotonic()-start,
            first_contact_ns=first_contact,first_full_release_ns=released,minimum_gap_nm=sim.contact.minimum_gap_nm if np.isfinite(sim.contact.minimum_gap_nm) else None,
            max_burgers_residual_b=max_res,max_plastic_closure=max_closure,max_normalized_contact_violation=max_normal,
            max_glide_normal_velocity_m_s=max_glide,max_node_displacement_nm=max_position_step,swept_area_nm2=sweep,
            stage_checks=identity.stage_counts,diagnostics=diagnostic,contact_enabled=spec.get('contact_enabled',True),
            qualification='Single affine-slip source per side, ideal I0; unsupported junctions fail closed. No VCoNi calibration or CTB transmission.')
        save_json(out/'result.json',result);print('G1_CASE',status,population,step,error,flush=True)
    if status!='PASS':raise RuntimeError('G1 case '+status+': '+str(out))
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--g0-root',type=Path,required=True);p.add_argument('--repo',required=True)
    p.add_argument('--spec',type=Path,required=True);p.add_argument('--population',choices=['M','T','MT'],required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--initial',type=Path);p.add_argument('--sources',type=Path)
    p.add_argument('--negative-penetration',action='store_true');a=p.parse_args();bootstrap(a.g0_root)
    from native import imports
    spec=json.loads(a.spec.read_text());assert spec['schema']=='W4.G1.SharedEvolution.v1'
    os.environ['OMP_NUM_THREADS']=str(spec['common_native_updates']['threads']);os.environ['OPENBLAS_NUM_THREADS']='1'
    py,_=imports(a.repo);py.initialize(num_threads=spec['common_native_updates']['threads'])
    try:run_case(spec,a.population,a.repo,a.out,a.initial,None if a.sources is None else json.loads(a.sources.read_text()),a.negative_penetration)
    finally:py.finalize()


if __name__=='__main__':main()
