"""Gate-controlled real ExaDiS microbenchmarks with accepted-step evidence."""
import argparse,csv,json,os,sys,time,traceback,shutil
from pathlib import Path
import numpy as np
from configuration import load,save_json
from native import imports,Simulation
from geometry import displacements
from events import Lineage
from contact import Contact
from diagnostics import Diagnostics,Censored


def area_tensor(a):
    if not len(a['positions_nm']):return np.zeros((3,3))
    p=a['positions_nm'];p=p-p.mean(axis=0);ids=a['ids']
    area=.5*np.cross(p[ids[:,0]],p[ids[:,1]])
    return np.einsum('si,sj->ij',a['burgers'],area)


def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',required=True);p.add_argument('--config',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--a0',type=Path,required=True)
    p.add_argument('--initial',type=Path);a=p.parse_args()
    gate=json.loads(a.a0.read_text());assert gate['gate']=='A0' and gate['status']=='PASS','A0 must PASS first'
    for name in ('configuration.py','geometry.py','native.py'):
        assert (a.a0.parent/'production_snapshot'/name).read_bytes()==(Path(__file__).resolve().parent/name).read_bytes(),f'A0 production snapshot changed: {name}; rerun A0 before physics'
    c=load(a.config);out=a.out.resolve()
    if out.exists() and any(out.iterdir()):raise RuntimeError('Refusing nonempty output')
    out.mkdir(parents=True);save_json(out/'config.json',c)
    (out/'code_snapshot').mkdir()
    for source in Path(__file__).resolve().parent.glob('*.py'):shutil.copy2(source,out/'code_snapshot'/source.name)
    save_json(out/'command.json',dict(argv=sys.argv,python=sys.executable,cwd=os.getcwd(),config_source=str(a.config.resolve())))
    os.environ['OMP_NUM_THREADS']=str(c['threads']);os.environ['OMP_PROC_BIND']='false'
    py,exa=imports(a.repo);py.initialize(num_threads=c['threads']);start=time.monotonic();sim=None;step=0;now=0.
    try:
        restart=None
        if a.initial:
            raw=dict(np.load(a.initial));bnm=c['b_m']*1e9
            restart=dict(cell=dict(h=raw['cell_nm']/bnm,origin=np.zeros(3),is_periodic=c['periodic']),
                nodes=dict(tags=raw['tags'],positions=raw['positions_nm']/bnm,constraints=raw['constraints']),
                segs=dict(nodeids=raw['ids'],burgers=raw['burgers'],planes=raw['planes']))
        sim=Simulation(c,a.repo,initial=restart);save_json(out/'sources.json',sim.sources)
        if c['walls']:
            sim.contact=Contact(c);sim.contact.validate(sim,'initialization')
        save_json(out/'initial_segment_metadata.json',dict(roles=sim.initial_roles,source_ids=sim.initial_source_ids))
        init=sim.arrays();np.savez_compressed(out/'initial.npz',**init);sim.g.write_data(str(out/'initial.data'))
        tracker=Lineage(sim,out);diag=Diagnostics(out,c);last_legal=init
        initial_area=area_tensor(init);initial_measures=sim.measures();plastic=np.zeros((3,3));vol=np.prod(c['box_nm'])
        pins={tuple(t):x for t,x,k in zip(init['tags'],init['positions_nm'],init['constraints'].ravel()) if k==7}
        last_frame=-1.;max_res=0.;max_pin_error=0.;min_dt=1.;max_closure=0.;sum_sweep=0.;history=[]
        with (out/'history.csv').open('w',newline='',encoding='utf-8') as stream:
            writer=None
            while True:
                if sim.g.num_nodes()==0:
                    print('NETWORK_ANNIHILATED',step,flush=True);break
                f,v=sim.evaluate();current=sim.arrays(False);meas=sim.measures();now=sim.state.get('time',0.)
                for axis,periodic in enumerate(c['periodic']):
                    if not periodic and (np.any(current['positions_nm'][:,axis]<0) or np.any(current['positions_nm'][:,axis]>c['box_nm'][axis])):
                        raise RuntimeError('Network escaped nonperiodic computational box; enlarge domain before qualification')
                max_res=max(max_res,meas['max_burgers_residual_b'])
                current_tags={tuple(t):i for i,t in enumerate(current['tags'])};lost=0
                for tag,pos in pins.items():
                    if tag not in current_tags:lost+=1;continue
                    error=current['positions_nm'][current_tags[tag]]-pos
                    error-=np.round(error/np.array(c['box_nm']))*np.array(c['box_nm'])*np.array(c['periodic'])
                    max_pin_error=max(max_pin_error,float(np.linalg.norm(error)))
                if max_res>1e-8 or max_pin_error>1e-4:raise RuntimeError('Burgers or fixed anchor violation')
                if lost:raise RuntimeError(f'{lost} original fixed source anchors lost')
                global_geom=-sim.bnm*(area_tensor(current)-initial_area)/vol;eps_geom=.5*(global_geom+global_geom.T)
                closure=float(np.linalg.norm(plastic-eps_geom));max_closure=max(max_closure,closure)
                if not any(c['periodic']) and closure>1e-12+1e-8*np.linalg.norm(eps_geom):
                    raise RuntimeError('Native plastic strain / geometric swept-area closure violation')
                row=dict(step=step,time_ns=now*1e9,**meas,confirmed_emissions=tracker.confirmed,
                         force_max_N=float(np.linalg.norm(f,axis=1).max())*sim.b**2,
                         velocity_max_m_s=float(np.linalg.norm(v,axis=1).max())*sim.b,
                         pin_error_nm=max_pin_error,lost_anchors=lost,
                         swept_area_nm2=sum_sweep,plastic_area_closure_absolute=closure,
                         wall_seconds=time.monotonic()-start)
                if sim.contact is not None:row.update(sim.contact.measures(sim))
                if writer is None:writer=csv.DictWriter(stream,fieldnames=list(row));writer.writeheader()
                writer.writerow(row)
                if step%20==0:stream.flush()
                if now*1e9-last_frame>=c['frame_interval_ns']-1e-12 or last_frame<0:
                    np.savez_compressed(out/f'frame_{step:07d}.npz',**sim.arrays())
                    history.append(dict(step=step,time_ns=now*1e9));last_frame=now*1e9
                    print('STEP',step,'T_NS',round(now*1e9,6),'N',meas['nodes'],'EMIT',tracker.confirmed,'WALL_S',round(time.monotonic()-start,2),flush=True)
                if c['stop_emissions'] and tracker.confirmed>=c['stop_emissions']:break
                if now>=c['end_ns']*1e-9-1e-24:break
                if step>=c['max_steps']:raise Censored('Accepted-step resource bound')
                if time.monotonic()-start>c['max_wall_s']:raise Censored('Wall-time resource bound')
                if meas['nodes']>c['max_nodes'] or meas['segments']>c['max_segments']:raise Censored('Network-size resource bound')
                old_fields={k:v.copy() for k,v in sim.velocity_diagnostic.items()};old_force_N=f.copy()*sim.b**2
                old=current;old_native=sim.g.get_nodes_data();sim.state['oldnodes_dict']=old_native
                remaining=c['end_ns']*1e-9-now
                if c['integrator']=='EulerForward':
                    dt=min(c['maxdt_ps']*1e-12,remaining,c['max_node_step_nm']/max(np.linalg.norm(v,axis=1).max()*sim.bnm,1e-99))
                    if sim.contact is not None:dt=min(dt,sim.contact.event_dt(sim,v))
                    sim.advance(dt)
                else:
                    if remaining<c['maxdt_ps']*1e-12:
                        state=dict(sim.state,maxdt=remaining,nextdt=remaining)
                        sim.integrator=exa.TimeIntegration(state=state,integrator='Trapezoid',force=sim.fadapter,mobility=sim.mobility)
                    sim.advance()
                dt=sim.state['dt'];min_dt=min(min_dt,dt);step+=1
                if dt<1e-20:raise RuntimeError('Native timestep collapse')
                if sim.contact is not None:sim.contact.validate(sim,'integration_substep')
                last_legal=sim.arrays(False)
                moved=sim.arrays(False);delta=moved['positions_nm']-old['positions_nm']
                delta-=np.round(delta/np.array(c['box_nm']))*np.array(c['box_nm'])*np.array(c['periodic'])
                dr=displacements(old['positions_nm'],old['ids'],c['box_nm'],c['periodic'])
                area=.5*np.cross(dr,delta[old['ids'][:,0]]+delta[old['ids'][:,1]])+.5*np.cross(delta[old['ids'][:,1]],delta[old['ids'][:,0]])
                sum_sweep+=float(np.linalg.norm(area,axis=1).sum())
                previous_area=area_tensor(moved);diag.area('integration',previous_area-area_tensor(old))
                tracker.update(sim,'integration',step,sim.state['time']*1e9)
                if c['collision']:
                    sim.collision.HandleCol(sim.net,sim.state);tracker.update(sim,'collision',step,sim.state['time']*1e9)
                if sim.contact is not None:sim.contact.validate(sim,'collision')
                changed=sim.arrays(False);new_area=area_tensor(changed);diag.area('collision',new_area-previous_area);previous_area=new_area;last_legal=changed
                sim.force.PreCompute(sim.net,sim.state)
                if c['topology']:
                    sim.topology.Handle(sim.net,sim.state);tracker.update(sim,'topology',step,sim.state['time']*1e9)
                if sim.contact is not None:sim.contact.validate(sim,'topology')
                changed=sim.arrays(False);new_area=area_tensor(changed);diag.area('topology',new_area-previous_area);previous_area=new_area;last_legal=changed
                if c['remesh']:
                    sim.remesh.Remesh(sim.net,sim.state);tracker.update(sim,'remesh',step,sim.state['time']*1e9)
                if sim.contact is not None:sim.contact.validate(sim,'remesh')
                changed=sim.arrays(False);new_area=area_tensor(changed);diag.area('remesh',new_area-previous_area);last_legal=changed
                dep,_,_=sim.g.net.get_plastic_strain();plastic+=np.array(dep).reshape(3,3)
                diag.external_work_J+=float(np.einsum('ij,ij->',np.array(c['stress_pa']),np.array(dep).reshape(3,3))*vol*1e-27)
                diag.append(step,now,dt,old,old_fields,old_force_N)
                tracker.accepted(step,sim.state['time']*1e9)
        final=sim.arrays(sim.g.num_nodes()>0);np.savez_compressed(out/'final.npz',**final);sim.g.write_data(str(out/'final.data'));tracker.save()
        save_json(out/'frames.json',history)
        diag.save()
        save_json(out/'contact_events.json',sim.contact.events if sim.contact else [])
        source_status='PASS' if c['kind']=='FR_source' and tracker.confirmed>=3 else ('NOT_RUN' if c['stop_emissions']==1 else 'FAIL') if c['kind']=='FR_source' else 'NOT_APPLICABLE'
        save_json(out/'result.json',dict(status='PASS',gate=c['gate'],source_qualification=source_status,steps=step,time_ns=sim.state.get('time',0.)*1e9,
            confirmed_emissions=tracker.confirmed,max_burgers_residual_b=max_res,pin_error_nm=max_pin_error,
            first_emission_observed=tracker.confirmed>=1,target_emissions=c['stop_emissions'],
            min_dt_s=min_dt,plastic_strain=plastic,swept_area_nm2=sum_sweep,
            max_plastic_area_closure_absolute=max_closure,initial=initial_measures,
            final=sim.measures(),wall_seconds=time.monotonic()-start,material_qualification='NOT_QUALIFIED',
            contact_checks=sim.contact.checks if sim.contact else {},
            minimum_contact_gap_nm=sim.contact.minimum_gap_nm if sim.contact and sim.contact.walls else None,
            stop_reason='network_annihilated' if sim.g.num_nodes()==0 else 'qualified_emissions' if c['stop_emissions'] and tracker.confirmed>=c['stop_emissions'] else 'end_time'))
        print('NUMERICAL_PASS SOURCE_3_EMISSION',source_status,'STEPS',step,'T_NS',sim.state.get('time',0.)*1e9,'EMISSIONS',tracker.confirmed,flush=True)
    except Censored as exc:
        np.savez_compressed(out/'checkpoint.npz',**sim.arrays(False),cumulative_plastic_strain=plastic,
            time_s=sim.state.get('time',0.),accepted_step=step,cumulative_swept_area_nm2=sum_sweep)
        tracker.save();diag.save();save_json(out/'frames.json',history)
        save_json(out/'contact_events.json',sim.contact.events if sim.contact else [])
        save_json(out/'result.json',dict(status='CENSORED',gate='W4-G0',error=str(exc),steps=step,
            time_ns=sim.state.get('time',0.)*1e9,wall_seconds=time.monotonic()-start,
            checkpoint='checkpoint.npz',restart_status='NOT_RUN: exact cumulative-ledger continuation must be qualified before combining segments',material_qualification='NOT_QUALIFIED'))
        print('CENSORED',str(exc),flush=True)
    except Exception as exc:
        if sim is not None:np.savez_compressed(out/'failure.npz',**sim.arrays(False))
        if 'last_legal' in locals():np.savez_compressed(out/'last_legal.npz',**last_legal)
        if 'diag' in locals():diag.save()
        if sim is not None:save_json(out/'contact_events.json',sim.contact.events if sim.contact else [])
        if 'tracker' in locals():tracker.save()
        save_json(out/'result.json',dict(status='FAIL',gate=c['gate'],error=repr(exc),steps=step,
             time_ns=sim.state.get('time',0.)*1e9 if sim is not None else 0.,wall_seconds=time.monotonic()-start,
             material_qualification='NOT_QUALIFIED'))
        traceback.print_exc();raise
    finally:py.finalize()


if __name__=='__main__':main()
