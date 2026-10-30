"""Accepted-step raw diagnostics; no invented velocities or energy field.

Ragged native node arrays are concatenated per chunk with step offsets. A sample
is the actual left state and its accepted dt, never an interpolated saved frame.
"""
from pathlib import Path
import numpy as np
from configuration import save_json


class Censored(RuntimeError):pass


def node_weights(a):
    ij=a['ids'];w=np.zeros(len(a['positions_nm']))
    length=np.linalg.norm(a['positions_nm'][ij[:,1]]-a['positions_nm'][ij[:,0]],axis=1)
    for end in (0,1):np.add.at(w,ij[:,end],length/2)
    return w


class Diagnostics:
    def __init__(self,out,c):
        self.root=Path(out)/'accepted_velocity';self.root.mkdir();self.c=c;self.buffer=[];self.chunks=[]
        self.weight_time=0.;self.cap_weight_time=0.;self.capped_steps=0;self.steps=0
        self.max_raw=0.;self.max_constrained=0.;self.max_power_residual=0.;self.max_normal_residual=0.
        self.reaction_work_J=0.;self.bulk_dissipation_left_J=0.;self.external_work_J=0.;self.stage_area={}
        self.integration_area_tensor=np.zeros((3,3));self.stage_area_tensors={}

    def append(self,step,time_s,dt,a,fields,force_N):
        w=node_weights(a);pre=np.linalg.norm(fields['constrained_velocity_m_s'],axis=1)
        cap=self.c['vmax_m_s'];hit=np.zeros(len(w),bool) if cap is None else pre>cap
        self.weight_time+=float(w.sum()*dt);self.cap_weight_time+=float(w@hit*dt)
        self.capped_steps+=bool(hit.any());self.steps+=1
        self.max_raw=max(self.max_raw,float(np.linalg.norm(fields['raw_velocity_m_s'],axis=1).max(initial=0)))
        self.max_constrained=max(self.max_constrained,float(pre.max(initial=0)))
        reaction=fields['reaction_N'];vc=fields['constrained_velocity_m_s'];v=fields['capped_velocity_m_s']
        velocity_scale=np.maximum(np.linalg.norm(fields['raw_velocity_m_s'],axis=1),pre)
        scale=np.maximum(np.linalg.norm(force_N,axis=1),np.linalg.norm(reaction,axis=1))*velocity_scale
        power=np.abs(np.einsum('ij,ij->i',reaction,vc))/np.maximum(scale,1e-30)
        self.max_power_residual=max(self.max_power_residual,float(power.max(initial=0)))
        self.reaction_work_J+=float(np.einsum('ij,ij->',reaction,v)*dt)
        self.bulk_dissipation_left_J+=float(np.einsum('ij,ij->',force_N+reaction,v)*dt)
        if self.max_power_residual>1e-8:raise RuntimeError('Normalized contact power exceeds frozen G0 tolerance')
        if np.any(np.einsum('ij,ij->i',force_N+reaction,v)<-1e-18):raise RuntimeError('Negative native mobility dissipation')
        for wall in self.c['walls']:
            if not wall['enabled']:continue
            gap=wall['inward']*(a['positions_nm'][:,2]-wall['z_nm']);active=gap<=self.c['contact_tol_nm']
            err=np.maximum(0,-wall['inward']*vc[:,2])/np.maximum(velocity_scale,1e-30)
            self.max_normal_residual=max(self.max_normal_residual,float(err[active].max(initial=0)))
        if self.max_normal_residual>1e-8:raise RuntimeError('Contact normal residual exceeds frozen G0 tolerance')
        self.buffer.append(dict(step=step,time_s=time_s,dt_s=dt,tags=a['tags'],positions_nm=a['positions_nm'],
            node_line_weight_nm=w,force_N=force_N,**fields))
        if len(self.buffer)>=self.c['diagnostic_chunk_steps']:self.flush()

    def area(self,stage,delta_tensor):
        self.stage_area_tensors[stage]=self.stage_area_tensors.get(stage,np.zeros((3,3)))+delta_tensor
        self.stage_area[stage]=self.stage_area.get(stage,0.)+float(np.linalg.norm(delta_tensor))

    def flush(self):
        if not self.buffer:return
        offsets=np.cumsum([0]+[len(x['tags']) for x in self.buffer]);packed=dict(offsets=offsets)
        for k in self.buffer[0]:
            packed[k]=np.array([x[k] for x in self.buffer]) if k in ('step','time_s','dt_s') else np.concatenate([x[k] for x in self.buffer])
        name=f'chunk_{len(self.chunks):05d}.npz';np.savez_compressed(self.root/name,**packed)
        self.chunks.append(name);self.buffer=[]

    def save(self):
        self.flush()
        summary=dict(status='PASS' if self.steps else 'NOT_RUN',accepted_samples=self.steps,chunks=self.chunks,
            node_line_time_nm_s=self.weight_time,capped_node_line_time_nm_s=self.cap_weight_time,
            line_time_cap_fraction=self.cap_weight_time/self.weight_time if self.weight_time else None,
            any_node_capped_record_fraction=self.capped_steps/self.steps if self.steps else None,
            maximum_raw_velocity_m_s=self.max_raw,maximum_constrained_velocity_m_s=self.max_constrained,
            max_normalized_contact_power=self.max_power_residual,max_normalized_contact_normal_violation=self.max_normal_residual,
            contact_reaction_work_J=self.reaction_work_J,left_rule_bulk_dissipation_J=self.bulk_dissipation_left_J,
            external_work_from_native_plastic_J=self.external_work_J,
            stage_area_norm_sum_nm2=self.stage_area,stage_area_tensors_nm2=self.stage_area_tensors,
            elastic_potential_energy_status='NOT_RUN: current Python native force binding exposes forces, not elastic energy; work is not relabelled as an independent energy',
            convention='Three exact native left-state velocities, with actual accepted dt and half-incident-line weights; not frame-count or interpolated estimates')
        save_json(self.root.parent/'velocity_diagnostics.json',summary)
        return summary
