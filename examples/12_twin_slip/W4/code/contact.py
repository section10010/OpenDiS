"""Ideal unilateral CTB contact on a prescribed glide manifold (I0 benchmark).

No node pinning or post-step coordinate clipping. Real native GLIDE supplies
both v=Mf and M times the wall normal. Every accepted native operation is checked.
"""
import numpy as np


class Contact:
    def __init__(self,c):
        self.walls=[w for w in c['walls'] if w['enabled']]
        self.tol=c['contact_tol_nm'];self.checks={};self.minimum_gap_nm=float('inf')
        self.events=[];self.active_tags={w['interface_id']:set() for w in self.walls}

    def project_one(self,sim,N,state,tag,v,positions=None,native_index=None,diagnostic=None):
        if positions is None:
            d=N.get_disnet(sim.exa.ExaDisNet).get_nodes_data()
            tags=d['tags'];ix=np.flatnonzero(np.all(tags==tag,axis=1))
            if not len(ix):raise RuntimeError('Contact tag absent')
            z=float(d['positions'][ix[0],2])*sim.bnm
        else:z=positions
        v=np.array(v,float)
        for j,wall in enumerate(self.walls):
            sign=wall['inward'];gap=sign*(z-wall['z_nm'])
            if gap<=self.tol and sign*v[2]<0:
                normal=np.array([0.,0.,float(sign)])
                if native_index is None:mn=np.array(sim.glide.OneNodeMobility(N,state,tag,normal,False))
                else:
                    g=N.get_disnet(sim.exa.ExaDisNet)
                    mn=np.array(sim.py.compute_node_mobility(g.net,native_index,mobility=sim.glide.mobility,fi=normal))
                denom=normal@mn
                if denom>0:
                    multiplier=-(normal@v)/denom
                    v=v+mn*multiplier
                    if diagnostic is not None:
                        reaction,multipliers=diagnostic
                        reaction[native_index]+=multiplier*normal
                        multipliers[native_index,j]=multiplier
                elif abs(normal@v)>1e-8:raise RuntimeError('Singular contact mobility')
        return v

    def project_all(self,sim,N,state,velocities):
        d=N.get_disnet(sim.exa.ExaDisNet).get_nodes_data();v=np.array(velocities,copy=True)
        self.reaction_native=np.zeros_like(v);self.multipliers_native=np.zeros((len(v),len(self.walls)))
        # Mobility output follows native node ordering; verify against nodeveltags.
        if 'nodeveltags' in state and not np.array_equal(state['nodeveltags'],d['tags']):raise RuntimeError('Velocity tag order mismatch')
        for i,(tag,x) in enumerate(zip(d['tags'],d['positions'])):
            z=float(x[2])*sim.bnm
            if any(w['inward']*(z-w['z_nm'])<=self.tol and w['inward']*v[i,2]<0 for w in self.walls):
                v[i]=self.project_one(sim,N,state,tag,v[i],z,i,(self.reaction_native,self.multipliers_native))
        return v

    def event_dt(self,sim,v):
        x=sim.arrays(False)['positions_nm'];limit=float('inf')
        for w in self.walls:
            gap=w['inward']*(x[:,2]-w['z_nm']);speed=-w['inward']*v[:,2]*sim.bnm
            mask=(gap>self.tol)&(speed>0)
            if mask.any():limit=min(limit,float(np.min(gap[mask]/speed[mask])))
        return limit

    def validate(self,sim,stage):
        a=sim.arrays(False);x=a['positions_nm'];self.checks[stage]=self.checks.get(stage,0)+1
        for w in self.walls:
            gap=w['inward']*(x[:,2]-w['z_nm']);m=float(gap.min(initial=float('inf')))
            self.minimum_gap_nm=min(self.minimum_gap_nm,m)
            # With nonperiodic straight segments a convex half-space is also
            # segment-admissible iff its endpoints are admissible.
            if m < -self.tol:raise RuntimeError(f'Illegal CTB penetration after {stage}: {w["interface_id"]}, gap={m} nm')
            key=w['interface_id'];active={tuple(t) for t in a['tags'][gap<=self.tol]};previous=self.active_tags[key]
            if active!=previous:
                self.events.append(dict(event='contact_membership_change',stage=stage,time_ns=sim.state.get('time',0.)*1e9,
                    interface_id=key,entered_tags=sorted(active-previous),exited_tags=sorted(previous-active),
                    note='Remesh tag changes are membership changes, not independent physical impacts'))
                self.active_tags[key]=active

    def measures(self,sim):
        a=sim.arrays(False);ij=a['ids'];x=a['positions_nm'];l=np.linalg.norm(x[ij[:,1]]-x[ij[:,0]],axis=1)
        out={}
        for w in sim.c['walls']:
            gap=w['inward']*(x[:,2]-w['z_nm']);label=w['interface_id']
            active=np.abs(gap)<=self.tol
            out[label+'_gap_nm']=float(gap.min())
            out[label+'_contact_nodes']=int(active.sum()) if w['enabled'] else 0
            out[label+'_tail_nm']=float(l[np.all(active[ij],axis=1)].sum()) if w['enabled'] else 0.
            # Geometric segment-band intersection, not midpoint selection.
            z0=x[ij[:,0],2];z1=x[ij[:,1],2];dz=z1-z0
            lo=w['z_nm']-10.;hi=w['z_nm']+10.
            aa=np.divide(lo-z0,dz,out=np.zeros_like(dz),where=np.abs(dz)>1e-12)
            bb=np.divide(hi-z0,dz,out=np.zeros_like(dz),where=np.abs(dz)>1e-12)
            fraction=np.maximum(0,np.minimum(1,np.maximum(aa,bb))-np.maximum(0,np.minimum(aa,bb)))
            flat=np.abs(dz)<=1e-12;fraction[flat]=((z0[flat]>=lo)&(z0[flat]<=hi)).astype(float)
            out[label+'_band10_line_nm']=float(l@fraction)
        return out
