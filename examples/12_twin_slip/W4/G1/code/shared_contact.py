"""One physical CTB, opposite unilateral constraints on its two material sides."""
import numpy as np
from material_identity import network_arrays


class SharedContact:
    def __init__(self,identity,tol=1e-6,enabled=True):
        self.identity=identity;self.tol=tol;self.enabled=enabled
        self.z=identity.interface['z_nm'];self.walls=[identity.interface]
        self.reaction_native=np.empty((0,3));self.multipliers_native=np.empty((0,1))
        self.events=[];self.active=set();self.minimum_gap_nm=float('inf')

    def geometry(self,sim,N=None):
        a=network_arrays(sim,N);meta=self.identity.inspect(a)
        sign=np.where(meta['node_host']=='M',-1.,1.)
        return a,meta,sign,sign*(a['positions_nm'][:,2]-self.z)

    def project_one(self,sim,N,state,tag,v,context=None):
        a,meta,sign,gap=self.geometry(sim,N) if context is None else context
        matches=np.flatnonzero(np.all(a['tags']==tag,axis=1))
        if len(matches)!=1:raise RuntimeError('Ambiguous/missing native mobility address')
        i=int(matches[0]);v=np.array(v,float)
        if not self.enabled or gap[i]>self.tol or sign[i]*v[2]>=0:return v
        n=np.array([0.,0.,sign[i]])
        g=N.get_disnet(sim.exa.ExaDisNet)
        mn=np.array(sim.py.compute_node_mobility(g.net,i,mobility=sim.glide.mobility,fi=n))
        denom=n@mn
        if denom<=0:raise RuntimeError('Singular shared CTB mobility')
        multiplier=-(n@v)/denom;out=v+mn*multiplier
        if context is not None:
            self.reaction_native[i]=multiplier*n;self.multipliers_native[i,0]=multiplier
        return out

    def project_all(self,sim,N,state,velocities):
        context=self.geometry(sim,N);a,_,sign,gap=context;v=np.array(velocities,copy=True)
        if not np.array_equal(state['nodeveltags'],a['tags']):raise RuntimeError('Native mobility tag ordering mismatch')
        self.reaction_native=np.zeros_like(v);self.multipliers_native=np.zeros((len(v),1))
        if self.enabled:
            for i in np.flatnonzero((gap<=self.tol)&(sign*v[:,2]<0)):
                v[i]=self.project_one(sim,N,state,a['tags'][i],v[i],context)
        return v

    def event_dt(self,sim,v):
        if not self.enabled:return float('inf')
        _,_,sign,gap=self.geometry(sim);speed=-sign*v[:,2]*sim.bnm
        mask=(gap>self.tol)&(speed>0)
        return float(np.min(gap[mask]/speed[mask])) if mask.any() else float('inf')

    def validate(self,sim,stage):
        a,meta,sign,gap=self.geometry(sim)
        if self.enabled:
            self.minimum_gap_nm=min(self.minimum_gap_nm,float(gap.min(initial=float('inf'))))
            if np.any(gap < -self.tol):raise RuntimeError('Illegal shared CTB penetration after '+stage)
        self.identity.commit(a,stage,sim.state.get('time',0.)*1e9)
        active={(meta['node_host'][i],*map(int,a['tags'][i])) for i in np.flatnonzero(gap<=self.tol)} if self.enabled else set()
        if active!=self.active:
            self.events.append(dict(event='contact_membership',stage=stage,time_ns=sim.state.get('time',0.)*1e9,
                canonical_interface_id=self.identity.interface['id'],entered=sorted(active-self.active),exited=sorted(self.active-active),
                note='Address membership, including remesh; not an independent emission/impact count'))
        self.active=active

    def metrics(self,sim):
        a,meta,sign,gap=self.geometry(sim);out={};ij=a['ids']
        lengths=np.linalg.norm(a['positions_nm'][ij[:,1]]-a['positions_nm'][ij[:,0]],axis=1)
        active=(gap<=self.tol)&self.enabled
        for host in ('M','T'):
            node=meta['node_host']==host;edge=meta['segment_host']==host
            out.update({host+'_line_nm':float(lengths[edge].sum()),host+'_gap_nm':float(gap[node].min()) if node.any() else None,
                host+'_contact_nodes':int((active&node).sum()),host+'_tail_nm':float(lengths[edge&np.all(active[ij],axis=1)].sum())})
        return out
