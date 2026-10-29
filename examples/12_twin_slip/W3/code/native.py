"""Real ExaDiS force/mobility/integration/collision/topology/remesh adapters."""
import os,sys
from pathlib import Path
import numpy as np
from geometry import build,residual,displacements


def imports(repo):
    for rel in ('python','lib','core/pydis/python','core/exadis/python'):
        sys.path.insert(0,str(Path(repo)/rel))
    import pyexadis
    import pyexadis_base as exa
    return pyexadis,exa


class Simulation:
    def __init__(self,c,repo,initial=None):
        self.c=c;self.b=c['b_m'];self.bnm=self.b*1e9
        self.py,self.exa=imports(repo)
        self.state=dict(crystal='fcc',burgmag=self.b,mu=c['mu_pa'],nu=c['nu'],a=c['core_a_b'],
            maxseg=c['maxseg_nm']/self.bnm,minseg=c['minseg_nm']/self.bnm,rann=c['rann_nm']/self.bnm,
            rtol=c['rtol_nm']/self.bnm,maxdt=c['maxdt_ps']*1e-12,nextdt=c['firstdt_ps']*1e-12,
            use_glide_planes=0,enforce_glide_planes=0,
            applied_stress=np.array(c['stress_pa']).ravel()[[0,4,8,5,2,1]])
        if initial is None:
            rn,lk,self.sources,self.initial_roles,self.initial_source_ids=build(c)
            rn[:,:3]/=self.bnm
            self.g=self.exa.ExaDisNet(self.py.Cell(h=np.diag(np.array(c['box_nm'])/self.bnm),is_periodic=c['periodic']),rn,lk)
        else:
            self.g=self.exa.ExaDisNet();self.g.import_data(initial)
            self.sources=[];self.initial_roles=np.array([]);self.initial_source_ids=np.array([],int)
        self.net=self.exa.DisNetManager(self.g)
        if c['force_mode']=='direct':
            self.force=self.exa.CalForce(state=self.state,force_mode='GLOBAL_MODEL',cell=self.g.cell,
                force_list={'FORCE_CORE_SELF_PKEXT':{},'FORCE_N2_ISO':{}})
        else:
            self.force=self.exa.CalForce(state=self.state,force_mode='DDD_FFT_MODEL',cell=self.g.cell,Ngrid=c['grid'])
            cutoff=float(self.force.force.neighbor_cutoff)*self.bnm
            if cutoff+c['maxseg_nm']>min(c['box_nm'])/3:raise ValueError('FFT neighbor completeness bound failed')
        self.glide=self.exa.MobilityLaw(state=self.state,mobility_law='GLIDE',Medge=c['Medge'],Mscrew=c['Mscrew'])
        self.contact=None
        owner=self
        class Force:
            def PreCompute(self,N,state):return owner.force.PreCompute(N,state)
            def NodeForce(self,N,state,pre_compute=True):return owner.force.NodeForce(N,state,pre_compute)
            def OneNodeForce(self,N,state,tag,update_state=True):return owner.force.OneNodeForce(N,state,tag,update_state,match_global=True)
        class Mobility:
            non_linear=True
            def Mobility(self,N,state):
                owner.glide.Mobility(N,state)
                if owner.contact is not None:state['nodevels']=owner.contact.project_all(owner,N,state,state['nodevels'])
                state['nodevels']=owner.cap(state['nodevels']);return state
            def OneNodeMobility(self,N,state,tag,f,update_state=True):
                v=owner.glide.OneNodeMobility(N,state,tag,f,False)
                if owner.contact is not None:v=owner.contact.project_one(owner,N,state,tag,v)
                return owner.cap(v)
        self.fadapter=Force();self.mobility=Mobility()
        self.integrator=self.exa.TimeIntegration(state=self.state,integrator=c['integrator'],dt=c['maxdt_ps']*1e-12,
            force=self.fadapter,mobility=self.mobility)
        self.collision=self.exa.Collision(state=self.state,collision_mode='Retroactive')
        self.topology=self.exa.Topology(state=self.state,topology_mode='TopologySerial',force=self.fadapter,mobility=self.mobility)
        self.remesh=self.exa.Remesh(state=self.state,remesh_rule='LengthBased',coarsen_mode=1)

    def cap(self,v):
        limit=self.c['vmax_m_s']
        if limit is None:return v
        return v*np.minimum(1.,limit/self.b/np.maximum(np.linalg.norm(v,axis=-1),1e-99))[...,None]

    def evaluate(self):
        self.force.NodeForce(self.net,self.state);self.mobility.Mobility(self.net,self.state)
        if not np.isfinite(self.state['nodeforces']).all() or not np.isfinite(self.state['nodevels']).all():raise FloatingPointError('Nonfinite native force/velocity')
        return self.state['nodeforces'],self.state['nodevels']

    def advance(self,dt=None):
        self.state['oldnodes_dict']=self.g.get_nodes_data()
        if dt is not None:self.integrator.dt=dt
        self.integrator.Update(self.net,self.state)

    def operations(self):
        if self.c['collision']:self.collision.HandleCol(self.net,self.state)
        if self.contact is not None:self.contact.validate(self,'collision')
        self.force.PreCompute(self.net,self.state)
        if self.c['topology']:self.topology.Handle(self.net,self.state)
        if self.contact is not None:self.contact.validate(self,'topology')
        if self.c['remesh']:self.remesh.Remesh(self.net,self.state)
        if self.contact is not None:self.contact.validate(self,'remesh')

    def snapshot(self):return self.g.export_data()

    def arrays(self,with_fields=True):
        d=self.g.export_data();n=d['nodes'];s=d['segs']
        a=dict(positions_nm=n['positions']*self.bnm,tags=n['tags'],constraints=n['constraints'],
            ids=s['nodeids'],burgers=s['burgers'],planes=s['planes'],cell_nm=np.array(d['cell']['h'])*self.bnm)
        if with_fields:
            f,v=self.evaluate();a.update(force_N=f*self.b**2,velocity_m_s=v*self.b)
        return a

    def measures(self):
        a=self.arrays(False);dr=displacements(a['positions_nm'],a['ids'],self.c['box_nm'],self.c['periodic']);length=np.linalg.norm(dr,axis=1)
        free=a['constraints'].ravel()!=7;mobile=free[a['ids']].any(axis=1)
        return dict(nodes=len(free),segments=len(length),line_nm=float(length.sum()),mobile_line_nm=float(length[mobile].sum()),
            max_burgers_residual_b=float(residual(len(free),a['ids'],a['burgers']).max(initial=0)))
