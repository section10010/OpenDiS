"""Translation-invariant native direct force evaluated in a nearby origin.

Only the coordinates supplied to the native force kernel are translated. The
physical graph, CTB, mobility, integration, stored coordinates and parameters
are unchanged. Applicable to nonperiodic infinite-medium direct forces only.
No force projection or balancing correction is performed.
"""
import copy
import numpy as np


class CenteredForce:
    def __init__(self,sim,raw,origin_nm):
        if sim.c['force_mode']!='direct' or any(sim.c['periodic']):raise NotImplementedError('Centered force is qualified only for nonperiodic direct infinite-medium forces')
        self.sim=sim;self.raw=raw;self.shift=np.array(origin_nm,float)/sim.bnm
        if self.shift.shape!=(3,) or not np.isfinite(self.shift).all():raise ValueError('Invalid force evaluation origin')

    def network(self,N):
        d=copy.deepcopy(N.get_disnet(self.sim.exa.ExaDisNet).export_data())
        d['nodes']['positions']-=self.shift;d['cell']['origin']-=self.shift
        g=self.sim.exa.ExaDisNet();g.import_data(d)
        return self.sim.exa.DisNetManager(g)

    def PreCompute(self,N,state):
        # Both admitted direct kernels precompute on each actual force call.
        return state

    def NodeForce(self,N,state,pre_compute=True):
        return self.raw.NodeForce(self.network(N),state,pre_compute=True)

    def OneNodeForce(self,N,state,tag,update_state=True,match_global=True):
        return self.raw.OneNodeForce(self.network(N),state,tag,update_state,match_global=True)


def install(sim,origin_nm):
    sim.force=CenteredForce(sim,sim.force,origin_nm)
    return sim
