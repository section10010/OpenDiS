"""Stage-wise identity for the admitted single-slip-per-source I0 benchmark.

Native tags are addresses, not material IDs. A segment inherits a unique source
from its Burgers, glide normal AND affine source plane. Unrecognised junctions
are rejected, never guessed from z or a recycled tag. General reactions are W5.
"""
import numpy as np


class IdentityError(RuntimeError):
    pass


def network_arrays(sim, N=None):
    g=sim.g if N is None else N.get_disnet(sim.exa.ExaDisNet)
    d=g.export_data();n=d['nodes'];s=d['segs']
    return dict(positions_nm=np.asarray(n['positions'])*sim.bnm,tags=np.asarray(n['tags']),
        constraints=np.asarray(n['constraints']),ids=np.asarray(s['nodeids'],int),
        burgers=np.asarray(s['burgers']),planes=np.asarray(s['planes']))


class MaterialIdentity:
    def __init__(self,sources,interface,plane_tol_nm=1e-4):
        self.sources=sources;self.interface=interface;self.plane_tol_nm=plane_tol_nm
        self.previous={};self.serial=0;self.events=[];self.stage_counts={}

    def inspect(self,a):
        x=a['positions_nm'];ij=a['ids'];m=len(ij);edge=np.full(m,-1,int)
        for k,source in enumerate(self.sources):
            slip=source['selected_slip'];b=np.array(slip['b']);n=np.array(slip['n']);center=np.array(slip['center_nm'])
            # Signs of a segment and its normal may flip during legal remeshing.
            good=(np.linalg.norm(np.cross(a['burgers'],b),axis=1)<1e-8)&(np.abs(np.linalg.norm(a['burgers'],axis=1)-np.linalg.norm(b))<1e-8)
            good&=(np.linalg.norm(np.cross(a['planes'],n),axis=1)<1e-8)&(np.abs(np.linalg.norm(a['planes'],axis=1)-1)<1e-8)
            good&=np.all(np.abs((x[ij]-center)@n)<self.plane_tol_nm,axis=1)
            if np.any(good&(edge>=0)):raise IdentityError('Ambiguous source fingerprint; require native provenance before accepting')
            edge[good]=k
        if np.any(edge<0):raise IdentityError('Unrecognised segment reaction/source plane; no silent material fallback')
        # Exact integer incidence reduction; no force, geometry or floating
        # criterion changes. Avoid repeated Python edge loops at fine meshes.
        minimum=np.full(len(x),len(self.sources),int);maximum=np.full(len(x),-1,int)
        for end in (0,1):
            np.minimum.at(minimum,ij[:,end],edge);np.maximum.at(maximum,ij[:,end],edge)
        if np.any(maximum<0):raise IdentityError('Orphan node has no material provenance')
        if np.any(minimum!=maximum):raise IdentityError('Cross-source/material junction is illegal in I0')
        node=maximum
        hosts=np.array([s['region'] for s in self.sources],dtype='U1')
        stable=np.array([s['initial_component_id'] for s in self.sources],int)
        return dict(node_source_id=stable[node],segment_source_id=stable[edge],node_host=hosts[node],segment_host=hosts[edge])

    def commit(self,a,stage,time_ns):
        meta=self.inspect(a);current={}
        for i,tag in enumerate(a['tags']):
            key=tuple(map(int,tag));source=int(meta['node_source_id'][i]);old=self.previous.get(key)
            # Observation IDs are unique every stage. Native tag lifetime inside
            # one C++ operation is not exposed and is NOT asserted as genealogy.
            self.serial+=1;current[key]=(source,self.serial)
            if old is None or old[0]!=source:
                self.events.append(dict(event='material_observation',stage=stage,time_ns=time_ns,tag=key,
                    source_id=source,host_layer=meta['node_host'][i],observation_id=self.serial,
                    prior_address_source=(None if old is None else old[0]),reason='new address or address reused across source'))
        for key in self.previous.keys()-current.keys():
            self.events.append(dict(event='address_retired',stage=stage,time_ns=time_ns,tag=key,source_id=self.previous[key][0]))
        self.previous=current;self.stage_counts[stage]=self.stage_counts.get(stage,0)+1
        return meta

    def fields(self,a):
        meta=self.inspect(a);host=meta['node_host'];z=a['positions_nm'][:,2];z0=self.interface['z_nm']
        meta.update(host_layer=host,origin_region=host,current_region=np.where(np.abs(z-z0)<=1e-6,'CTB0',np.where(z<z0,'M','T')),
            lineage_id=meta['node_source_id'].copy(),canonical_interface_id=np.full(len(host),self.interface['id']),
            image_index=np.zeros(len(host),int))
        return meta
