"""Accepted native-topology lineage evidence; never generates dislocations."""
import numpy as np
from geometry import components
from configuration import save_json


def topology(a):
    ids=a['ids'];groups=components(len(a['positions_nm']),ids);degree=np.bincount(ids.ravel(),minlength=len(a['positions_nm']))
    result=[]
    for ix in groups:
        free=a['constraints'].ravel()[ix]!=7
        result.append(dict(ix=ix,tags={tuple(t) for t in a['tags'][ix]},anchored=bool(np.any(~free)),
                           free_cycle=bool(free.all() and np.all(degree[ix]==2))))
    return result


class Lineage:
    def __init__(self,sim,out):
        self.out=out;self.kind=sim.c['kind'];self.gauge=sim.c['source_gauge_nm'];self.initial=sim.arrays(False);self.events=[];self.children={};self.next_id=1
        self.tag_roots={};self.tag_lineages={}
        for source in sim.sources:
            for i in source['initial_node_indices']:
                tag=tuple(self.initial['tags'][i]);self.tag_roots[tag]={source['source_id']};self.tag_lineages[tag]={0}
        self.last=self.initial;self.groups=topology(self.last)
        for group in self.groups:group['lineage_id']=0
        (out/'events').mkdir(exist_ok=True)

    def update(self,sim,stage,step,time_ns):
        current=sim.arrays(False);groups=topology(current)
        previous_anchor={t for g in self.groups if g['anchored'] for t in g['tags']}
        prior_free={t for g in self.groups if not g['anchored'] for t in g['tags']}
        for group in groups:
            roots=set().union(*(self.tag_roots.get(t,set()) for t in group['tags']))
            # All W3 dynamical source tests are isolated single-source cases.
            # Multi-parent collisions retain a set; no invented unique ancestor.
            if not roots:roots={-1}
            for tag in group['tags']:self.tag_roots[tag]=set(roots)
            # Recycled native tags can overlap an old debris loop. Assign the
            # predecessor with the largest actual tag overlap, never an any-tag veto.
            predecessor=max(self.groups,key=lambda g:len(g['tags']&group['tags']),default=None)
            group['lineage_id']=predecessor.get('lineage_id',0) if predecessor is not None else 0
            if self.kind=='FR_source' and not group['anchored'] and predecessor is not None and predecessor['anchored'] and (group['tags']&previous_anchor):
                cid=self.next_id;self.next_id+=1
                radius=self.radius(current,group)
                child=dict(lineage_id=cid,parent_source_ids=sorted(roots),birth_step=step,birth_time_ns=time_ns,
                           birth_stage=stage,confirmed=False,last_tags=set(group['tags']),persistence_steps=0,
                           birth_area_radius_nm=radius,max_area_radius_nm=radius,last_radius_nm=radius,
                           birth_free_cycle=group['free_cycle'])
                self.children[cid]=child
                group['lineage_id']=cid
                event=dict(event='native_detachment_candidate',lineage_id=cid,parent_source_ids=sorted(roots),
                           step=step,time_ns=time_ns,stage=stage,node_tags=sorted(group['tags']))
                self.events.append(event)
                stem=f'events/detachment_{cid:03d}_{step:07d}_{stage}'
                np.savez_compressed(self.out/(stem+'_before.npz'),**self.last)
                np.savez_compressed(self.out/(stem+'_after.npz'),**current)
                event['before_snapshot']=stem+'_before.npz';event['after_snapshot']=stem+'_after.npz'
                for tag in group['tags']:self.tag_lineages[tag]={cid}
            else:
                lineages=set().union(*(self.tag_lineages.get(t,set()) for t in group['tags'])) or {0}
                for tag in group['tags']:self.tag_lineages[tag]=set(lineages)
        self.last=current;self.groups=groups

    def accepted(self,step,time_ns):
        anchored=any(g['anchored'] and np.any(self.last['constraints'].ravel()[g['ix']]!=7) for g in self.groups)
        for cid,ch in self.children.items():
            # Component lineage is inherited from its largest predecessor. A
            # three-node debris loop absorbed by a large loop never inherits
            # that large loop's growth, even if native node tags are recycled.
            matches=[g for g in self.groups if g['free_cycle'] and g['lineage_id']==cid]
            if len(matches)==1 and anchored:
                group=matches[0];radius=self.radius(self.last,group)
                ch['last_tags']=set(group['tags']);ch['persistence_steps']+=1
                ch['max_area_radius_nm']=max(ch['max_area_radius_nm'],radius);ch['last_radius_nm']=radius
                if not ch['confirmed'] and ch['persistence_steps']>=3 and radius>=self.gauge/2 and radius>=1.1*ch['birth_area_radius_nm']:
                    ch['confirmed']=True
                    self.events.append(dict(event='native_detachment_confirmed',lineage_id=cid,parent_source_ids=ch['parent_source_ids'],
                                            step=step,time_ns=time_ns,birth_step=ch['birth_step'],persistence_steps=ch['persistence_steps']))
            elif not ch['confirmed']:ch['persistence_steps']=0

    @property
    def confirmed(self):
        active={g['lineage_id'] for g in self.groups if g['free_cycle']}
        return sum(c['confirmed'] and cid in active for cid,c in self.children.items())

    def save(self):
        children=[{k:sorted(v) if isinstance(v,set) else v for k,v in c.items()} for c in self.children.values()]
        save_json(self.out/'source_lineage.json',dict(events=self.events,children=children,confirmed_detachments=self.confirmed,
            interpretation='Largest tag-overlap predecessor assigns a unique component lineage. A debris component absorbed by another loop cannot borrow its area growth. Confirmation requires free closed cycle, mobile parent, >=3 accepted steps, area-equivalent radius >= half gauge, and >=10 percent radius growth. Count reports currently distinct qualified free components; debris and expired candidates remain visible.',
            initialization_injected_daughters=0))

    @staticmethod
    def radius(a,group):
        mask=np.isin(a['ids'][:,0],group['ix']);ij=a['ids'][mask]
        x=a['positions_nm']-a['positions_nm'][group['ix']].mean(axis=0)
        return float(np.sqrt(np.linalg.norm(.5*np.cross(x[ij[:,0]],x[ij[:,1]]).sum(axis=0))/np.pi))
