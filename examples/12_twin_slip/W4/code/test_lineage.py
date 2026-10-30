"""Graph-only regression fixtures, never input to the dynamical solver.

Tests the bookkeeping failures found in retained real native trajectories:
junction-bearing detachments and microscopic debris absorbed into old loops.
"""
import tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from events import Lineage


def network(cycles,radii,pins=(0,1)):
    tags=sorted(set(t for cycle in cycles for t in cycle));index={t:i for i,t in enumerate(tags)}
    x=np.zeros((len(tags),3));edges=[]
    for ci,(cycle,radius) in enumerate(zip(cycles,radii)):
        for j,t in enumerate(cycle):x[index[t]]=[radius*np.cos(2*np.pi*j/len(cycle)),radius*np.sin(2*np.pi*j/len(cycle))+ci*1000,0]
        edges.extend((index[cycle[j]],index[cycle[(j+1)%len(cycle)]]) for j in range(len(cycle)))
    return dict(positions_nm=x,tags=np.array([[0,t] for t in tags]),constraints=np.array([[7 if t in pins else 0] for t in tags]),
        ids=np.array(edges),burgers=np.tile([1.,0,0],(len(edges),1)),planes=np.tile([0.,0,1],(len(edges),1)),cell_nm=np.eye(3)*10000)


class LineageRegression(unittest.TestCase):
    def test_junction_birth_and_absorbed_debris_do_not_double_count(self):
        holder=[network([[0,1,2,3,4,5,7]],[200])]
        sim=SimpleNamespace(c=dict(kind='FR_source',source_gauge_nm=240.),arrays=lambda _:holder[0],sources=[dict(source_id=0,initial_node_indices=list(range(7)))])
        with tempfile.TemporaryDirectory() as tmp:
            track=Lineage(sim,Path(tmp))
            # Free component has a degree-four junction at tag 2.
            holder[0]=network([[0,1,6],[2,3,4],[2,5,7]],[20,150,150]);track.update(sim,'collision',1,0.)
            self.assertEqual(len(track.children),1);self.assertFalse(track.children[1]['birth_free_cycle'])
            for step in (2,3,4):
                holder[0]=network([[0,1,6],[2,3,4,5,7]],[20,250+step]);track.update(sim,'topology',step,step*.01);track.accepted(step,step*.01)
            self.assertEqual(track.confirmed,1)
            holder[0]=network([[0,1,30],[2,3,4,5,7],[6,20,21]],[20,260,1]);track.update(sim,'topology',5,.05);track.accepted(5,.05)
            self.assertEqual(len(track.children),2)
            for step in (6,7,8,9):
                holder[0]=network([[0,1,30],[2,3,4,5,7,6,20,21]],[20,300+step]);track.update(sim,'remesh',step,step*.01);track.accepted(step,step*.01)
            self.assertEqual(track.confirmed,1);self.assertFalse(track.children[2]['confirmed'])


if __name__=='__main__':unittest.main()
