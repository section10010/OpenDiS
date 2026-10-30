"""Regression of actual SerialDisNet pin-preserving merges (not a mock)."""
import os,sys,unittest
from pathlib import Path
import numpy as np

repo=Path(os.environ.get('OPENDIS_ROOT','/home/sakura/opendis/src'))
for rel in ('python','lib','core/pydis/python','core/exadis/python'):sys.path.insert(0,str(repo/rel))
import pyexadis
from pyexadis_base import ExaDisNet


class PinMerge(unittest.TestCase):
    def fixture(self):
        p=np.array([[100,100,100],[120,100,100],[140,100,100],
                    [160,100,100],[140,120,100],[120,120,100]],float)
        rn=np.column_stack([p,[7,0,0,7,0,0]])
        lk=np.array([[i,(i+1)%6,1,0,0,0,0,1] for i in range(6)],float)
        g=ExaDisNet(pyexadis.Cell(h=1000*np.eye(3),is_periodic=[True]*3),rn,lk)
        return g

    def test_reject_move_fixed_pin(self):
        g=self.fixture();s=g.net._get_serial_network()
        before=g.get_positions().copy()
        self.assertTrue(s.merge_nodes_position(0,1,np.array([105.,100,100]),np.zeros((3,3))))
        np.testing.assert_array_equal(g.get_positions(),before)

    def test_pin_tag_survives_when_second_argument(self):
        g=self.fixture();s=g.net._get_serial_network()
        self.assertFalse(s.merge_nodes_position(1,0,np.array([100.,100,100]),np.zeros((3,3))))
        s.purge_network()
        tags=g.get_tags();ix=np.flatnonzero((tags==[0,0]).all(axis=1))
        self.assertEqual(len(ix),1)
        self.assertEqual(g.get_nodes_data()['constraints'][ix[0],0],7)
        np.testing.assert_array_equal(g.get_positions()[ix[0]],[100,100,100])

    def test_distinct_pins_cannot_merge(self):
        g=self.fixture();s=g.net._get_serial_network()
        self.assertTrue(s.merge_nodes_position(0,3,np.array([130.,100,100]),np.zeros((3,3))))


if __name__=='__main__':
    pyexadis.initialize(num_threads=2)
    try:unittest.main()
    finally:pyexadis.finalize()
