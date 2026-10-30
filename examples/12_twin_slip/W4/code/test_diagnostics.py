"""Bookkeeping unit tests, not synthetic DDD trajectories."""
import tempfile,unittest
from pathlib import Path
import numpy as np
from configuration import config
from diagnostics import Diagnostics,node_weights


class DiagnosticTests(unittest.TestCase):
    def test_exact_line_time_weighting_and_ragged_roundtrip(self):
        a=dict(positions_nm=np.array([[0.,0,0],[2.,0,0],[2.,6,0]]),ids=np.array([[0,1],[1,2]]),tags=np.array([[0,0],[0,1],[0,2]]))
        np.testing.assert_allclose(node_weights(a),[1,4,3])
        with tempfile.TemporaryDirectory() as t:
            d=Diagnostics(Path(t),config(vmax_m_s=1.,walls=[],diagnostic_chunk_steps=2))
            v=np.array([[2.,0,0],[0.,0,0],[0.,0,0]])
            fields=dict(raw_velocity_m_s=v,constrained_velocity_m_s=v,capped_velocity_m_s=np.minimum(v,1.),reaction_N=np.zeros((3,3)),multipliers_N=np.empty((3,0)))
            d.append(1,0.,2.,a,fields,np.ones((3,3)))
            fields={**fields,'raw_velocity_m_s':v*0,'constrained_velocity_m_s':v*0,'capped_velocity_m_s':v*0}
            d.append(2,2.,6.,a,fields,np.ones((3,3)));r=d.save()
            self.assertAlmostEqual(r['line_time_cap_fraction'],2/(8*8));self.assertEqual(r['any_node_capped_record_fraction'],.5)
            with np.load(Path(t)/'accepted_velocity/chunk_00000.npz') as z:
                np.testing.assert_array_equal(z['offsets'],[0,3,6]);np.testing.assert_array_equal(z['dt_s'],[2.,6.])

    def test_unknown_budget_is_rejected(self):
        with self.assertRaises(ValueError):config(max_nodes=0)
        with self.assertRaises(ValueError):config(magic_budget=1)


if __name__=='__main__':unittest.main()
