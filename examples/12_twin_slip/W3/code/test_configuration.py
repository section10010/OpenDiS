import unittest
from copy import deepcopy
import numpy as np
from configuration import config
from geometry import build,systems,pure_shear,residual


class ConfigurationTests(unittest.TestCase):
    def test_geometry_parameters_change_actual_arrays(self):
        c=config();rn,lk,ss,roles,_=build(c)
        for key,value in [('source_gauge_nm',360.),('backbone_offset_nm',160.),('source_count',2)]:
            d=config(**{key:value});rr,ll,sd,_,_=build(d)
            self.assertFalse(rn.shape==rr.shape and np.array_equal(rn,rr),key)
            self.assertEqual(len(sd),d['source_count'])
            s=sd[0];pts=rr[s['initial_node_indices'],:3];center=np.array(s['center_nm'])
            self.assertAlmostEqual(float(np.ptp((pts-center)@s['q'])),d['source_gauge_nm'])
            self.assertAlmostEqual(float(np.ptp((pts-center)@s['n'])),d['backbone_offset_nm'])
            self.assertLess(residual(len(rr),ll[:,:2].astype(int),ll[:,2:5]).max(),1e-12)
    def test_global_quotas_and_500nm(self):
        for width in (200,500,1000,3000):
            c=config(kind='multilamellar',box_nm=[12000.,12000.,18000.],hT_nm=float(width),hM_nm=2.*width)
            rn,lk,ss,_,_=build(c);self.assertEqual(len(ss),180)
            for phase,total in [('twin',60),('matrix',120)]:
                subset=[s for s in ss if s['phase']==phase];self.assertEqual(len(subset),total)
                self.assertEqual([sum(s['slip_index']==i for s in subset) for i in range(12)],[total//12]*12)
            self.assertLess(residual(len(rn),lk[:,:2].astype(int),lk[:,2:5]).max(),1e-12)
        c=config(kind='multilamellar',box_nm=[12000.,12000.,18000.],hT_nm=500.,hM_nm=1000.,sources_twin=120)
        self.assertEqual(len(build(c)[2]),240)
    def test_stress_controls_burgers_and_rss(self):
        for phase in ('matrix','twin'):
            for mode in ('HI','HII','soft'):
                sigma=pure_shear(phase,mode,100e6)
                c=config(phase=phase,mode=mode,stress_pa=sigma.tolist());_,_,s,_,_=build(c)
                d=config(phase=phase,mode=mode,stress_pa=(-sigma).tolist());_,_,r,_,_=build(d)
                np.testing.assert_allclose(s[0]['b'],-np.array(r[0]['b']),atol=1e-12)
                self.assertAlmostEqual(s[0]['rss_pa'],np.array(s[0]['b'])@sigma@np.array(s[0]['n']))
                self.assertLess(s[0]['rss_pa'],0.)
    def test_reject_unused_or_unsupported_settings(self):
        with self.assertRaises(ValueError):config(gauge=120)
        with self.assertRaises(ValueError):config(friction_pa=1e6)
        with self.assertRaises(ValueError):config(Medge=0)
    def test_crystallography(self):
        for phase in ('matrix','twin'):
            s=systems(phase,np.eye(3));self.assertEqual(len(s),12)
            self.assertEqual([sum(v['mode']==m for v in s) for m in ('HI','HII','soft')],[6,3,3])
            for v in s:self.assertLess(abs(v['b']@v['n']),1e-12)


if __name__=='__main__':unittest.main()
