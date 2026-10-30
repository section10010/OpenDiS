"""Non-mock geometry/identity regressions; dynamics are separately native."""
import copy,json,sys,unittest
from pathlib import Path
import numpy as np
from material_identity import MaterialIdentity,IdentityError
from shared_force_probe import initial_data
from specification import validate_spec

if 'configuration' not in sys.modules:sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'code'))
from configuration import config
from geometry import build


class Tests(unittest.TestCase):
    def setUp(self):
        specpath=Path(__file__).resolve().parents[1]/'configs/shared_evolution_v1.json'
        self.spec=json.loads(specpath.read_text());self.c=config(**self.spec['common_native_updates'])
        seed,_,self.sources=initial_data(self.c,'MT',build,self.spec)
        self.a=dict(positions_nm=seed['nodes']['positions']*self.c['b_m']*1e9,tags=seed['nodes']['tags'],constraints=seed['nodes']['constraints'],
            ids=seed['segs']['nodeids'],burgers=seed['segs']['burgers'],planes=seed['segs']['planes'])
        self.identity=MaterialIdentity(self.sources,self.spec['interface'])

    def test_original_source_ids(self):
        m=self.identity.inspect(self.a);self.assertEqual(set(m['node_source_id']),{0,1})
        self.assertTrue(np.all(m['node_host'][m['node_source_id']==1]=='T'))

    def test_wrong_source_plane_rejected(self):
        self.a['positions_nm'][0]+=np.array(self.sources[0]['selected_slip']['n'])*.1
        with self.assertRaises(IdentityError):self.identity.inspect(self.a)

    def test_unrecognised_burgers_rejected(self):
        self.a['burgers'][0]*=2
        with self.assertRaises(IdentityError):self.identity.inspect(self.a)

    def test_ambiguous_provenance_rejected(self):
        duplicate=copy.deepcopy(self.sources[0]);duplicate['region']='T'
        with self.assertRaises(IdentityError):MaterialIdentity([*self.sources,duplicate],self.spec['interface']).inspect(self.a)

    def test_tag_reuse_does_not_control_host(self):
        self.identity.commit(self.a,'initial',0)
        seed,_,sources=initial_data(self.c,'T',build,self.spec)
        a=dict(positions_nm=seed['nodes']['positions']*self.c['b_m']*1e9,tags=seed['nodes']['tags'].copy(),constraints=seed['nodes']['constraints'],
            ids=seed['segs']['nodeids'],burgers=seed['segs']['burgers'],planes=seed['segs']['planes'])
        a['tags'][0]=self.a['tags'][0]
        result=self.identity.commit(a,'remesh_fixture',0)
        self.assertEqual(result['node_host'][0],'T');self.assertEqual(result['node_source_id'][0],1)
        self.assertTrue(any(e.get('prior_address_source')==0 and e.get('host_layer')=='T' for e in self.identity.events))

    def test_origin_is_not_inferred_from_current_z(self):
        meta=self.identity.inspect(self.a);source=self.sources[1];n=np.array(source['selected_slip']['n'])
        direction=np.array([0.,0.,1.])-n*n[2]
        self.a['positions_nm'][meta['node_host']=='T']-=400*direction
        result=self.identity.inspect(self.a)
        self.assertTrue(np.all(result['node_host'][meta['node_host']=='T']=='T'))
        self.assertTrue(np.all(self.a['positions_nm'][meta['node_host']=='T',2]<2000))

    def test_unimplemented_physics_refused(self):
        validate_spec(self.spec,self.c)
        bad=copy.deepcopy(self.spec);bad['transmission_rule']='invented'
        with self.assertRaises(ValueError):validate_spec(bad,self.c)
        with self.assertRaises(NotImplementedError):validate_spec(self.spec,config(**{**self.c,'kind':'FR_source'}))
        with self.assertRaises(NotImplementedError):validate_spec(self.spec,config(**{**self.c,'source_count':2}))
        bad=copy.deepcopy(self.spec);bad['regions'][1]['allowed_normal']=[0,0,-1]
        with self.assertRaises(ValueError):validate_spec(bad,self.c)


if __name__=='__main__':unittest.main()
