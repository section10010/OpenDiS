"""Native origin repair qualification on the exact formerly failed geometry."""
import argparse,json,copy,os
from pathlib import Path
import numpy as np
from run_shared import bootstrap,seed_from_npz
from material_identity import MaterialIdentity
from centered_force import install,CenteredForce


def main():
    p=argparse.ArgumentParser();p.add_argument('--historical-root',type=Path,required=True);p.add_argument('--g0-root',type=Path,required=True);p.add_argument('--repo',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();bootstrap(a.g0_root)
    from native import Simulation,imports
    from configuration import config,save_json
    if a.out.exists():raise RuntimeError('Refuse overwrite')
    a.out.mkdir(parents=True);os.environ.update(OMP_NUM_THREADS='2',OMP_PROC_BIND='false');py,_=imports(a.repo);py.initialize(num_threads=2)
    try:
        parent=a.historical_root/'cases/mesh2p5_MT_budget_r2';c=config(**json.loads((parent/'spec.json').read_text())['common_native_updates'])
        seed=seed_from_npz(parent/'final.npz',c);sources=json.loads((parent/'sources.json').read_text());origin=[2000.,2000.,2000.]
        sim=install(Simulation(c,a.repo,initial=seed),origin);before=sim.snapshot();f,_=sim.evaluate();f=f.copy();after=sim.snapshot()
        for k in ('tags','positions','constraints'):np.testing.assert_array_equal(before['nodes'][k],after['nodes'][k])
        for k in ('nodeids','burgers','planes'):np.testing.assert_array_equal(before['segs'][k],after['segs'][k])
        translated=copy.deepcopy(seed);translated['nodes']['positions']-=np.array(origin)/sim.bnm;translated['cell']['origin']-=np.array(origin)/sim.bnm
        manual=Simulation(c,a.repo,initial=translated);fm,_=manual.evaluate();np.testing.assert_allclose(f,fm,rtol=1e-12,atol=1e-8)
        meta=MaterialIdentity(sources,dict(id='CTB0',z_nm=2000.)).inspect(sim.arrays(False));single=np.zeros_like(f)
        for host in ('M','T'):
            node=meta['node_host']==host;edge=meta['segment_host']==host;ix=np.flatnonzero(node);remap=np.full(len(node),-1);remap[ix]=np.arange(len(ix))
            sub=dict(cell=seed['cell'],nodes={k:v[node] for k,v in seed['nodes'].items()},segs=dict(nodeids=remap[seed['segs']['nodeids'][edge]],burgers=seed['segs']['burgers'][edge],planes=seed['segs']['planes'][edge]))
            one=install(Simulation(c,a.repo,initial=sub),origin);sf,_=one.evaluate();single[node]=sf
        cross=(f-single)*sim.b**2;balance=float(np.linalg.norm(cross.sum(axis=0))/np.linalg.norm(cross,axis=1).sum());assert balance<1e-8
        np.savez_compressed(a.out/'formerly_failed_native_state.npz',**sim.arrays(),centered_cross_force_N=cross)
        rejected=False
        try:CenteredForce(type('Dummy',(),dict(c={'force_mode':'fft','periodic':[True]*3}))(),None,origin)
        except NotImplementedError:rejected=True
        assert rejected
        checks=[dict(test='physical_graph_unchanged',status='PASS'),dict(test='independent_native_translated_graph_force',status='PASS'),
            dict(test='formerly_failed_reciprocity',status='PASS',residual=balance,unchanged_threshold=1e-8),dict(test='periodic_FFT_refused',status='PASS')]
        save_json(a.out/'result.json',dict(status='PASS',checks=checks,repair='Uniform translation of force kernel input only; no force balancing projection or changes to source/contact/load/mobility',historical_geometry=str(parent/'final.npz')))
        print('CENTERED_FORCE_NATIVE_PASS',len(checks),balance)
    finally:py.finalize()


if __name__=='__main__':main()
