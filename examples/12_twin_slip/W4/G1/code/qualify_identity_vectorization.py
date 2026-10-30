"""Exact source-label equivalence against immutable pre-optimization snapshots."""
import argparse,importlib.util,json,time
from pathlib import Path
import numpy as np
from material_identity import MaterialIdentity


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--historical-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise RuntimeError('Refuse overwrite')
    a.out.mkdir(parents=True);path=a.root/'cases/base_MT/code_snapshot/material_identity.py'
    module_spec=importlib.util.spec_from_file_location('identity_reference',path);reference=importlib.util.module_from_spec(module_spec);module_spec.loader.exec_module(reference)
    checked=0;nodes=0;start=time.monotonic()
    for root in (a.historical_root,a.root):
        for folder in sorted((root/'cases').iterdir()):
            if not (folder/'sources.json').exists() or not (folder/'result.json').exists():continue
            sources=json.loads((folder/'sources.json').read_text());interface=json.loads((folder/'spec.json').read_text())['interface']
            old=reference.MaterialIdentity(sources,interface);new=MaterialIdentity(sources,interface)
            for file in [folder/'initial.npz',*sorted(folder.glob('frame_*.npz')),folder/'final.npz']:
                if not file.exists():continue
                with np.load(file) as z:
                    data={k:z[k] for k in ('positions_nm','tags','constraints','ids','burgers','planes')}
                x=old.inspect(data);y=new.inspect(data)
                for key in x:np.testing.assert_array_equal(x[key],y[key])
                checked+=1;nodes+=len(data['tags'])
    assert checked>100
    (a.out/'result.json').write_text(json.dumps(dict(status='PASS',real_saved_graphs=checked,node_observations=nodes,
        comparison='All source/host arrays exactly equal; unchanged floating predicates and source signatures',reference=str(path),wall_s=time.monotonic()-start),indent=2));print('IDENTITY_INTEGER_EQUIVALENCE_PASS',checked,nodes)


if __name__=='__main__':main()
