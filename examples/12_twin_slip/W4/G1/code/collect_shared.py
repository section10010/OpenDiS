"""Complete raw archive and small review view; size inventories, no hashing."""
import argparse,json,tarfile,zipfile
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--dest',type=Path,required=True);a=p.parse_args()
    label=a.root.name.rsplit('_',1)[-1]
    if label not in ('002','003'):raise ValueError('Explicit G1 run root 002/003 required')
    full=a.dest/f'G1_native_evidence_20261030_{label}.tar';view=a.dest/f'G1_native_view_20261030_{label}.zip'
    if full.exists() or view.exists():raise RuntimeError('Refuse overwrite')
    files=[p for p in a.root.rglob('*') if p.is_file() and '__pycache__' not in p.parts and 'report_cache' not in p.parts]
    inventory={str(p.relative_to(a.root)):p.stat().st_size for p in files}
    (a.root/'delivery_inventory.json').write_text(json.dumps(inventory,indent=2));files.append(a.root/'delivery_inventory.json')
    with tarfile.open(full,'w') as tar:
        for p in files:tar.add(p,arcname=str(p.relative_to(a.root)),recursive=False)
    small=[p for p in files if 'accepted_velocity' not in p.parts]
    with zipfile.ZipFile(view,'w',zipfile.ZIP_DEFLATED) as z:
        for p in small:z.write(p,p.relative_to(a.root))
    print('G1_COMPLETE_ARCHIVE',len(files),full.stat().st_size,'VIEW',len(small),view.stat().st_size)


if __name__=='__main__':main()
