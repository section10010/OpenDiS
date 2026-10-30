"""Archive native evidence once; a small view excludes only dense velocity chunks.

No hashes, deletion, rewriting of original results, or fallback fabricated data.
"""
import argparse,json,tarfile,zipfile
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--tar',type=Path,required=True);p.add_argument('--view-zip',type=Path,required=True);a=p.parse_args();r=a.root.resolve()
    if a.tar.exists() or a.view_zip.exists():raise RuntimeError('Refuse archive overwrite')
    files=[x for x in sorted(r.rglob('*')) if x.is_file() and '__pycache__' not in x.parts and 'reference_tmp' not in x.relative_to(r).parts]
    manifest=dict(schema='W4.G0.FileInventory.v1',files=[dict(path=str(x.relative_to(r)),bytes=x.stat().st_size) for x in files],
        total_bytes=sum(x.stat().st_size for x in files),integrity='Size inventory; NPZ archives independently read during native velocity audit. No SHA256 inventory.',
        view_note='View ZIP omits accepted_velocity chunks only. Full TAR contains them; do not claim a local full velocity audit without extracting the full TAR.')
    mp=r/'delivery_inventory.json';mp.write_text(json.dumps(manifest,indent=2),encoding='utf-8');files.append(mp)
    with tarfile.open(a.tar,'w') as archive:
        for path in files:archive.add(path,arcname=path.relative_to(r),recursive=False)
    with zipfile.ZipFile(a.view_zip,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for path in files:
            if 'accepted_velocity' not in path.relative_to(r).parts:archive.write(path,path.relative_to(r))
        archive.writestr('VIEW_README.txt',manifest['view_note'])
    print('NATIVE_ARCHIVE_SAVED',len(files),'TAR_BYTES',a.tar.stat().st_size,'VIEW_BYTES',a.view_zip.stat().st_size)


if __name__=='__main__':main()
