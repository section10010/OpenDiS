"""Requalify a fresh runtime without borrowing another installation's A0 PASS."""
import argparse,os,shutil,sys
from pathlib import Path
from record import run_logged


def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--build',type=Path,required=True);p.add_argument('--root',type=Path,required=True);a=p.parse_args()
    root=a.root.resolve()
    if root.exists() and any(root.iterdir()):raise RuntimeError('A0 requires a new empty root')
    root.mkdir(parents=True,exist_ok=True);code=root/'code';shutil.copytree(Path(__file__).resolve().parent,code,ignore=shutil.ignore_patterns('__pycache__'))
    out=root/'A0';out.mkdir();os.environ.update(OPENDIS_ROOT=str(a.repo.resolve()),OMP_NUM_THREADS='2',OMP_PROC_BIND='false')
    os.environ['PATH']=str(Path(sys.executable).parent)+os.pathsep+os.environ['PATH']
    jobs=[('contracts',[sys.executable,'-m','unittest','test_configuration','-v'],code,'contracts_final.log'),
      ('native_v2',[sys.executable,'gate_a0_native.py','--repo',str(a.repo.resolve()),'--out',str(out/'native_v2')],code,'native_v2.log'),
      ('upstream_units',['ctest','--test-dir',str(a.build.resolve()),'-L','unit','--output-on-failure'],a.build,'upstream_units.log'),
      ('upstream_FR',['ctest','--test-dir',str(a.build.resolve()),'-L','full_run','-R','02_frank_read_src.*pydis_exadis_elast.*coarsen_mode','--output-on-failure'],a.build,'upstream_FR.log')]
    for name,cmd,cwd,log in jobs:
        r=run_logged(cmd,str(cwd),out/(name+'_receipt.json'),out/log,600)
        if r['exit_code']!=0 or r['timed_out']:raise RuntimeError('A0 blocked: '+name)
    import subprocess
    subprocess.run([sys.executable,str(code/'accept_a0.py'),'--root',str(root)],check=True)


if __name__=='__main__':main()
