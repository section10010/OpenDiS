"""Execute an exact argv, retaining command, output, exit code and timeout."""
import argparse,subprocess,time,os,signal,sys
from pathlib import Path
from configuration import save_json


def run_logged(cmd,cwd,receipt,log,timeout):
    receipt=Path(receipt);log=Path(log);receipt.parent.mkdir(parents=True,exist_ok=True);log.parent.mkdir(parents=True,exist_ok=True)
    if receipt.exists() or log.exists():raise RuntimeError('Refusing to overwrite command evidence')
    start=time.monotonic();timed_out=False
    with log.open('w',encoding='utf-8') as stream:
        p=subprocess.Popen(cmd,cwd=cwd,stdout=stream,stderr=subprocess.STDOUT,start_new_session=os.name!='nt')
        try:code=p.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out=True
            if os.name!='nt':os.killpg(p.pid,signal.SIGTERM)
            else:p.terminate()
            try:code=p.wait(timeout=5)
            except subprocess.TimeoutExpired:
                if os.name!='nt':os.killpg(p.pid,signal.SIGKILL)
                else:p.kill()
                code=p.wait()
    data=dict(command=cmd,cwd=str(cwd),exit_code=code,timed_out=timed_out,wall_seconds=time.monotonic()-start,
              log=str(log),omp_num_threads=os.environ.get('OMP_NUM_THREADS'))
    save_json(receipt,data);print('COMMAND_EXIT',code,'TIMEOUT',timed_out,'WALL_S',round(data['wall_seconds'],2),flush=True)
    return data


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cwd',required=True);p.add_argument('--receipt',required=True);p.add_argument('--log',required=True);p.add_argument('--timeout',type=float,default=600)
    p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args();cmd=a.command[1:] if a.command[:1]==['--'] else a.command
    r=run_logged(cmd,a.cwd,a.receipt,a.log,a.timeout);sys.exit(r['exit_code'] or (124 if r['timed_out'] else 0))
