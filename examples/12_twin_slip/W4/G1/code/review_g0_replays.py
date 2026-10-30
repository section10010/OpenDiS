"""Finish interrupted G0 review against its original frozen scientific contract.

The auxiliary 1e-10 full force/velocity replay assertion is retained as FAIL.
It is NOT relaxed or silently renamed PASS. Original G0 observables, graph and
invariants receive a separate verdict; their thresholds are read from G0.
"""
import argparse,json,sys,subprocess,os,time,csv
from pathlib import Path
import numpy as np


def main():
    p=argparse.ArgumentParser();p.add_argument('--g0-root',type=Path,required=True);p.add_argument('--run-root',type=Path,required=True);p.add_argument('--repo',required=True);a=p.parse_args()
    old=a.g0_root;new=a.run_root;criteria=json.loads((old/'frozen_criteria.json').read_text());os.environ.update(OMP_NUM_THREADS='2',OMP_PROC_BIND='false',OPENBLAS_NUM_THREADS='1')
    assert json.loads((new/'A0/gate_result.json').read_text())['status']=='PASS'
    for n in ('contact_admission','diagnostic_admission'):assert json.loads((new/n/'result.json').read_text())['status']=='PASS'
    code=list((old/'code').glob('*.py'));assert all(p.read_bytes()==(new/'code'/p.name).read_bytes() for p in code)
    checks=[];strict=[];primary_tol=criteria['primary_refinement_relative'];event_tol=criteria['event_time_refinement_relative']
    for family in ('loop_HI','threader_HII'):
        for wall in ('none','both'):
            name=f'G0_{family}_{wall}_uncapped_diagnostic';dest=new/'cases'/name
            if not (dest/'result.json').exists():
                cmd=[sys.executable,str(new/'code/run.py'),'--repo',a.repo,'--config',str(old/'configs'/(name+'.json')),'--out',str(dest),'--a0',str(new/'A0/gate_result.json')];t=time.monotonic()
                with (new/(name+'_r2.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=300)
                (new/(name+'_r2_receipt.json')).write_text(json.dumps(dict(command=cmd,exit_code=r.returncode,wall_s=time.monotonic()-t),indent=2));assert r.returncode==0
            x=json.loads((old/'cases'/name/'result.json').read_text());y=json.loads((dest/'result.json').read_text());assert y['status']=='PASS'
            errors={};absolute={}
            with np.load(old/'cases'/name/'final.npz') as u,np.load(dest/'final.npz') as v:
                for key in ('tags','ids','constraints','burgers','planes'):np.testing.assert_array_equal(u[key],v[key])
                for key in ('positions_nm','force_N','velocity_m_s'):
                    errors[key]=float(np.linalg.norm(u[key]-v[key])/max(np.linalg.norm(u[key]),1e-25));absolute[key]=float(np.abs(u[key]-v[key]).max(initial=0))
            strict.append(dict(case_id=name,criterion='Original auxiliary full-field relative error < 1e-10; unchanged',
                relative=errors,absolute_max=absolute,status='PASS' if max(errors.values())<1e-10 else 'FAIL'))
            area=abs(x['swept_area_nm2']-y['swept_area_nm2'])/abs(x['swept_area_nm2'])
            line=abs(x['final']['line_nm']-y['final']['line_nm'])/x['final']['line_nm']
            # Read first geometric contact from the original and replay histories.
            def first_contact(folder):
                for row in csv.DictReader((folder/'history.csv').open()):
                    keys=[k for k in row if k.endswith('_contact_nodes')]
                    if any(int(row[k])>0 for k in keys):return float(row['time_ns'])
                return None
            e0=first_contact(old/'cases'/name);e1=first_contact(dest)
            event_error=0. if e0 is None and e1 is None else abs(e1-e0)/e0 if e0 and e1 is not None else float('inf')
            assert area<=primary_tol and line<=primary_tol and event_error<=event_tol
            assert y['max_burgers_residual_b']<=criteria['burgers_residual_b'] and y['max_plastic_area_closure_absolute']<1e-12
            checks.append(dict(case_id=name,status='PASS',initial_source_code='BYTE_IDENTICAL',graph_and_crystal='EXACT',
                swept_area_relative=area,line_relative=line,event_relative=event_error,steps_old=x['steps'],steps_replay=y['steps'],
                plastic_closure=y['max_plastic_area_closure_absolute']))
    strict_status='PASS' if all(c['status']=='PASS' for c in strict) else 'FAIL'
    (new/'strict_full_field_replay.json').write_text(json.dumps(dict(status=strict_status,checks=strict,
        conclusion='No 1e-10 full force/velocity replay qualification. The native floating difference source is not independently isolated. Historical assertion failure is retained.'),indent=2))
    (new/'result.json').write_text(json.dumps(dict(status='PASS',gate_scope='Original frozen G0 scientific observables/invariants; NOT a bitwise or 1e-10 full-field qualification',
        A0='PASS',contact_native='PASS',diagnostic_native='PASS',source_byte_identity='PASS',four_native_trajectory_replays=checks,
        frozen_criteria_source=str(old/'frozen_criteria.json'),primary_threshold=primary_tol,event_threshold=event_tol,
        auxiliary_strict_full_field_replay=strict_status,auxiliary_evidence='strict_full_field_replay.json',
        preserved_failure='Parent g0_regression.log and its original 1e-10 assertion remain unchanged'),indent=2))
    print('G0_ORIGINAL_CONTRACT_PASS; AUXILIARY_STRICT_FULL_FIELD',strict_status)


if __name__=='__main__':main()
