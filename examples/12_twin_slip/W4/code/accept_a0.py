"""Freeze A0 only after actual configuration and native regressions pass."""
import argparse,json,shutil
from pathlib import Path
from configuration import save_json


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();folder=a.root/'A0'
    receipts=[json.loads((folder/n).read_text()) for n in ('contracts_receipt.json','native_v2_receipt.json','upstream_units_receipt.json','upstream_FR_receipt.json')]
    assert all(r['exit_code']==0 and not r['timed_out'] for r in receipts)
    native=json.loads((folder/'native_v2/result.json').read_text());assert native['status']=='PASS'
    for f,count in [('upstream_units.log',12),('upstream_FR.log',2)]:
        s=(folder/f).read_text();assert '100% tests passed' in s and f'out of {count}' in s
    snapshot=folder/'production_snapshot';snapshot.mkdir(exist_ok=True)
    for name in ('configuration.py','geometry.py','native.py'):
        shutil.copy2(a.root/'code'/name,snapshot/name)
    save_json(folder/'gate_result.json',dict(gate='A0',status='PASS',pure_contracts=5,native_response_checks=len(native['checks']),
        native_pin_cases=3,upstream_unit_cases=12,upstream_FR_cases=2,physics_execution_unlocked=True,
        no_VCoNi_material_qualification=True))
    print('GATE_A0_PASS')


if __name__=='__main__':main()
