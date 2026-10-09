import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from urllib.request import urlopen
from rdkit import Chem, rdBase
from rdkit.Chem import GraphDescriptors, Descriptors

ROOT = Path(__file__).resolve().parents[1]
URL = 'https://raw.githubusercontent.com/molecularinformatics/Computational-ADME/b00df003de117ce9e5b381afd886095c5f2af2d5/ADME_public_set_3521.csv'
SHA = '2cfabc2667740c224487876c33b23124159ef43294e0f9e4d926cb6276c95a3b'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,help='Pinned official CSV; if omitted fetch the public source in memory')
    args = parser.parse_args()
    data = args.source.read_bytes() if args.source else urlopen(URL,timeout=60).read()
    assert hashlib.sha256(data).hexdigest()==SHA
    rows = list(csv.DictReader(io.StringIO(data.decode('utf-8-sig'))))
    case = rows[2945]
    assert case['Internal ID']=='Mol2946'
    mol = Chem.MolFromSmiles(case['SMILES'])
    ipc,avg = GraphDescriptors.Ipc(mol),GraphDescriptors.AvgIpc(mol)
    assert abs(ipc/1.9516295e17-1)<1e-6
    trace=list(csv.DictReader((ROOT/'artifacts/revision_002/case_trace.csv').open(encoding='utf-8')))
    maximum=0.
    for row in trace:
        r,p,m,loss,n = [float(row[k]) for k in ['case_R','case_p','case_m','case_loss','n']]
        pseudo=(1-r/p)*m+r*loss/p
        observed=float(row['case_pseudo_contribution'])*n
        delta=abs(pseudo-observed)/max(1.,abs(observed))
        maximum=max(maximum,delta)
        assert delta<1e-10
    print(json.dumps(dict(source_url=URL,source_sha256=SHA,source_zero_based_row=2945,internal_id=case['Internal ID'],
        rdkit=rdBase.rdkitVersion,atoms=mol.GetNumAtoms(),heavy_atoms=mol.GetNumHeavyAtoms(),rings=Descriptors.RingCount(mol),
        molecular_weight=Descriptors.MolWt(mol),Ipc=ipc,AvgIpc=avg,trace_rows=len(trace),max_relative_identity_difference=maximum,
        scope='Recomputes the descriptor from pinned public structure and replays saved single-row AIPW contributions; does not refit models'),indent=2))


if __name__=='__main__':
    main()
