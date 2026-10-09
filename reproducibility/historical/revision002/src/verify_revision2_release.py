import csv
import hashlib
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KEYS = ['endpoint','anchor','mechanism','model','method','propensity','loss_variant']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    assert __debug__, 'Assertions must be enabled'
    entries = [line.split('  ',1) for line in (ROOT/'MANIFEST.sha256').read_text().splitlines()]
    manifest = {name:h for h,name in entries}
    assert len(entries)==len(manifest)
    files = {p.relative_to(ROOT).as_posix():p for p in ROOT.rglob('*') if p.is_file() and p.name!='MANIFEST.sha256'}
    assert set(files)==set(manifest), 'Package inventory differs from manifest; use a fresh extraction'
    forbidden_columns = {'smiles','canonical_smiles','scaffold_group','ocnt_id','molecule_id'}
    for name,path in files.items():
        assert not Path(name).is_absolute() and '..' not in Path(name).parts
        assert sha(path)==manifest[name], name
        assert path.suffix not in {'.pyc','.npz','.zip','.docx'}, name
        assert not any(p in {'row_vectors','__pycache__','.venv','input','backups'} for p in path.parts), name
        if path.suffix=='.parquet':
            assert name=='artifacts/revision_002/alignment/replicates.parquet'
            assert not forbidden_columns.intersection(pd.read_parquet(path).columns.str.lower())
        if path.suffix=='.csv':
            with path.open(encoding='utf-8-sig',newline='') as stream:
                assert not forbidden_columns.intersection(x.lower() for x in next(csv.reader(stream))), name
        if path.suffix in {'.json','.md','.txt','.py','.csv','.svg','.bib'} and path!=Path(__file__).resolve():
            content = path.read_text(encoding='utf-8')
            assert not re.search(r'(?<![A-Za-z0-9])[A-Za-z]:[\\/]',content), name
            assert not re.search(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----',content), name
            assert not re.search(r'(?i:api_key)\s*=\s*[\"\x27]?[A-Za-z0-9_-]{16,}',content), name
    contract = json.loads((ROOT/'release_contract.json').read_text())
    for name,value in contract['protected_original_code_sha256'].items():
        assert sha(ROOT/name)==value
    for name,value in contract['retained_parent_files'].items():
        assert sha(ROOT/name)==value
    for row in json.loads((ROOT/'verification/distribution_redactions.json').read_text()):
        assert sha(ROOT/row['path'])==row['distributed_sha256']
    assert {p.name:sha(p) for p in (ROOT/'paper/figures').iterdir() if p.suffix!='.json'} == contract['figure_files']
    assert len(contract['figure_files'])==28
    base = ROOT/'artifacts/revision_002'
    counts = {}
    for phase,n,total in [('original',30400,304),('bridge',14400,144)]:
        data = pd.read_csv(base/phase/'estimates.csv',float_precision='round_trip')
        summary = pd.read_csv(base/phase/'summary.csv',float_precision='round_trip')
        assert len(data)==n and len(summary)==total
        grouped = data.groupby(KEYS)
        for row in summary.itertuples():
            group = grouped.get_group(tuple(getattr(row,k) for k in KEYS))
            assert len(group)==100 and group.seed.nunique()==100
            err = group.estimate.to_numpy()-row.truth
            expected = {'signed_bias':err.mean(),'rmse':np.sqrt(np.mean(err**2)),
                        'median_abs_error':np.median(abs(err)),'mc_sd':np.std(err,ddof=1)}
            assert all(np.isclose(v,getattr(row,k),rtol=1e-10,atol=1e-12) for k,v in expected.items())
        counts[phase]={'estimates':n,'summaries':total}
    raw = pd.read_csv(base/'original/estimates.csv',float_precision='round_trip')
    raw = raw[raw.method.eq('aipw')&raw.loss_variant.eq('raw')]
    assert raw.groupby('propensity').size().eq(1600).all() and raw.largest_source_row.eq(2945).all()
    assert raw[raw.propensity.eq('low')].groupby('case_R').size().to_dict()=={0:820,1:780}
    rep = pd.read_parquet(base/'alignment/replicates.parquet')
    summary = pd.read_csv(base/'alignment/summary.csv',float_precision='round_trip')
    grouped = rep.groupby(['endpoint','anchor','model','pairing'])
    assert len(rep)==320000 and len(summary)==192
    for row in summary.itertuples():
        part = grouped.get_group((row.endpoint,row.anchor,row.model,row.pairing))
        err = part[row.method].to_numpy()-row.truth
        assert len(err)==10000
        cc = np.median(abs(part.complete_case-row.truth))
        assert np.isclose(np.median(abs(err)),row.median_abs_error,rtol=1e-11,atol=1e-12)
        assert np.isclose(np.sqrt(np.mean(err**2)),row.rmse,rtol=1e-11,atol=1e-12)
        assert np.isclose(np.median(abs(err))/cc,row.median_abs_error_ratio_vs_cc,rtol=1e-11)
    design = pd.read_csv(base/'alignment/design.csv')
    assert design.groupby(['endpoint','anchor','model']).sorted_probability_sha256.nunique().eq(1).all()
    trace = pd.read_csv(base/'case_trace.csv')
    assert len(trace)==64 and set(trace.case_R)=={0,1}
    contribution = ((1-trace.case_R/trace.case_p)*trace.case_m+trace.case_R*trace.case_loss/trace.case_p)/trace.n
    assert np.allclose(contribution,trace.case_pseudo_contribution,rtol=1e-10,atol=1e-12)
    text = (ROOT/'paper/manuscript/JPA_original_article.md').read_text(encoding='utf-8')
    assert all(s in text for s in ['11/16','10/16','1/16','1,600','Table 2'])
    ga = (ROOT/'paper/figures/Graphical_Abstract.svg').read_text(encoding='utf-8')
    assert 'Single-seed Ipc diagnostic' not in ga and '323.07' in ga and '11/16' in ga
    report = {'status':'PASS_LOCAL_PACKAGE_CHECKS','manifest_entries':len(manifest),'figure_files':28,
              'fitted_summaries_checked':448,'alignment_summaries_checked':192,'case_trace_rows':64,
              'counts':counts,'protected_code_unchanged':True,'retained_parent_files':len(contract['retained_parent_files']),
              'scope':'File integrity, distribution scan, aggregate recomputation and saved case arithmetic; not a full predictor/auxiliary retraining or external scientific validation'}
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
