import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
result = subprocess.run([sys.executable, '-B', str(ROOT / 'verify_candidate.py')], cwd=ROOT,
                        capture_output=True, text=True, encoding='utf-8')
assert result.returncode == 0, (result.stdout, result.stderr)
print(result.stdout.strip())
read = lambda name: pd.read_csv(ROOT / name, float_precision='round_trip')
keys = ['endpoint', 'anchor', 'mechanism', 'model', 'numerator', 'denominator']
figure = read('paper/source_data/Figure_3_fitted_hajek.csv').sort_values(keys).reset_index(drop=True)
table = read('paper/tables/Table_S19_paired_mc_ratio_intervals.csv')
table = table[table.denominator.eq('hajek_low_none')].sort_values(keys).reset_index(drop=True)
pd.testing.assert_frame_equal(figure, table, check_exact=False, rtol=1e-14, atol=1e-15)
assert len(figure) == 32 and not figure.duplicated(keys).any()
assert not (figure.conditional_mc_upper < 1).any()
assert (figure.conditional_mc_lower > 1).sum() == 2
multi = read('addenda/multisplit/results/estimates.csv')
unit = ['split_seed', 'endpoint', 'anchor', 'mechanism', 'mask_seed', 'model']
assert len(multi) == 80000 and len(multi[unit].drop_duplicates()) == 16000
expected = {('complete_case', 'none'), ('hajek', 'none'), ('aipw', 'low'), ('aipw', 'raw'), ('aipw', 'drop')}
assert set(multi.groupby(['method', 'loss_variant']).size().items()) == {(k, 16000) for k in expected}
assert multi.groupby(unit).size().eq(5).all()
assert not multi.duplicated(unit + ['method', 'loss_variant']).any()
original = read('paper/tables/Table_1_cohorts_models_and_availability.csv')
pxr = pd.read_csv(ROOT / 'paper/tables/Table_S5C_pxr_cohort_display.csv', dtype=str)
assert len(pxr) == 2 and pxr.iloc[0]['Value'] == '10,747'
assert np.isclose(float(pxr.iloc[1]['Value']), 2722 / 10747)
assert original[original.section.eq('PXR boundary')].shape[0] == 2
for ext in ['png', 'pdf', 'svg']:
    a = ROOT / f'paper/figures/Figure_S5_fitted_oracle.{ext}'
    b = ROOT / f'paper/figures/Figure_3_fitted_bridge.{ext}'
    assert a.read_bytes() == b.read_bytes()
manifest = json.loads((ROOT / 'MANIFEST.json').read_text())
assert set(manifest['files']) | {'MANIFEST.json'} == {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()}
for name, entry in manifest['files'].items():
    data = (ROOT / name).read_bytes()
    assert len(data) == entry['bytes'] and hashlib.sha256(data).hexdigest() == entry['sha256'], name
print(json.dumps({'status': 'PASS_MINOR_REVISION_SAVED_RESULTS_AND_FIGURE_INPUTS',
                  'figure_points_and_intervals': 32, 'multisplit_units': 16000,
                  'multisplit_records_per_unit': 5, 'multisplit_estimates': 80000,
                  'oracle_figure_bytes_preserved': True, 'PXR_rows_preserved': 2,
                  'new_model_training': False, 'public_access_created': False}, indent=2))
