import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
KEYS = ['endpoint', 'anchor', 'mechanism', 'model']
PATHS = [ROOT / 'artifacts/revision_002/bridge' / name for name in ['summary.csv', 'estimates.csv']]
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
before = {str(p.relative_to(ROOT)): sha(p) for p in PATHS}
s, x = [pd.read_csv(p, float_precision='round_trip') for p in PATHS]
for frame in [s, x]:
    frame['label'] = frame.method + '_' + frame.propensity + '_' + frame.loss_variant
assert len(s) == 144 and len(x) == 14400
assert not s.duplicated(KEYS + ['label']).any()
assert not x.duplicated(KEYS + ['label', 'seed']).any()
assert s.replicates.eq(100).all()
med = s.pivot(index=KEYS, columns='label', values='median_abs_error')
avg = s.pivot(index=KEYS, columns='label', values='mean_abs_error')
assert med.shape == avg.shape == (16, 9)
recomputed = x.groupby(KEYS + ['label']).abs_error.agg(['mean', 'median'])
matched = s.set_index(KEYS + ['label']).join(recomputed)
maximum_difference = float(max((matched.mean_abs_error - matched['mean']).abs().max(),
                               (matched.median_abs_error - matched['median']).abs().max()))
assert maximum_difference < 1e-12
wide = med.add_suffix('_median_abs_error').join(avg.add_suffix('_mean_abs_error'))
pairs = [(a, b) for b in ['complete_case_none_none', 'hajek_low_none']
         for a in ['aipw_low_low', 'aipw_low_drop']]
pairs += [('hajek_low_none', 'complete_case_none_none'), ('ht_low_none', 'complete_case_none_none')]
comparison_summaries = []
for a, b in pairs:
    ratio = med[a] / med[b]
    wide[a + '_over_' + b + '_median_error_ratio'] = ratio
    comparison_summaries.append({'numerator': a, 'denominator': b, 'point_wins': int((ratio < 1).sum()),
                                'ratio_min': float(ratio.min()), 'ratio_max': float(ratio.max()),
                                'median_setting_ratio': float(ratio.median())})
wide.to_csv(OUT / 'fitted_comparator_16settings.csv', float_format='%.17g')
indices = np.random.default_rng(923).integers(0, 100, size=(10000, 100))
rows = []
for key, group in x.groupby(KEYS):
    p = group.pivot(index='seed', columns='label', values='abs_error').sort_index()
    assert p.index.tolist() == list(range(71001, 71101)) and p.shape == (100, 9)
    for a, b in pairs[:4]:
        ratios = np.median(p[a].to_numpy()[indices], axis=1) / np.median(p[b].to_numpy()[indices], axis=1)
        lo, hi = np.quantile(ratios, [.025, .975], method='linear')
        rows.append({**dict(zip(KEYS, key)), 'numerator': a, 'denominator': b,
                     'ratio': float(med.loc[key, a] / med.loc[key, b]),
                     'conditional_mc_lower': float(lo), 'conditional_mc_upper': float(hi)})
b = pd.DataFrame(rows)
b.to_csv(OUT / 'paired_mc_ratio_intervals.csv', index=False, float_format='%.17g')
counts = []
for (a, denominator), g in b.groupby(['numerator', 'denominator']):
    counts.append({'numerator': a, 'denominator': denominator, 'settings': len(g),
                   'point_wins': int(g.ratio.lt(1).sum()), 'upper_below_one': int(g.conditional_mc_upper.lt(1).sum()),
                   'interval_contains_one': int((g.conditional_mc_lower.le(1) & g.conditional_mc_upper.ge(1)).sum())})
assert before == {str(p.relative_to(ROOT)): sha(p) for p in PATHS}
report = {'scope': 'Review-only reaggregation and conditional Monte Carlo precision; no training, manuscript edit, new confirmatory experiment or population inference.',
          'settings': 16, 'masks_per_setting': 100, 'bootstrap_replicates': 10000, 'bootstrap_seed': 923,
          'resampling': 'Same paired seed index matrix applied to all settings and methods; seeds sorted 71001 through 71100; percentile intervals, linear quantiles.',
          'index_matrix_sha256': hashlib.sha256(indices.tobytes()).hexdigest(),
          'aggregation': 'Equal-setting arithmetic mean of each setting mean absolute error; not median of all errors, not molecule weighting. Equivalent to pooling 1600 estimates per method only because all 16 settings have 100 masks.',
          'method_labels': {'hajek_low_none': 'Hajek self-normalized IPW with fitted cross-fitted low-dimensional logistic propensity; no loss model.',
                            'ht_low_none': 'Horvitz-Thompson with fitted cross-fitted low-dimensional propensity; denominator is full target cohort size.',
                            'aipw_low_low': 'AIPW with fitted low-dimensional propensity and low-dimensional conditional-loss ridge.',
                            'aipw_low_drop': 'AIPW with same fitted low-dimensional propensity and full conditional-loss ridge after raw Ipc removal.',
                            'complete_case_none_none': 'Observed-case mean loss.'},
          'equal_setting_mean_absolute_errors': {k: float(v) for k, v in avg.mean().items()},
          'comparisons': comparison_summaries, 'bootstrap_counts': counts,
          'summary_to_estimates_max_difference': maximum_difference, 'inputs_unchanged': True,
          'input_sha256': before, 'script_sha256': sha(Path(__file__)),
          'versions': {'numpy': np.__version__, 'pandas': pd.__version__}}
(OUT / 'existing_comparison_reanalysis.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps({'status': 'PASS_REVIEW_ONLY_REAGGREGATION', 'summary_max_difference': maximum_difference,
                  'bootstrap_counts': counts, 'input_files_unchanged': True}, indent=2))
