import hashlib
import json
from pathlib import Path
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
C1 = ROOT / 'artifacts/revision_003/incremental'
REVIEW = ROOT / 'paper/review/peer_review_20260923'
paths = [C1 / 'summary.json', REVIEW / 'existing_comparison_reanalysis.json',
         REVIEW / 'fitted_comparator_16settings.csv', REVIEW / 'paired_mc_ratio_intervals.csv']
names = ['choice_by_setting', 'choice_rows', 'diagnostic_by_setting', 'verified_utility_contrasts']
paths += [C1 / (name + '.csv') for name in names]
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
before = {str(p.relative_to(ROOT)): sha(p) for p in paths}
read = lambda p: pd.read_csv(p, float_precision='round_trip')
write = lambda frame, name: frame.to_csv(OUT / name, index=False, float_format='%.17g')
c1 = json.loads(paths[0].read_text(encoding='utf-8'))
review = json.loads(paths[1].read_text(encoding='utf-8'))
utility, display, diagnostic, contrasts = [], [], [], []
for phase, label in [('bridge', 'Molecular weight'), ('original', 'Original masks')]:
    v = c1[phase]
    means = v['mean_abs_errors']
    change = 100 * (means['choice'] / means['low'] - 1)
    lo, hi = v['conditional_bootstrap_95ci']
    utility.append({'phase': phase, **means, 'change_percent': change,
                    'choice_minus_low': v['choice_minus_low_mean_abs_error'], 'conditional_mc_lower': lo,
                    'conditional_mc_upper': hi, 'settings': 16, 'masks_per_setting': 100,
                    'choice_better_low_settings': v['choice_better_low_settings'],
                    'choice_drop_fraction': v['choice_drop_fraction'],
                    'paired_improved_fraction': v['paired_improved_fraction'],
                    'paired_worsened_fraction': v['paired_worsened_fraction']})
    display.append({'Mask set': label, 'Always low': f"{means['low']:.6f}", 'Choice': f"{means['choice']:.6f}",
                    'Change (%)': f'{change:+.2f}',
                    'Difference (95% MC interval)': f"{v['choice_minus_low_mean_abs_error']:.6f} [{lo:.6f}, {hi:.6f}]",
                    'Improved settings': f"{v['choice_better_low_settings']}/16"})
    diagnostic += [{'phase': phase, **row} for row in v['diagnostic_summary']]
    contrasts += [{'phase': phase, **{k: val for k, val in row.items() if k != 'conditional_95ci'},
                    'conditional_mc_lower': row['conditional_95ci'][0],
                    'conditional_mc_upper': row['conditional_95ci'][1]} for row in v['auc_contrasts']]
for rows, suffix in [(utility, 'utility_summary'), (display, 'display'),
                     (diagnostic, 'diagnostic_summary'), (contrasts, 'auc_contrasts')]:
    write(pd.DataFrame(rows), 'Table_S18_' + suffix + '.csv')
for name, expected in zip(names, [128, 3200, 320, 32]):
    frame = read(C1 / (name + '.csv'))
    assert len(frame) == expected
    write(frame, 'Table_S18_' + name.replace('verified_', '') + '.csv')
comparators = read(REVIEW / 'fitted_comparator_16settings.csv')
intervals = read(REVIEW / 'paired_mc_ratio_intervals.csv')
assert len(comparators) == 16 and len(intervals) == 64
write(comparators, 'Table_S19_fitted_comparators_16settings.csv')
write(intervals, 'Table_S19_paired_mc_ratio_intervals.csv')
write(pd.DataFrame(review['bootstrap_counts']), 'Table_S19_mc_precision_summary.csv')
write(pd.DataFrame(review['comparisons']), 'Table_S19_comparison_summary.csv')
display = []
keys = ['endpoint', 'anchor', 'mechanism', 'model']
ci = intervals.set_index(keys + ['numerator', 'denominator'])
for row in comparators.to_dict('records'):
    key = tuple(row[k] for k in keys)
    label = ' / '.join([row['endpoint'].removeprefix('LOG_'), row['anchor'].removeprefix('LOG_').replace('_ER', ''),
                        'MW+' if row['mechanism'].endswith('positive') else 'MW-',
                        'XGB' if row['model'] == 'ecfp_desc_xgb' else 'RF'])
    entry = {'Setting': label, 'Fitted Hajek error': f"{row['hajek_low_none_median_abs_error']:.6f}"}
    for variant, column in [('low', 'Low AIPW / Hajek'), ('drop', 'Dropped Ipc / Hajek')]:
        r = ci.loc[key + ('aipw_low_' + variant, 'hajek_low_none')]
        entry[column] = f"{r.ratio:.3f} [{r.conditional_mc_lower:.3f}, {r.conditional_mc_upper:.3f}]"
    display.append(entry)
write(pd.DataFrame(display), 'Table_S19_display.csv')
assert before == {str(p.relative_to(ROOT)): sha(p) for p in paths}
print(json.dumps({'status': 'PASS_EXISTING_EVIDENCE_EXPORT', 'source_files_unchanged': True,
                  'S18_display_rows': 2, 'S19_display_rows': 16, 'S19_interval_rows': 64,
                  'source_sha256': before}, indent=2))
