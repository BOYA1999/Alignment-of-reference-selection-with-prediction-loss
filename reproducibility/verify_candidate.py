import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
KEYS = ['endpoint', 'anchor', 'mechanism', 'model']
read = lambda p: pd.read_csv(ROOT / p, float_precision='round_trip')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
close = lambda a, b: np.allclose(a, b, rtol=1e-10, atol=1e-12, equal_nan=True)


def check_files(manifest):
    files = {p.relative_to(ROOT).as_posix(): p for p in ROOT.rglob('*') if p.is_file() and p != ROOT / 'MANIFEST.json'}
    assert set(files) == set(manifest['files']), 'Inventory differs; use a fresh extraction'
    prohibited = {'smiles', 'canonical_smiles', 'scaffold_group', 'ocnt_id', 'molecule_id'}
    findings = []
    for name, path in files.items():
        expected = manifest['files'][name]
        assert path.stat().st_size == expected['bytes'] and sha(path) == expected['sha256'], name
        assert not path.is_symlink() and not Path(name).is_absolute() and '..' not in Path(name).parts and ':' not in name, name
        assert path.suffix.lower() not in {'.zip', '.npz', '.pyc', '.docx', '.sdf', '.smi'}, name
        assert not {'row_vectors', '__pycache__', '.venv', 'raw', 'input', 'backups'} & set(Path(name).parts), name
        if path.suffix == '.csv':
            with path.open(encoding='utf-8-sig', newline='') as stream:
                assert not prohibited.intersection(v.lower() for v in next(csv.reader(stream))), name
        if path.suffix == '.parquet':
            assert name == 'historical/revision002/artifacts/revision_002/alignment/replicates.parquet', name
        if path.suffix in {'.json', '.md', '.txt', '.py', '.csv', '.svg', '.bib', '.cff'}:
            content = path.read_text(encoding='utf-8-sig')
            for tag, pattern in [('local_drive_path', r'(?<![A-Za-z0-9])[A-Za-z]:[\\/]'),
                                 ('private_key', r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
                                 ('credential_assignment', r'(?i:api_key)\s*=\s*["\x27]?[A-Za-z0-9_-]{16,}')]:
                if re.search(pattern, content):
                    findings.append({'path': name, 'finding': tag})
    assert not findings, findings
    return len(files)


def main():
    assert __debug__, 'Do not use python -O'
    manifest = json.loads((ROOT / 'MANIFEST.json').read_text(encoding='utf-8'))
    manifest_hash = sha(ROOT / 'MANIFEST.json')
    count = check_files(manifest)
    inherited = {}
    for prefix, name in [('historical/revision002', 'MANIFEST.sha256'), ('addenda/multisplit', 'manifest-sha256.txt')]:
        expected = {}
        for line in (ROOT / prefix / name).read_text().splitlines():
            value, rel = line.split('  ', 1)
            assert rel not in expected
            expected[rel] = value
        actual = {p.relative_to(ROOT / prefix).as_posix(): sha(p) for p in (ROOT / prefix).rglob('*') if p.is_file() and p != ROOT / prefix / name}
        assert actual == expected
        inherited[prefix] = {'manifest_payload_files': len(actual), 'hashes_match': True}
    script = ROOT / 'addenda/multisplit/verify_package.py'
    result = subprocess.run([sys.executable, '-B', str(script)], cwd=script.parent, capture_output=True, text=True, encoding='utf-8')
    assert result.returncode == 0, (result.stdout, result.stderr)
    inherited['multisplit_saved_result_replay'] = result.stdout.strip()
    source = {phase: read(f'artifacts/revision_002/{phase}/estimates.csv') for phase in ['bridge', 'original']}
    fitted_count = 0
    for phase, estimates in source.items():
        assert close(estimates.abs_error, abs(estimates.estimate - estimates.truth))
        summary = read(f'artifacts/revision_002/{phase}/summary.csv')
        key = KEYS + ['method', 'propensity', 'loss_variant']
        grouped = estimates.groupby(key)
        for row in summary.itertuples():
            group = grouped.get_group(tuple(getattr(row, k) for k in key))
            assert len(group) == 100 and group.seed.nunique() == 100
            err = group.estimate.to_numpy() - row.truth
            assert close([err.mean(), np.sqrt((err ** 2).mean()), np.median(abs(err)), np.std(err, ddof=1)],
                         [row.signed_bias, row.rmse, row.median_abs_error, row.mc_sd])
            fitted_count += 1
    assert fitted_count == 448
    trace = read('historical/revision002/artifacts/revision_002/case_trace.csv')
    assert len(trace) == 64 and set(trace.case_R) == {0, 1}
    contribution = ((1 - trace.case_R / trace.case_p) * trace.case_m + trace.case_R * trace.case_loss / trace.case_p) / trace.n
    assert close(contribution, trace.case_pseudo_contribution)
    choices = read('artifacts/revision_003/incremental/choice_rows.csv')
    settings = read('artifacts/revision_003/incremental/choice_by_setting.csv')
    report3 = json.loads((ROOT / 'artifacts/revision_003/incremental/summary.json').read_text())
    assert len(choices) == 3200 and len(settings) == 128 and not choices.duplicated(['phase'] + KEYS + ['seed']).any()
    rng3 = np.random.default_rng(93001).integers(0, 100, (2000, 100))
    utility, intervals = {}, {}
    for phase, groups in choices.groupby('phase'):
        deltas = []
        src = source[phase]
        for key, rows in groups.groupby(KEYS):
            rows = rows.sort_values('seed')
            subset = src.loc[(src[KEYS] == pd.Series(key, index=KEYS)).all(axis=1)]
            assert len(rows) == 100 and rows.seed.nunique() == 100
            for variant in ['low', 'drop']:
                original = subset[subset.method.eq('aipw') & subset.propensity.eq('low') & subset.loss_variant.eq(variant)].sort_values('seed')
                assert np.array_equal(rows.seed, original.seed) and close(rows.truth, original.truth)
                assert close(rows[variant], original.estimate) and close(rows['share_' + variant], original.top1_abs_share)
            cc = subset[subset.method.eq('complete_case')].sort_values('seed')
            assert close(rows.complete_case, cc.estimate)
            assert np.array_equal(rows.choice_drop, rows.share_drop < rows.share_low)
            assert close(rows.choice, np.where(rows.share_drop < rows.share_low, rows['drop'], rows.low))
            deltas.append((abs(rows.choice - rows.truth) - abs(rows.low - rows.truth)).to_numpy())
            for method in ['low', 'drop', 'choice', 'complete_case']:
                line = settings.loc[settings.phase.eq(phase) & settings.method.eq(method) & (settings[KEYS] == pd.Series(key, index=KEYS)).all(axis=1)].iloc[0]
                err = rows[method].to_numpy() - rows.truth.to_numpy()
                assert close([abs(err).mean(), np.median(abs(err)), err.mean(), np.sqrt((err ** 2).mean())],
                             [line.mean_abs_error, line.median_abs_error, line.bias, line.rmse])
        delta = np.stack(deltas)
        assert delta.shape == (16, 100)
        ci = np.quantile(delta.mean(axis=0)[rng3].mean(axis=1), [.025, .975])
        assert close(ci, report3[phase]['conditional_bootstrap_95ci'])
        assert close(delta.mean(), report3[phase]['choice_minus_low_mean_abs_error'])
        means = {m: float(abs(groups[m] - groups.truth).mean()) for m in ['low', 'drop', 'choice', 'complete_case']}
        assert close(list(means.values()), [report3[phase]['mean_abs_errors'][m] for m in means])
        utility[phase] = 100 * (means['choice'] / means['low'] - 1)
        table = read('paper/tables/Table_S18_utility_summary.csv').set_index('phase').loc[phase]
        assert close([utility[phase], *ci], [table.change_percent, table.conditional_mc_lower, table.conditional_mc_upper])
        display = pd.read_csv(ROOT / 'paper/tables/Table_S18_display.csv', dtype=str, keep_default_na=False)
        row = display.iloc[0 if phase == 'bridge' else 1]
        assert row['Always low'] == f"{means['low']:.6f}" and row['Choice'] == f"{means['choice']:.6f}"
        assert row['Change (%)'] == f"{utility[phase]:+.2f}"
        assert row['Difference (95% MC interval)'] == f"{delta.mean():.6f} [{ci[0]:.6f}, {ci[1]:.6f}]"
    bridge = source['bridge'].copy()
    bridge['label'] = bridge.method + '_' + bridge.propensity + '_' + bridge.loss_variant
    s19 = read('paper/tables/Table_S19_paired_mc_ratio_intervals.csv').set_index(KEYS + ['numerator', 'denominator'])
    comparator = read('paper/tables/Table_S19_fitted_comparators_16settings.csv')
    comparator_index = comparator.set_index(KEYS)
    rng19 = np.random.default_rng(923).integers(0, 100, (10000, 100))
    pairs = [(a, b) for b in ['complete_case_none_none', 'hajek_low_none'] for a in ['aipw_low_low', 'aipw_low_drop']]
    values, wins, upper, contains = {}, {}, {}, {}
    for key, group in bridge.groupby(KEYS):
        p = group.pivot(index='seed', columns='label', values='abs_error').sort_index()
        assert p.index.tolist() == list(range(71001, 71101)) and p.shape == (100, 9)
        for label in p.columns:
            assert close([p[label].median(), p[label].mean()], [comparator_index.loc[key, label + '_median_abs_error'], comparator_index.loc[key, label + '_mean_abs_error']])
        for a, b in pairs:
            ratio = np.median(p[a]) / np.median(p[b])
            ci = np.quantile(np.median(p[a].to_numpy()[rng19], axis=1) / np.median(p[b].to_numpy()[rng19], axis=1), [.025, .975])
            line = s19.loc[key + (a, b)]
            assert close([ratio, *ci], [line.ratio, line.conditional_mc_lower, line.conditional_mc_upper])
            values[key + (a, b)] = (ratio, *ci)
            label = a + '/' + b
            wins[label] = wins.get(label, 0) + int(ratio < 1)
            upper[label] = upper.get(label, 0) + int(ci[1] < 1)
            contains[label] = contains.get(label, 0) + int(ci[0] <= 1 <= ci[1])
    assert len(values) == 64
    assert wins['aipw_low_low/hajek_low_none'] == 3 and wins['aipw_low_drop/hajek_low_none'] == 1
    for row in read('paper/tables/Table_S19_mc_precision_summary.csv').itertuples():
        label = row.numerator + '/' + row.denominator
        assert (row.point_wins, row.upper_below_one, row.interval_contains_one) == (wins[label], upper[label], contains[label])
    display = pd.read_csv(ROOT / 'paper/tables/Table_S19_display.csv', dtype=str, keep_default_na=False)
    for i, row in enumerate(comparator.to_dict('records')):
        key = tuple(row[k] for k in KEYS)
        assert display.iloc[i]['Fitted Hajek error'] == f"{row['hajek_low_none_median_abs_error']:.6f}"
        for variant, column in [('low', 'Low AIPW / Hajek'), ('drop', 'Dropped Ipc / Hajek')]:
            ratio, lo, hi = values[key + ('aipw_low_' + variant, 'hajek_low_none')]
            assert display.iloc[i][column] == f'{ratio:.3f} [{lo:.3f}, {hi:.3f}]'
    copies = [('Table_S18_' + name + '.csv', 'artifacts/revision_003/incremental/' + src + '.csv') for name, src in [
        ('choice_rows', 'choice_rows'), ('choice_by_setting', 'choice_by_setting'),
        ('diagnostic_by_setting', 'diagnostic_by_setting'), ('utility_contrasts', 'verified_utility_contrasts')]]
    copies += [(name, 'paper/review/peer_review_20260923/' + src) for name, src in [
        ('Table_S19_paired_mc_ratio_intervals.csv', 'paired_mc_ratio_intervals.csv'),
        ('Table_S19_fitted_comparators_16settings.csv', 'fitted_comparator_16settings.csv')]]
    for name, origin in copies:
        pd.testing.assert_frame_equal(read('paper/tables/' + name), read(origin), check_exact=False, rtol=1e-12, atol=1e-14)
    diagnostics = pd.DataFrame([{'phase': phase, **row} for phase, data in report3.items() for row in data['diagnostic_summary']])
    pd.testing.assert_frame_equal(read('paper/tables/Table_S18_diagnostic_summary.csv'), diagnostics, check_exact=False)
    auc = pd.DataFrame([{'phase': phase, **{k: v for k, v in row.items() if k != 'conditional_95ci'},
                         'conditional_mc_lower': row['conditional_95ci'][0], 'conditional_mc_upper': row['conditional_95ci'][1]}
                        for phase, data in report3.items() for row in data['auc_contrasts']])
    pd.testing.assert_frame_equal(read('paper/tables/Table_S18_auc_contrasts.csv'), auc, check_exact=False)
    reanalysis = json.loads((ROOT / 'paper/review/peer_review_20260923/existing_comparison_reanalysis.json').read_text())
    pd.testing.assert_frame_equal(read('paper/tables/Table_S19_comparison_summary.csv'), pd.DataFrame(reanalysis['comparisons']), check_exact=False)
    assert len(list((ROOT / 'paper/tables').glob('Table_S18_*.csv'))) == 8
    assert len(list((ROOT / 'paper/tables').glob('Table_S19_*.csv'))) == 5
    assert check_files(manifest) == count and sha(ROOT / 'MANIFEST.json') == manifest_hash
    print(json.dumps({'status': 'PASS_LOCAL_CANDIDATE_SAVED_RESULT_QA', 'payload_files': count,
                      'exact_inventory_sizes_sha256': True, 'distribution_scan_findings': [], 'files_unchanged_after_checks': True,
                      'inherited_readonly_checks': inherited, 'S18_choice_rows': 3200, 'S18_setting_summaries': 128,
                      'fitted_summaries_recomputed': fitted_count, 'case_contributions_recomputed': 64,
                      'alignment_parquet_arithmetic_repeated': False,
                      'S18_choice_change_percent': utility, 'S19_paired_intervals': 64, 'S19_point_wins': wins,
                      'S18_S19_electronic_tables': 13, 'full_model_retraining': False, 'external_validation': False,
                      'versions': {'python': sys.version, 'numpy': np.__version__, 'pandas': pd.__version__},
                      'scope': 'Saved aggregate and choice arithmetic, conditional Monte Carlo intervals, exact file and limited pattern checks; no raw-data or full model replay'}, indent=2))


if __name__ == '__main__':
    main()
