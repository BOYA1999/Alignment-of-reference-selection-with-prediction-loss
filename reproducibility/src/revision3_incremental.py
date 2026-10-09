import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'artifacts/revision_003/incremental'
KEYS = ['endpoint', 'anchor', 'mechanism', 'model']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bootstrap_auc(score, label, samples):
    values, positive = np.asarray(score)[samples], np.asarray(label)[samples]
    n = positive.sum(axis=1)
    valid = (n > 0) & (n < values.shape[1])
    result = np.full(len(samples), np.nan)
    ranks = rankdata(values, axis=1)
    result[valid] = ((ranks * positive).sum(axis=1)[valid] - n[valid] * (n[valid] + 1) / 2) / (n[valid] * (values.shape[1] - n[valid]))
    return result


def nullable(value):
    return float(value) if np.isfinite(value) else None


def main():
    OUT.mkdir(parents=True, exist_ok=False)
    paths = [ROOT / f'artifacts/revision_002/{p}/estimates.csv' for p in ['bridge', 'original']]
    paths += [ROOT / 'artifacts/revision_003/protocol.md', Path(__file__)]
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    diagnostics, choices, summaries, auc_inputs = [], [], [], []
    for phase, path in zip(['bridge', 'original'], paths):
        frame = pd.read_csv(path, float_precision='round_trip')
        aux = frame[frame.method.eq('aipw') & frame.propensity.eq('low') & frame.loss_variant.isin(['low', 'drop'])]
        complete = frame[frame.method.eq('complete_case')].set_index(KEYS + ['seed'])
        for key, group in aux.groupby(KEYS):
            ident = {'phase': phase, **dict(zip(KEYS, key))}
            variants = {v: g.sort_values('seed').set_index('seed') for v, g in group.groupby('loss_variant')}
            low, drop = variants['low'], variants['drop']
            assert len(low) == len(drop) == 100 and low.index.equals(drop.index)
            assert np.array_equal(low.truth, drop.truth)
            cc = complete.loc[key].reindex(low.index)
            assert len(cc) == 100 and cc.estimate.notna().all()
            choice_drop = drop.top1_abs_share.to_numpy() < low.top1_abs_share.to_numpy()
            est = np.where(choice_drop, drop.estimate, low.estimate)
            truth = low.truth.to_numpy()
            for pos, seed in enumerate(low.index):
                choices.append({**ident, 'seed': int(seed), 'truth': truth[pos], 'choice_drop': bool(choice_drop[pos]),
                                'low': low.estimate.iloc[pos], 'drop': drop.estimate.iloc[pos], 'choice': est[pos],
                                'complete_case': cc.estimate.iloc[pos], 'share_low': low.top1_abs_share.iloc[pos],
                                'share_drop': drop.top1_abs_share.iloc[pos]})
            for name, values in [('low', low.estimate), ('drop', drop.estimate), ('choice', est), ('complete_case', cc.estimate)]:
                error = np.asarray(values) - truth
                summaries.append({**ident, 'method': name, 'n': 100, 'mean_abs_error': np.abs(error).mean(),
                                  'median_abs_error': np.median(np.abs(error)), 'bias': error.mean(),
                                  'rmse': np.sqrt(np.mean(error**2)), 'choice_drop_fraction': choice_drop.mean()})
            for variant, values in variants.items():
                target = values.abs_error.to_numpy()
                harmful = target > np.abs(cc.estimate.to_numpy() - truth)
                scores = {'concentration': values.top1_abs_share.to_numpy(),
                          'max_weight': values.max_normalized_weight.to_numpy(),
                          'inverse_ess': -values.ess_fraction.to_numpy(),
                          'smd': values.max_abs_smd.to_numpy(), 'clipping': values.clip_fraction.to_numpy()}
                auc_inputs.append((phase, variant, harmful, scores))
                for name, score in scores.items():
                    rho = spearmanr(score, target).statistic if np.ptp(score) > 0 and np.ptp(target) > 0 else np.nan
                    auc = roc_auc_score(harmful, score) if len(np.unique(harmful)) == 2 else np.nan
                    diagnostics.append({**ident, 'variant': variant, 'diagnostic': name, 'n': 100,
                                        'harmful_count': int(harmful.sum()), 'spearman': rho, 'auc': auc})
    c, d, s = pd.DataFrame(choices), pd.DataFrame(diagnostics), pd.DataFrame(summaries)
    c.to_csv(OUT / 'choice_rows.csv', index=False, float_format='%.17g')
    d.to_csv(OUT / 'diagnostic_by_setting.csv', index=False, float_format='%.17g')
    s.to_csv(OUT / 'choice_by_setting.csv', index=False, float_format='%.17g')
    report = {}
    rng = np.random.default_rng(93001)
    samples = rng.integers(0, 100, (2000, 100))
    for phase, g in c.groupby('phase'):
        groups = [z.sort_values('seed') for _, z in g.groupby(KEYS)]
        assert len(groups) == 16
        error = {m: np.stack([np.abs(z[m] - z.truth) for z in groups]) for m in ['low', 'drop', 'choice', 'complete_case']}
        delta = error['choice'] - error['low']
        boot = delta.mean(axis=0)[samples].mean(axis=1)
        by_diagnostic = d[d.phase.eq(phase)].groupby(['variant', 'diagnostic'])
        rank = [{ 'variant': v, 'diagnostic': name, 'settings': len(z), 'auc_defined': int(z.auc.notna().sum()),
                  'median_auc': nullable(z.auc.median()), 'median_spearman': nullable(z.spearman.median()),
                  'positive_spearman_settings': int(z.spearman.gt(0).sum())} for (v, name), z in by_diagnostic]
        report[phase] = {'settings': 16, 'masks_per_setting': 100,
                         'mean_abs_errors': {m: float(e.mean()) for m, e in error.items()},
                         'oracle_min_mean_abs_error': float(np.minimum(error['low'], error['drop']).mean()),
                         'choice_minus_low_mean_abs_error': float(delta.mean()),
                         'conditional_bootstrap_95ci': np.quantile(boot, [.025, .975]).tolist(),
                         'choice_better_low_settings': int((delta.mean(axis=1) < 0).sum()),
                         'choice_drop_fraction': float(g.choice_drop.mean()),
                         'paired_improved_fraction': float((delta < 0).mean()), 'paired_worsened_fraction': float((delta > 0).mean()),
                         'diagnostic_summary': rank}
        contrasts = []
        for variant in ['low', 'drop']:
            cells = [x for x in auc_inputs if x[0] == phase and x[1] == variant]
            base_auc = [bootstrap_auc(scores['concentration'], label, samples) for _, _, label, scores in cells]
            for baseline in ['max_weight', 'inverse_ess', 'smd', 'clipping']:
                boot_delta = np.nanmean(np.stack([a - bootstrap_auc(scores[baseline], label, samples)
                                                 for a, (_, _, label, scores) in zip(base_auc, cells)]), axis=0)
                z = d[d.phase.eq(phase) & d.variant.eq(variant)]
                primary = z[z.diagnostic.eq('concentration')].set_index(KEYS).auc
                ordinary = z[z.diagnostic.eq(baseline)].set_index(KEYS).auc
                contrasts.append({'variant': variant, 'baseline': baseline,
                                  'mean_auc_difference': nullable((primary - ordinary).mean()),
                                  'conditional_95ci': np.nanquantile(boot_delta, [.025, .975]).tolist(),
                                  'defined_bootstrap_replicates': int(np.isfinite(boot_delta).sum())})
        report[phase]['auc_contrasts'] = contrasts
    (OUT / 'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    assert all(sha(ROOT / p) == h for p, h in hashes.items())
    (OUT / 'manifest.json').write_text(json.dumps({'input_hashes': hashes, 'output_hashes': {p.name: sha(p) for p in OUT.iterdir()},
                                                  'scope': 'Post-result diagnostic test, no retraining or independent confirmation'}, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({p: {k: v for k, v in x.items() if k != 'diagnostic_summary'} for p, x in report.items()}, indent=2))


if __name__ == '__main__':
    main()
