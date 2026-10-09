import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'artifacts/revision_003/incremental'
KEYS = ['endpoint', 'anchor', 'mechanism', 'model']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pair_auc(label, score):
    a, b = score[label], score[~label]
    return float(((a[:, None] > b).sum() + .5 * (a[:, None] == b).sum()) / (len(a) * len(b)))


def main():
    manifest = json.loads((OUT / 'manifest.json').read_text())
    assert all(digest(ROOT / p) == h for p, h in manifest['input_hashes'].items())
    assert all(digest(OUT / p) == h for p, h in manifest['output_hashes'].items())
    protected = ROOT / 'artifacts/confirmatory/confirm_001/full_recovery_001'
    protected_hashes = {p.name: digest(p) for p in protected.iterdir()}
    assert len(protected_hashes) == 19
    c = pd.read_csv(OUT / 'choice_rows.csv', float_precision='round_trip')
    d = pd.read_csv(OUT / 'diagnostic_by_setting.csv', float_precision='round_trip')
    s = pd.read_csv(OUT / 'choice_by_setting.csv', float_precision='round_trip')
    report = json.loads((OUT / 'summary.json').read_text())
    assert len(c) == 3200 and len(d) == 320 and len(s) == 128
    assert not c.duplicated(['phase'] + KEYS + ['seed']).any()
    assert np.array_equal(c.choice_drop, c.share_drop < c.share_low)
    assert np.array_equal(c.choice, np.where(c.share_drop < c.share_low, c['drop'], c.low))
    predictions = pd.read_parquet(protected / 'predictions.parquet')
    models = ['ecfp_rf', 'ecfp_desc_xgb']
    losses = {ep: np.column_stack([g[g.model.eq(m)].sort_values('source_row_index').loss for m in models])
              for ep, g in predictions.groupby('endpoint')}
    counts, max_difference, contrasts = 0, 0., []
    bootstrap = np.random.default_rng(93001).integers(0, 100, (2000, 100))
    for phase, phase_rows in c.groupby('phase'):
        source = pd.read_csv(ROOT / f'artifacts/revision_002/{phase}/estimates.csv', float_precision='round_trip')
        deltas = []
        for key, rows in phase_rows.groupby(KEYS):
            rows = rows.sort_values('seed')
            assert len(rows) == 100
            delta = abs(rows.choice - rows.truth) - abs(rows.low - rows.truth)
            deltas.append(delta.to_numpy())
            med_low, med_choice = np.median(abs(rows.low - rows.truth)), np.median(abs(rows.choice - rows.truth))
            contrasts.append({'phase': phase, **dict(zip(KEYS, key)), 'choice_minus_low_mean_abs_error': delta.mean(),
                              'choice_to_low_median_error_ratio': med_choice / med_low})
        matrix = np.stack(deltas)
        boot_values = np.array([matrix[:, take].mean() for take in bootstrap])
        assert abs(matrix.mean() - report[phase]['choice_minus_low_mean_abs_error']) < 1e-14
        assert np.allclose(np.quantile(boot_values, [.025, .975]), report[phase]['conditional_bootstrap_95ci'], rtol=0, atol=1e-14)
        for key, rows in phase_rows.groupby(['endpoint', 'anchor', 'mechanism', 'seed']):
            endpoint, anchor, mechanism, seed = key
            vector = ROOT / f'artifacts/revision_002/{phase}/row_vectors/{endpoint}/{anchor}_{mechanism}/{seed}.npz'
            with np.load(vector) as a:
                r, p = a['R'], a['p_low']
                observed_loss = np.where(r[:, None].astype(bool), losses[endpoint], 0.)
                changed_loss = np.where(r[:, None].astype(bool), losses[endpoint], 12345.)
                w = r[:, None] / p[:, None]
                for variant in ['low', 'drop']:
                    m = a[f'm_{variant}']
                    psi = (1 - w) * m + w * observed_loss
                    assert np.array_equal(psi, (1 - w) * m + w * changed_loss)
                    share = np.max(abs(psi), axis=0) / np.sum(abs(psi), axis=0)
                    for k, model in enumerate(models):
                        row = rows[rows.model.eq(model)].iloc[0]
                        max_difference = max(max_difference, abs(share[k] - row[f'share_{variant}']), abs(psi[:, k].mean() - row[variant]))
                counts += 1
        assert max_difference < 1e-12
        for (endpoint, anchor, mechanism, model, variant), records in d[d.phase.eq(phase)].groupby(KEYS + ['variant']):
            z = source[source.endpoint.eq(endpoint) & source.anchor.eq(anchor) & source.mechanism.eq(mechanism) & source.model.eq(model)]
            aux = z[z.method.eq('aipw') & z.propensity.eq('low') & z.loss_variant.eq(variant)].sort_values('seed')
            cc = z[z.method.eq('complete_case')].sort_values('seed')
            harmful = aux.abs_error.to_numpy() > cc.abs_error.to_numpy()
            scores = {'concentration': aux.top1_abs_share, 'max_weight': aux.max_normalized_weight,
                      'inverse_ess': -aux.ess_fraction, 'smd': aux.max_abs_smd, 'clipping': aux.clip_fraction}
            for _, record in records.iterrows():
                assert abs(pair_auc(harmful, scores[record.diagnostic].to_numpy()) - record.auc) < 1e-14
    assert counts == 1600
    assert protected_hashes == {p.name: digest(p) for p in protected.iterdir()}
    assert all(digest(ROOT / p) == h for p, h in manifest['input_hashes'].items())
    pd.DataFrame(contrasts).to_csv(OUT / 'verified_utility_contrasts.csv', index=False, float_format='%.17g')
    result = {'status': 'PASS_SPECIFIED_READBACK', 'choice_rows': len(c), 'summary_rows': len(s), 'diagnostic_rows': len(d),
              'nuisance_files_read': counts, 'max_expanded_formula_difference': max_difference,
              'unobserved_reference_perturbation_invariant': True, 'pairwise_auc_checks': len(d),
              'conditional_choice_bootstrap_recomputed': True, 'prior_inputs_unchanged': True,
              'verifier_sha256': digest(Path(__file__)),
              'scope': 'Independent arithmetic/readback code, not independent datasets or human scientific validation'}
    (OUT / 'verification.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
