import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import rankdata


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'artifacts/confirmatory/confirm_001/full_recovery_001'
OUT = ROOT / 'artifacts/revision_001/risk'
KEYS = ['endpoint', 'anchor', 'model', 'mechanism']
REPS = ['propensity_representation', 'loss_representation']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(frame, name):
    frame.to_csv(OUT / name, index=False, float_format='%.17g')


def metrics(values, truth, se=None):
    values = np.asarray(values)
    error = values - truth
    sd = values.std(ddof=1)
    result = {'replicates': len(values), 'truth': truth, 'mean_estimate': values.mean(),
              'signed_bias': error.mean(), 'mc_sd': sd, 'variance': sd ** 2,
              'rmse': np.sqrt(np.mean(error ** 2)), 'median_abs_error': np.median(np.abs(error)),
              'mean_abs_error': np.mean(np.abs(error)), 'bias_mcse': sd / np.sqrt(len(values))}
    if se is not None:
        se = np.broadcast_to(se, values.shape)
        covered = np.abs(error) <= 1.96 * se + 1e-14
        proportion = covered.mean()
        denom = 1 + 1.96 ** 2 / len(values)
        center = (proportion + 1.96 ** 2 / (2 * len(values))) / denom
        half = 1.96 * np.sqrt(proportion * (1 - proportion) / len(values) + 1.96 ** 2 / (4 * len(values) ** 2)) / denom
        result.update(mean_se=se.mean(), mean_interval_width=3.92 * se.mean(),
                      median_interval_width=np.median(3.92 * se), coverage=covered.mean(),
                      coverage_mcse=np.sqrt(covered.mean() * (1 - covered.mean()) / len(values)),
                      coverage_wilson95_lower=center - half, coverage_wilson95_upper=center + half,
                      mean_se_over_mc_sd=se.mean() / sd if sd > 0 else np.nan)
    return result


def probability_loss(loss, p):
    truth = loss.mean()
    cov = np.mean((p - p.mean()) * (loss - truth))
    shift = np.average(loss, weights=p) - truth
    assert abs(shift - cov / p.mean()) < 1e-12
    corr = cov / (p.std() * loss.std()) if p.std() > 1e-12 else 0.0
    return {'n_total': len(loss), 'truth': truth, 'mean_probability': p.mean(),
            'p_min': p.min(), 'p_max': p.max(), 'p_loss_covariance_population': cov,
            'p_loss_correlation': corr, 'probability_weighted_risk': truth + shift,
            'probability_weighted_shift': shift, 'shift_relative_to_truth': shift / truth}


def stored_diagnostics(scenarios, pred, gen, replay):
    summaries = []
    for key, group in scenarios.groupby(KEYS + REPS):
        base = dict(zip(KEYS + REPS, key))
        truth = group.full_truth.iloc[0]
        for method in ['cc', 'or_mean', 'fitted_ipw', 'fitted_aipw', 'oracle_ipw', 'oracle_aipw']:
            se = group[method.replace('aipw', 'se')] if method.endswith('aipw') else None
            summaries.append({**base, 'method': method, **metrics(group[method], truth, se)})
    save(pd.DataFrame(summaries), 'original_estimator_summary.csv')
    associations = []
    for (endpoint, model), pp in pred.groupby(['endpoint', 'model']):
        pp = pp.sort_values('source_row_index')
        for (anchor, mechanism), gg in gen[gen.endpoint.eq(endpoint)].groupby(['anchor', 'mechanism']):
            gg = gg.sort_values('source_row_index')
            assert np.array_equal(pp.source_row_index, gg.source_row_index)
            associations.append({'endpoint': endpoint, 'anchor': anchor, 'model': model,
                                 'mechanism': mechanism, **probability_loss(pp.loss.to_numpy(), gg.p_oracle.to_numpy())})
    associations = pd.DataFrame(associations)
    save(associations, 'original_probability_loss_association.csv')
    candidate = scenarios[scenarios.propensity_representation.eq('lowdim_11') & scenarios.loss_representation.eq('lowdim_11')]
    gates = []
    for key, group in candidate.groupby(KEYS):
        base = dict(zip(KEYS, key))
        threshold = 1.25 if key[-1] == 'mcar' else .75
        cc = group.cc_abs_error.median()
        for prefix in ['fitted', 'oracle']:
            ratio = group[f'{prefix}_aipw_abs_error'].median() / cc
            clip = group.fitted_clip_fraction.to_numpy() if prefix == 'fitted' else np.zeros(len(group))
            checks = {
                'error_ratio_pass': ratio <= threshold,
                'coverage_pass': group[f'{prefix}_covered'].mean() >= .9,
                'ess_pass': group[f'{prefix}_ess_fraction'].quantile(.05) >= .5,
                'median_clipping_pass': np.median(clip) <= .1,
                'q95_clipping_pass': np.quantile(clip, .95) <= .2,
                'smd_pass': group[f'{prefix}_max_abs_smd'].quantile(.95) <= .1,
                'fold_support_pass': group.fold_class_support.all(),
                'finite_pass': group.all_finite.all(),
            }
            gates.append({**base, 'method': f'{prefix}_aipw_lowM', 'original_candidate': prefix == 'fitted',
                          'median_cc_abs_error': cc, 'median_abs_error': group[f'{prefix}_aipw_abs_error'].median(),
                          'error_ratio': ratio, 'error_ratio_threshold': threshold,
                          'coverage': group[f'{prefix}_covered'].mean(),
                          'ess_fraction_q05': group[f'{prefix}_ess_fraction'].quantile(.05),
                          'max_abs_smd_q95': group[f'{prefix}_max_abs_smd'].quantile(.95),
                          **checks, 'all_support_diagnostics_pass': all(v for k, v in checks.items() if k not in ['error_ratio_pass', 'coverage_pass']),
                          'all_original_gate_components_pass': all(checks.values())})
    gates = pd.DataFrame(gates).merge(associations, on=KEYS, validate='many_to_one')
    cc_summary = pd.DataFrame(summaries)
    cc_summary = cc_summary[cc_summary.method.eq('cc') & cc_summary.propensity_representation.eq('lowdim_11') & cc_summary.loss_representation.eq('lowdim_11')]
    gates = gates.merge(cc_summary[KEYS + ['signed_bias', 'mc_sd', 'rmse']].rename(columns={'signed_bias': 'cc_signed_bias', 'mc_sd': 'cc_mc_sd', 'rmse': 'cc_rmse'}), on=KEYS, validate='many_to_one')
    save(gates, 'original_same_cell_gate_chain.csv')
    save(gates[gates.original_candidate], 'original_candidate_16_cells.csv')
    subset_rows = []
    for (endpoint, model), pp in pred.groupby(['endpoint', 'model']):
        pp = pp.sort_values('source_row_index')
        keep = pp.reference.to_numpy() > pp.reference.min()
        for (anchor, mechanism), gg in gen[gen.endpoint.eq(endpoint)].groupby(['anchor', 'mechanism']):
            gg = gg.sort_values('source_row_index')
            p = gg.p_oracle.to_numpy()
            masks = np.array([np.random.default_rng(seed).random(len(p)) < p for seed in range(41001, 41101)])
            for label, select in [('all_measured', np.ones(len(p), dtype=bool)), ('strictly_above_released_minimum', keep)]:
                loss, pr, r = pp.loss.to_numpy()[select], p[select], masks[:, select]
                truth = loss.mean()
                cc = (r @ loss) / r.sum(axis=1)
                w = r / pr
                ipw = (w @ loss) / w.sum(axis=1)
                for method, vals in [('cc', cc), ('oracle_hajek', ipw)]:
                    subset_rows.append({'endpoint': endpoint, 'anchor': anchor, 'model': model, 'mechanism': mechanism,
                                        'subset': label, 'method': method, 'released_minimum': pp.reference.min(),
                                        'n_at_minimum_full_cohort': (~keep).sum(), **probability_loss(loss, pr), **metrics(vals, truth)})
                if label == 'all_measured':
                    original = candidate[candidate.endpoint.eq(endpoint) & candidate.model.eq(model) & candidate.anchor.eq(anchor) & candidate.mechanism.eq(mechanism)].sort_values('seed')
                    assert np.max(np.abs(cc - original.cc.to_numpy())) < 1e-12
                    assert np.max(np.abs(ipw - original.oracle_ipw.to_numpy())) < 1e-12
    save(pd.DataFrame(subset_rows), 'nonminimum_sensitivity.csv')
    return gates


def matched_probability(score, coverage):
    lo, hi = -50., 50.
    for _ in range(100):
        mid = (lo + hi) / 2
        if np.clip(expit(score + mid), .1, .9).mean() < coverage:
            lo = mid
        else:
            hi = mid
    probability = np.clip(expit(score + (lo + hi) / 2), .1, .9)
    assert abs(probability.mean() - coverage) < 1e-12
    return probability


def calibration(pred, gen, replay):
    from scipy.sparse import csr_matrix
    n_replicates = 10000
    summaries, mechanisms, p1_rows, raw = [], [], [], []
    low_replay = replay[replay.propensity_representation.eq('lowdim_11') & replay.loss_representation.eq('lowdim_11') & replay.mechanism.eq('mar')]
    for endpoint_i, (endpoint, ep) in enumerate(pred.groupby('endpoint')):
        n = ep.cohort_row_index.nunique()
        uniforms = np.random.default_rng(61001 + endpoint_i).random((n_replicates, n))
        groups = low_replay[low_replay.endpoint.eq(endpoint)].sort_values('source_row_index').drop_duplicates('source_row_index').scaffold_group.to_numpy()
        code = pd.factorize(groups, sort=False)[0]
        g = code.max() + 1
        incidence = csr_matrix((np.ones(n), (np.arange(n), code)), shape=(n, g))
        group_n = np.bincount(code)

        def cluster_se(pseudo):
            estimates = pseudo.mean(axis=1)
            totals = (incidence.T @ pseudo.T).T - estimates[:, None] * group_n
            return np.sqrt(g / (g - 1) * np.square(totals).sum(axis=1)) / n

        for model, pp in ep.groupby('model'):
            pp = pp.sort_values('source_row_index')
            loss = pp.loss.to_numpy()
            truth = loss.mean()
            ranks = rankdata(loss, method='average')
            score = (ranks - ranks.mean()) / ranks.std()
            p1_se = cluster_se(loss[None, :])[0]
            p1_rows.append({'endpoint': endpoint, 'model': model, 'n_total': n, 'truth': truth,
                            'p': 1., 'true_mask_variance': 0., 'true_mask_se': 0.,
                            'legacy_centered_pseudo_cluster_se': p1_se,
                            'legacy_interval_width': 3.92 * p1_se})
            for anchor, gg in gen[gen.endpoint.eq(endpoint) & gen.mechanism.eq('mar')].groupby('anchor'):
                gg = gg.sort_values('source_row_index')
                coverage = gg.p_oracle.mean()
                fixed_m = low_replay[low_replay.endpoint.eq(endpoint) & low_replay.model.eq(model) & low_replay.anchor.eq(anchor)].sort_values('source_row_index').mhat.to_numpy()
                assert len(fixed_m) == n
                probability_specs = [('original_mcar', np.full(n, coverage)),
                                     ('original_mar', gg.p_oracle.to_numpy()),
                                     ('loss_rank_beta_0.5', matched_probability(.5 * score, coverage)),
                                     ('loss_rank_beta_2.0', matched_probability(2.0 * score, coverage))]
                for mechanism, p in probability_specs:
                    ident = {'endpoint': endpoint, 'anchor': anchor, 'model': model, 'mechanism': mechanism}
                    mechanisms.append({**ident, **probability_loss(loss, p), 'rng_seed': 61001 + endpoint_i,
                                       'replicates': n_replicates, 'probability_source': 'frozen_original' if mechanism.startswith('original_') else 'post_result_loss_directed_not_observational_MAR'})
                    vals = {name: [] for name in ['cc', 'oracle_ht', 'oracle_hajek', 'true_p_fixed_mean_aipw', 'constant_p_fixed_mean_aipw', 'true_p_fixed_replay_m_aipw']}
                    ses = {name: [] for name in ['true_p_fixed_mean_aipw', 'true_p_fixed_replay_m_aipw']}
                    residual = loss - truth
                    se_design = np.sqrt(np.sum((1 - p) / p * residual ** 2)) / n
                    se_ht = np.sqrt(np.sum((1 - p) / p * loss ** 2)) / n
                    se_replay = np.sqrt(np.sum((1 - p) / p * (loss - fixed_m) ** 2)) / n
                    se_wrong_actual = np.sqrt(np.sum(p * (1 - p) / coverage ** 2 * residual ** 2)) / n
                    se_wrong_claimed = np.sqrt(np.sum((1 - coverage) / coverage * residual ** 2)) / n
                    for start in range(0, n_replicates, 500):
                        r = uniforms[start:start + 500] < p
                        w = r / p
                        cc = (r @ loss) / r.sum(axis=1)
                        ht = (w @ loss) / n
                        hajek = ht * n / w.sum(axis=1)
                        pseudo = truth + w * residual
                        learned_pseudo = fixed_m + w * (loss - fixed_m)
                        computed = {'cc': cc, 'oracle_ht': ht, 'oracle_hajek': hajek,
                                    'true_p_fixed_mean_aipw': pseudo.mean(axis=1),
                                    'constant_p_fixed_mean_aipw': truth + (r @ residual) / (n * coverage),
                                    'true_p_fixed_replay_m_aipw': learned_pseudo.mean(axis=1)}
                        for method, value in computed.items():
                            vals[method].append(value)
                        ses['true_p_fixed_mean_aipw'].append(cluster_se(pseudo))
                        ses['true_p_fixed_replay_m_aipw'].append(cluster_se(learned_pseudo))
                    vals = {k: np.concatenate(v) for k, v in vals.items()}
                    ses = {k: np.concatenate(v) for k, v in ses.items()}
                    cc_med = np.median(np.abs(vals['cc'] - truth))
                    base_raw = {**ident, 'replicate': np.arange(n_replicates), **vals,
                                'legacy_cluster_se_fixed_mean': ses['true_p_fixed_mean_aipw'],
                                'legacy_cluster_se_fixed_replay_m': ses['true_p_fixed_replay_m_aipw']}
                    raw.append(pd.DataFrame(base_raw))
                    for method, value in vals.items():
                        exact_se = {'oracle_ht': se_ht, 'true_p_fixed_mean_aipw': se_design,
                                    'true_p_fixed_replay_m_aipw': se_replay,
                                    'constant_p_fixed_mean_aipw': se_wrong_actual}.get(method)
                        stat = metrics(value, truth, exact_se)
                        ratio = stat['median_abs_error'] / cc_med
                        theoretical_bias = np.sum(p * residual) / (n * coverage) if method == 'constant_p_fixed_mean_aipw' else (0. if method in ['oracle_ht', 'true_p_fixed_mean_aipw', 'true_p_fixed_replay_m_aipw'] else np.nan)
                        threshold = 1.25 if mechanism == 'original_mcar' else .75
                        row = {**ident, 'method': method, **stat, 'median_cc_abs_error': cc_med,
                               'median_abs_error_ratio_vs_cc': ratio, 'original_error_component_threshold': threshold,
                               'original_error_component_pass': ratio <= threshold,
                               'full_original_gate_evaluated': False, 'analytical_mask_se': exact_se,
                               'analytical_signed_bias': theoretical_bias,
                               'normal_ci_uses_true_design_variance': exact_se is not None,
                               'original_gate_error_component_only': True}
                        if method in ses:
                            legacy = metrics(value, truth, ses[method])
                            row.update({f'legacy_cluster_{k}': legacy[k] for k in ['mean_se', 'mean_interval_width', 'coverage', 'mean_se_over_mc_sd']})
                        if method == 'constant_p_fixed_mean_aipw':
                            claimed = metrics(value, truth, se_wrong_claimed)
                            row.update(misspecified_assumption_se=se_wrong_claimed,
                                       misspecified_assumption_ci_coverage=claimed['coverage'])
                        summaries.append(row)
                    print(f'Completed {endpoint} {model} {anchor} {mechanism}', flush=True)
    save(pd.DataFrame(summaries), 'mechanism_calibration_summary.csv')
    save(pd.DataFrame(mechanisms), 'mechanism_calibration_design.csv')
    save(pd.DataFrame(p1_rows), 'complete_observation_variance_check.csv')
    pd.concat(raw, ignore_index=True).to_parquet(OUT / 'mechanism_calibration_replicates.parquet', index=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((SOURCE / 'run_manifest.json').read_text())
    before = {p.name: sha(p) for p in SOURCE.iterdir() if p.is_file()}
    assert len(before) == 19 and not any(p.is_dir() for p in SOURCE.iterdir())
    for name, digest in manifest['output_hashes_before_manifest'].items():
        assert before[name] == digest
    inputs = ['semisynthetic_scenarios.csv', 'predictions.parquet', 'generator_vectors.parquet', 'row_replay_seed41001.parquet']
    scenarios = pd.read_csv(SOURCE / inputs[0], float_precision='round_trip')
    pred, gen, replay = [pd.read_parquet(SOURCE / name) for name in inputs[1:]]
    gates = stored_diagnostics(scenarios, pred, gen, replay)
    calibration(pred, gen, replay)
    after = {p.name: sha(p) for p in SOURCE.iterdir() if p.is_file()}
    assert before == after
    summary = {'analysis_id': 'revision_001_risk', 'status': 'POST_RESULT_DIAGNOSTIC_COMPLETE',
               'scope': 'Original outputs unchanged. New controlled masks calibrate descriptive benefit and interval behavior, not original gate validation.',
               'input_hashes': before, 'source_code_sha256': sha(Path(__file__)),
               'original_19_file_directory_unchanged': before == after,
               'original_mar_fitted_benefit_pass': int(gates.loc[gates.mechanism.eq('mar') & gates.method.eq('fitted_aipw_lowM'), 'error_ratio_pass'].sum()),
               'original_mar_oracle_benefit_pass': int(gates.loc[gates.mechanism.eq('mar') & gates.method.eq('oracle_aipw_lowM'), 'error_ratio_pass'].sum()),
               'original_mar_all_support_pass': int(gates.loc[gates.mechanism.eq('mar') & gates.method.eq('fitted_aipw_lowM'), 'all_support_diagnostics_pass'].sum()),
               'outputs': {p.name: sha(p) for p in OUT.iterdir() if p.is_file() and p.name != 'summary.json'}}
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in summary.items() if k not in ['input_hashes', 'outputs']}, indent=2))


if __name__ == '__main__':
    main()
