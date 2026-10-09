from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from scipy.linalg import solve
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / 'artifacts/confirmatory/confirm_001/full_recovery_001'
OUT = ROOT / 'artifacts/revision_001/forensic'
KEYS = ['endpoint', 'anchor', 'mechanism', 'model', 'propensity_representation', 'loss_representation']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def inspect():
    inputs = {p.name: sha(p) for p in FROZEN.iterdir() if p.is_file()}
    replay = pd.read_parquet(FROZEN / 'row_replay_seed41001.parquet')
    rows, influences = [], []
    for key, part in replay.groupby(KEYS, sort=True):
        meta = dict(zip(KEYS, key))
        m = part.mhat.to_numpy()
        row = {**meta, 'n': len(part), 'truth': part.loss.mean(), 'm_mean': m.mean(),
               'm_negative_fraction': (m < 0).mean(), 'm_abs_gt10_fraction': (abs(m) > 10).mean()}
        row.update({f'm_q{int(q * 100):02d}': v for q, v in zip([0, .01, .05, .5, .95, .99, 1], np.quantile(m, [0, .01, .05, .5, .95, .99, 1]))})
        for kind in ['fitted', 'oracle']:
            correction = part.R.to_numpy() * (part.loss.to_numpy() - m) / part[f'p_{kind}'].to_numpy()
            pseudo = m + correction
            assert np.array_equal(pseudo, part[f'pseudo_outcome_{kind}'].to_numpy())
            row.update({f'{kind}_correction_mean': correction.mean(), f'{kind}_aipw': pseudo.mean(),
                        f'{kind}_absolute_error': abs(pseudo.mean() - row['truth'])})
            centered = abs(pseudo - pseudo.mean())
            largest = np.argsort(centered)[-5:][::-1]
            for rank, i in enumerate(largest, 1):
                r = part.iloc[i]
                influences.append({**meta, 'probability': kind, 'rank': rank, 'source_row_index': int(r.source_row_index),
                                   'cohort_row_index': int(r.cohort_row_index), 'molecule_id': r.molecule_id,
                                   'fold': int(r.fold), 'R': int(r.R), 'loss': r.loss, 'mhat': r.mhat,
                                   'p': r[f'p_{kind}'], 'plugin_contribution': m[i] / len(part),
                                   'residual_contribution': correction[i] / len(part), 'pseudo_contribution': pseudo[i] / len(part),
                                   'absolute_pseudo_share': abs(pseudo[i]) / abs(pseudo).sum(),
                                   'absolute_centered_share': centered[i] / centered.sum()})
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / 'stored_decomposition.csv', index=False)
    pd.DataFrame(influences).to_csv(OUT / 'top_influence_rows.csv', index=False)

    cfg = json.loads((ROOT / 'artifacts/confirmatory/confirm_001/config.json').read_text())
    source = ROOT / cfg['source_freeze']['biogen']['analysis_derivative']['path']
    assert sha(source) == cfg['source_freeze']['biogen']['analysis_derivative']['sha256']
    data = pd.read_csv(source)
    structure = json.loads((FROZEN / 'structure_audit.json').read_text())
    names = structure['descriptor_names']
    funcs = dict(Descriptors._descList)
    fingerprint = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    unique = replay[['source_row_index']].drop_duplicates().sort_values('source_row_index').source_row_index.to_numpy()
    features = np.empty((len(unique), 2265), dtype=np.float64)
    for i, source_row in enumerate(unique):
        mol = Chem.MolFromSmiles(data.iloc[source_row].SMILES)
        DataStructs.ConvertToNumpyArray(fingerprint.GetFingerprint(mol), features[i, :2048])
        features[i, 2048:] = [funcs[name](mol) for name in names]
    assert np.isfinite(features).all()
    np.savez_compressed(OUT / 'local_feature_cache.npz', source_rows=unique, features=features)
    feature_names = [f'ECFP_{i}' for i in range(2048)] + names + ['max_train_similarity']
    feature_rows = []
    for endpoint, part in replay.drop_duplicates(['endpoint', 'source_row_index']).groupby('endpoint'):
        part = part.sort_values('source_row_index')
        X = np.column_stack([features[np.searchsorted(unique, part.source_row_index)], part.audit_max_train_similarity])
        for fold in sorted(part.fold.unique()):
            train, valid = part.fold.to_numpy() != fold, part.fold.to_numpy() == fold
            scaler = StandardScaler().fit(X[train])
            zvalid = scaler.transform(X[valid])
            for j in np.argsort(abs(zvalid).max(axis=0))[-20:][::-1]:
                local_i = np.argmax(abs(zvalid[:, j]))
                feature_rows.append({'endpoint': endpoint, 'fold': int(fold), 'feature': feature_names[j],
                                     'train_min': X[train, j].min(), 'train_max': X[train, j].max(),
                                     'train_mean': scaler.mean_[j], 'train_sd': np.sqrt(scaler.var_[j]),
                                     'train_scale': scaler.scale_[j], 'train_zero_variance': bool(scaler.var_[j] == 0),
                                     'valid_min': X[valid, j].min(), 'valid_max': X[valid, j].max(),
                                     'max_abs_valid_z': abs(zvalid[:, j]).max(),
                                     'max_z_source_row': int(part.loc[valid].iloc[local_i].source_row_index)})
    pd.DataFrame(feature_rows).to_csv(OUT / 'fold_feature_extremes.csv', index=False)
    save_json('inspection_manifest.json', {'role': 'post_result_diagnostic', 'seed': 41001, 'rdkit': rdBase.rdkitVersion,
                                         'frozen_input_sha256': inputs, 'source_sha256': sha(source),
                                         'script_sha256': sha(Path(__file__)), 'n_replay_rows': len(replay),
                                         'n_cells': len(rows), 'n_unique_source_rows': len(unique)})
    print(pd.DataFrame(feature_rows).sort_values('max_abs_valid_z', ascending=False).head(10).to_string(index=False))


def refit():
    replay = pd.read_parquet(FROZEN / 'row_replay_seed41001.parquet')
    cache = np.load(OUT / 'local_feature_cache.npz')
    names = json.loads((FROZEN / 'structure_audit.json').read_text())['descriptor_names']
    ipc = 2048 + names.index('Ipc')
    base = replay.query("propensity_representation == 'lowdim_11' and loss_representation == 'full_2266'")
    model_names = ['ecfp_rf', 'ecfp_desc_xgb']
    rows, fits, coefficients, estimates = [], [], [], []
    for key, part in base.groupby(['endpoint', 'anchor', 'mechanism'], sort=True):
        endpoint, anchor, mechanism = key
        meta = {'endpoint': endpoint, 'anchor': anchor, 'mechanism': mechanism}
        first = part[part.model.eq(model_names[0])].sort_values('source_row_index')
        ids = first.source_row_index.to_numpy()
        X = np.column_stack([cache['features'][np.searchsorted(cache['source_rows'], ids)], first.audit_max_train_similarity])
        Y = np.column_stack([part[part.model.eq(model)].sort_values('source_row_index').loss for model in model_names])
        original = np.column_stack([part[part.model.eq(model)].sort_values('source_row_index').mhat for model in model_names])
        r = first.R.to_numpy().astype(bool)
        fold_ids = first.fold.to_numpy()
        predictions = {name: np.zeros_like(Y) for name in ['lsqr_replay', 'independent_dual', 'ipc_train_range_clip']}
        for fold in range(3):
            train, valid = fold_ids != fold, fold_ids == fold
            scaler = StandardScaler().fit(X[train])
            xt = scaler.transform(X[train])[r[train]]
            xv = scaler.transform(X[valid])
            yt = Y[train][r[train]]
            ridge = Ridge(alpha=10, solver='lsqr', tol=1e-6).fit(xt, yt)
            predictions['lsqr_replay'][valid] = ridge.predict(xv)
            xc, yc = xt - xt.mean(axis=0), yt - yt.mean(axis=0)
            dual = solve(xc @ xc.T + 10 * np.eye(len(xt)), yc, assume_a='pos')
            coef = xc.T @ dual
            normal_residual = xc.T @ (xc @ coef - yc) + 10 * coef
            predictions['independent_dual'][valid] = (xv - xt.mean(axis=0)) @ coef + yt.mean(axis=0)
            clipped = X[valid].copy()
            clipped[:, ipc] = np.clip(clipped[:, ipc], X[train, ipc].min(), X[train, ipc].max())
            predictions['ipc_train_range_clip'][valid] = ridge.predict(scaler.transform(clipped))
            fits.append({**meta, 'fold': fold, 'n_observed_train': len(xt),
                         'dual_normal_residual_max': abs(normal_residual).max(),
                         'dual_relative_normal_residual': np.linalg.norm(normal_residual) / np.linalg.norm(xc.T @ yc),
                         'max_coef_difference_lsqr_dual': abs(ridge.coef_.T - coef).max(),
                         'max_prediction_difference_lsqr_dual': abs(predictions['lsqr_replay'][valid] - predictions['independent_dual'][valid]).max(),
                         'max_prediction_replay_abs_difference': abs(predictions['lsqr_replay'][valid] - original[valid]).max(),
                         'ipc_validation_clipped_count': int(np.count_nonzero(clipped[:, ipc] != X[valid, ipc]))})
            for model_i, model in enumerate(model_names):
                local = np.argmax(abs(predictions['lsqr_replay'][valid, model_i]))
                coefficients.append({**meta, 'fold': fold, 'model': model,
                                     'source_row_index': int(ids[valid][local]), 'Ipc': X[valid, ipc][local],
                                     'Ipc_scaled': xv[local, ipc], 'Ipc_ridge_coefficient': ridge.coef_[model_i, ipc],
                                     'Ipc_prediction_contribution': xv[local, ipc] * ridge.coef_[model_i, ipc],
                                     'whole_prediction': predictions['lsqr_replay'][valid, model_i][local],
                                     'independent_Ipc_coefficient': coef[ipc, model_i],
                                     'clipped_prediction': predictions['ipc_train_range_clip'][valid, model_i][local]})
        for model_i, model in enumerate(model_names):
            for i, source_row in enumerate(ids):
                rows.append({**meta, 'model': model, 'source_row_index': int(source_row), 'fold': int(fold_ids[i]),
                             'R': int(r[i]), 'loss': Y[i, model_i], 'stored_mhat': original[i, model_i],
                             **{name: value[i, model_i] for name, value in predictions.items()}})
            for propensity in ['lowdim_11', 'full_2266']:
                pr = replay[(replay.endpoint.eq(endpoint)) & replay.anchor.eq(anchor) & replay.mechanism.eq(mechanism)
                            & replay.model.eq(model) & replay.propensity_representation.eq(propensity)
                            & replay.loss_representation.eq('full_2266')].sort_values('source_row_index')
                for probability in ['fitted', 'oracle']:
                    p = pr[f'p_{probability}'].to_numpy()
                    for variant, pred in {'stored': original, **predictions}.items():
                        m = pred[:, model_i]
                        correction = r * (Y[:, model_i] - m) / p
                        estimate = np.mean(m + correction)
                        estimates.append({**meta, 'model': model, 'propensity_representation': propensity, 'probability': probability,
                                          'variant': variant, 'truth': Y[:, model_i].mean(), 'm_mean': m.mean(),
                                          'residual_mean': correction.mean(), 'aipw': estimate,
                                          'absolute_error': abs(estimate - Y[:, model_i].mean()),
                                          'm_min': m.min(), 'm_max': m.max(), 'm_negative_fraction': (m < 0).mean()})
    pd.DataFrame(rows).to_csv(OUT / 'numerical_refit_predictions.csv', index=False)
    pd.DataFrame(fits).to_csv(OUT / 'numerical_solver_checks.csv', index=False)
    pd.DataFrame(coefficients).to_csv(OUT / 'ipc_coefficient_contributions.csv', index=False)
    pd.DataFrame(estimates).to_csv(OUT / 'stable_sensitivity_estimates.csv', index=False)
    previous = json.loads((OUT / 'inspection_manifest.json').read_text())['frozen_input_sha256']
    current = {p.name: sha(p) for p in FROZEN.iterdir() if p.is_file()}
    assert previous == current
    output_hashes = {p.name: sha(p) for p in OUT.iterdir() if p.is_file() and p.name != 'refit_manifest.json'}
    save_json('refit_manifest.json', {'role': 'post_result_diagnostic', 'seed': 41001, 'cells': 8, 'fold_fits': len(fits),
                                    'frozen_files_unchanged': previous == current, 'frozen_file_count': len(current),
                                    'script_sha256': sha(Path(__file__)), 'output_sha256': output_hashes,
                                    'max_replay_prediction_difference': max(x['max_prediction_replay_abs_difference'] for x in fits),
                                    'max_dual_relative_normal_residual': max(x['dual_relative_normal_residual'] for x in fits)})
    print(pd.DataFrame(fits).to_string(index=False))
    print(pd.DataFrame(estimates).groupby(['probability', 'variant']).absolute_error.agg(['min', 'median', 'max']).to_string())


if __name__ == '__main__':
    with threadpool_limits(limits=2):
        {'inspect': inspect, 'refit': refit}[sys.argv[1]]()
