import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
import warnings

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from revision_risk_diagnostics import matched_probability, metrics, probability_loss

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / 'artifacts/confirmatory/confirm_001/full_recovery_001'
BASE = ROOT / 'artifacts/revision_001/forensic'
OUT = ROOT / 'artifacts/revision_002'
MODELS = ['ecfp_rf', 'ecfp_desc_xgb']
CASE = 2945


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fhash(a):
    return hashlib.sha256(np.asarray(a, dtype='<f8').tobytes()).hexdigest()


def contexts():
    manifest = json.loads((BASE / 'refit_manifest.json').read_text())
    assert sha(BASE / 'local_feature_cache.npz') == manifest['output_sha256']['local_feature_cache.npz']
    cache = np.load(BASE / 'local_feature_cache.npz')
    names = json.loads((FROZEN / 'structure_audit.json').read_text())['descriptor_names']
    assert 'AvgIpc' in names
    ipc = 2048 + names.index('Ipc')
    replay = pd.read_parquet(FROZEN / 'row_replay_seed41001.parquet')
    pred = pd.read_parquet(FROZEN / 'predictions.parquet')
    gen = pd.read_parquet(FROZEN / 'generator_vectors.parquet')
    old = pd.read_csv(FROZEN / 'semisynthetic_scenarios.csv', float_precision='round_trip')
    result = {}
    for endpoint, part in replay.groupby('endpoint', sort=True):
        first = part.sort_values('source_row_index').drop_duplicates('source_row_index')
        ids = first.source_row_index.to_numpy()
        low = first[[x for x in first.columns if x.startswith('audit_')]].to_numpy(float)
        full = np.column_stack([cache['features'][np.searchsorted(cache['source_rows'], ids)], first.audit_max_train_similarity])
        y = np.column_stack([pred[pred.endpoint.eq(endpoint) & pred.model.eq(m)].sort_values('source_row_index').loss for m in MODELS])
        folds = first.fold.to_numpy()
        prepared = {}
        for variant, x in [('low', low), ('raw', full), ('drop', np.delete(full, ipc, axis=1))]:
            prepared[variant] = []
            for fold in range(3):
                train, valid = np.flatnonzero(folds != fold), np.flatnonzero(folds == fold)
                scaler = StandardScaler()
                xt = scaler.fit_transform(x[train])
                xv = scaler.transform(x[valid])
                xc = None
                if variant == 'raw':
                    clipped = x[valid].copy()
                    clipped[:, ipc] = np.clip(clipped[:, ipc], x[train, ipc].min(), x[train, ipc].max())
                    xc = scaler.transform(clipped)
                prepared[variant].append((train, valid, xt, xv, xc))
        gg = gen[gen.endpoint.eq(endpoint)]
        generators = {(a, m): g.sort_values('source_row_index').p_oracle.to_numpy() for (a, m), g in gg.groupby(['anchor', 'mechanism'])}
        original = old[old.endpoint.eq(endpoint)].set_index(['anchor', 'mechanism', 'seed', 'model', 'propensity_representation', 'loss_representation'])
        result[endpoint] = {'ids': ids, 'low': low, 'y': y, 'folds': folds, 'prepared': prepared,
                            'generators': generators, 'old': original, 'case': int(np.flatnonzero(ids == CASE)[0])}
    return result


def propensity(prepared, r):
    raw = np.empty(len(r))
    for train, valid, xt, xv, _ in prepared:
        model = LogisticRegression(C=1, solver='liblinear', max_iter=1000, random_state=20260910)
        model.fit(xt, r[train])
        assert model.n_iter_.max() < 1000
        raw[valid] = model.predict_proba(xv)[:, 1]
    return raw, np.clip(raw, .05, .95)


def loss_fit(prepared, r, y):
    prediction, clipped = np.empty_like(y), np.empty_like(y)
    for train, valid, xt, xv, xc in prepared:
        observed = r[train].astype(bool)
        model = Ridge(alpha=10, solver='lsqr', tol=1e-6).fit(xt[observed], y[train][observed])
        prediction[valid] = model.predict(xv)
        if xc is not None:
            clipped[valid] = model.predict(xc)
    return prediction, clipped


def weighting(r, p, low):
    w = r / p
    sd = low.std(axis=0)
    valid = sd > 0
    smd = np.abs((w @ low / w.sum() - low.mean(axis=0))[valid] / sd[valid])
    return {'ess_fraction': w.sum() ** 2 / (w @ w) / r.sum(), 'max_abs_smd': smd.max(),
            'max_normalized_weight': w.max() / w.sum()}


def job(context, endpoint, anchor, mechanism, seed, p, phase, destination):
    started = time.perf_counter()
    y, ids, case = context['y'], context['ids'], context['case']
    r = (np.random.default_rng(seed).random(len(y)) < p).astype(int)
    support = all(min(r[train].sum(), len(train) - r[train].sum()) >= 20 for train, *_ in context['prepared']['low'])
    assert support
    props = {'oracle': (p, p), 'low': propensity(context['prepared']['low'], r)}
    variants = ['low', 'drop']
    if phase == 'original':
        props['full'] = propensity(context['prepared']['raw'], r)
        variants += ['raw']
    fitted = {}
    for variant in variants:
        fitted[variant], clipped = loss_fit(context['prepared'][variant], r, y)
        if variant == 'raw':
            fitted['heldout_clip'] = clipped
    assert all(np.isfinite(x).all() for x in fitted.values())
    rows = []
    base = {'endpoint': endpoint, 'anchor': anchor, 'mechanism': mechanism, 'seed': seed, 'phase': phase,
            'n': len(y), 'n_observed': int(r.sum()), 'case_source_row': CASE, 'case_R': int(r[case]), 'fold_support': support}
    for k, model in enumerate(MODELS):
        loss, truth = y[:, k], y[:, k].mean()
        cc = loss[r.astype(bool)].mean()
        rows.append({**base, 'model': model, 'method': 'complete_case', 'propensity': 'none', 'loss_variant': 'none',
                     'estimate': cc, 'truth': truth, 'error': cc-truth, 'abs_error': abs(cc-truth)})
        for pr, (raw_p, use_p) in props.items():
            w = r / use_p
            diagnostic = {**weighting(r, use_p, context['low']), 'clip_fraction': float(((raw_p < .05) | (raw_p > .95)).mean()) if pr != 'oracle' else 0.,
                          'probability_rmse': np.sqrt(np.mean((raw_p-p)**2)), 'case_p': use_p[case]}
            for method, estimate in [('ht', np.mean(w*loss)), ('hajek', np.sum(w*loss)/w.sum())]:
                rows.append({**base, **diagnostic, 'model': model, 'method': method, 'propensity': pr, 'loss_variant': 'none',
                             'estimate': estimate, 'truth': truth, 'error': estimate-truth, 'abs_error': abs(estimate-truth)})
            for variant, predictions in fitted.items():
                m = predictions[:, k]
                residual = w * (loss-m)
                pseudo = m + residual
                estimate = pseudo.mean()
                scale = np.max(np.abs(pseudo))
                expanded = (1-w)*m+w*loss
                assert np.max(np.abs(pseudo-expanded)) <= 1e-10*max(1., scale)
                share = np.abs(pseudo)/np.abs(pseudo).sum()
                row = {**base, **diagnostic, 'model': model, 'method': 'aipw', 'propensity': pr, 'loss_variant': variant,
                       'estimate': estimate, 'truth': truth, 'error': estimate-truth, 'abs_error': abs(estimate-truth),
                       'plugin_mean': m.mean(), 'residual_mean': residual.mean(), 'm_min': m.min(), 'm_max': m.max(),
                       'm_negative_fraction': np.mean(m < 0), 'm_abs_gt10_fraction': np.mean(np.abs(m)>10),
                       'case_m': m[case], 'case_loss': loss[case], 'case_m_coefficient': 1-w[case],
                       'case_plugin_contribution': m[case]/len(y), 'case_residual_contribution': residual[case]/len(y),
                       'case_pseudo_contribution': pseudo[case]/len(y), 'case_abs_share': share[case],
                       'largest_source_row': int(ids[np.argmax(share)]), 'top1_abs_share': share.max(),
                       'top5_abs_share': np.sort(share)[-5:].sum(), 'expanded_identity_max_difference': np.max(np.abs(pseudo-expanded))}
                if phase == 'original' and variant in ['raw', 'low'] and pr in ['low', 'full']:
                    historical = context['old'].loc[(anchor, mechanism, seed, model, 'lowdim_11' if pr == 'low' else 'full_2266',
                                                     'lowdim_11' if variant == 'low' else 'full_2266')]
                    row.update(frozen_estimate=historical.fitted_aipw, refit_absolute_difference=abs(estimate-historical.fitted_aipw),
                               original_p_hash_exact=fhash(use_p)==historical.p_fitted_sha256,
                               original_m_hash_exact=fhash(m)==historical.mhat_sha256)
                rows.append(row)
    folder = destination / 'row_vectors' / endpoint / f'{anchor}_{mechanism}'
    folder.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(folder / f'{seed}.npz', R=r, **{f'p_{name}': q for name, (_, q) in props.items()},
                        **{f'p_raw_{name}': q for name, (q, _) in props.items()}, **{f'm_{name}': a for name, a in fitted.items()})
    return rows, {**base, 'seconds': time.perf_counter()-started}


def run(phase, smoke):
    destination = OUT / (phase + ('_smoke' if smoke else ''))
    destination.mkdir(exist_ok=False)
    initial = {p.name: sha(p) for p in FROZEN.iterdir()}
    source_hashes = {'protocol': sha(OUT/'protocol.md'), 'script': sha(Path(__file__)),
                     'revision001_zip': sha(ROOT/'paper/reproducibility/JPA_selective_reference_revision_20260912.zip')}
    start = time.time()
    ctx = contexts()
    tasks, designs = [], []
    for endpoint, context in ctx.items():
        for anchor in sorted({a for a, _ in context['generators']}):
            if phase == 'original':
                mechanisms = [(m, context['generators'][(anchor,m)]) for m in ['mcar','mar']]
                seeds = range(41001, 41002 if smoke else 41101)
            else:
                z = StandardScaler().fit_transform(context['low'][:, :1]).ravel()
                coverage = context['generators'][(anchor,'mar')].mean()
                mechanisms = [(f'observable_MolWt_{direction}', matched_probability(sign*2*z, coverage))
                              for direction, sign in [('positive',1),('negative',-1)]]
                seeds = range(71001, 71002 if smoke else 71101)
            for mechanism, p in mechanisms:
                for k, model in enumerate(MODELS):
                    designs.append({'endpoint':endpoint,'anchor':anchor,'mechanism':mechanism,'model':model,
                                    **probability_loss(context['y'][:,k],p)})
                tasks.extend((context,endpoint,anchor,mechanism,seed,p,phase,destination) for seed in seeds)
    rows, timings = [], []
    with threadpool_limits(limits=1), warnings.catch_warnings():
        warnings.simplefilter('error', ConvergenceWarning)
        warnings.simplefilter('ignore', FutureWarning)
        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = [pool.submit(job,*task) for task in tasks]
            for i, future in enumerate(as_completed(futures),1):
                value,timing = future.result()
                rows.extend(value)
                timings.append(timing)
                if i % 50 == 0 or i == len(tasks):
                    message = f'{phase}: {i}/{len(tasks)} settings; {time.time()-start:.1f}s'
                    print(message,flush=True)
                    with (destination/'progress.log').open('a',encoding='utf-8') as handle:
                        handle.write(message+'\n')
    frame = pd.DataFrame(rows).sort_values(['endpoint','anchor','mechanism','seed','model','method','propensity','loss_variant'])
    frame.to_csv(destination/'estimates.csv',index=False,float_format='%.17g')
    pd.DataFrame(timings).to_csv(destination/'timings.csv',index=False)
    pd.DataFrame(designs).to_csv(destination/'design.csv',index=False,float_format='%.17g')
    summaries = []
    for key, group in frame.groupby(['endpoint','anchor','mechanism','model','method','propensity','loss_variant']):
        ident = dict(zip(['endpoint','anchor','mechanism','model','method','propensity','loss_variant'],key))
        summaries.append({**ident, **metrics(group.estimate,group.truth.iloc[0])})
    pd.DataFrame(summaries).to_csv(destination/'summary.csv',index=False,float_format='%.17g')
    current = {p.name:sha(p) for p in FROZEN.iterdir()}
    assert initial == current and len(initial) == 19
    assert sha(ROOT/'paper/reproducibility/JPA_selective_reference_revision_20260912.zip') == source_hashes['revision001_zip']
    record = {'phase':phase,'smoke_only':smoke,'settings':len(tasks),'rows':len(frame),'duration_seconds':time.time()-start,
              'finished_utc':datetime.now(timezone.utc).isoformat(),'source_hashes':source_hashes,'frozen_files_unchanged':True,
              'frozen_sha256':initial,'output_hashes':{p.name:sha(p) for p in destination.glob('*.csv')}}
    if phase == 'original':
        compared = frame[frame.frozen_estimate.notna()]
        record.update(raw_grid_comparisons=len(compared),max_scalar_refit_difference=compared.refit_absolute_difference.max(),
                      materially_discrepant_count=int(((compared.refit_absolute_difference>1e-3)&
                          (compared.refit_absolute_difference>1e-3*np.abs(compared.frozen_estimate))).sum()),
                      p_vector_exact_count=int(compared.original_p_hash_exact.sum()),m_vector_exact_count=int(compared.original_m_hash_exact.sum()))
    (destination/'manifest.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k not in ['frozen_sha256','output_hashes','source_hashes']},indent=2),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--phase',choices=['original','bridge'],required=True)
    parser.add_argument('--smoke',action='store_true')
    args = parser.parse_args()
    run(args.phase,args.smoke)
