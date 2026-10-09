import warnings

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.special import expit, logit
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler


warnings.filterwarnings(
    "ignore",
    message="'penalty' was deprecated",
    category=FutureWarning,
    module="sklearn.linear_model._logistic",
)


def scaffold_groups(smiles):
    groups = []
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            raise ValueError(f"invalid SMILES at row {i}")
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(
            mol=mol, includeChirality=True
        )
        groups.append(scaffold or f"ACYCLIC_SINGLETON_{i}")
    return np.asarray(groups, dtype=object)


def _folds_and_scaled(X, groups, folds):
    split = list(GroupKFold(n_splits=folds).split(X, groups=groups))
    prepared = []
    sparse = hasattr(X, "tocsr")
    for train, valid in split:
        scaler = StandardScaler(with_mean=not sparse)
        x_train = scaler.fit_transform(X[train])
        x_valid = scaler.transform(X[valid])
        prepared.append((train, valid, x_train, x_valid))
    return prepared


def prepare_folds(X, groups, folds):
    return _folds_and_scaled(X, groups, folds)


def _logistic(cfg):
    return LogisticRegression(
        penalty=cfg["penalty"],
        C=float(cfg["C"]),
        solver=cfg["solver"],
        class_weight=cfg.get("class_weight"),
        max_iter=int(cfg["max_iter"]),
        random_state=int(cfg["random_state"]),
    )


def _weights(r, p):
    w = np.asarray(r, dtype=float) / np.asarray(p, dtype=float)
    ess = float(w.sum() ** 2 / np.square(w).sum())
    observed = int(np.asarray(r, dtype=bool).sum())
    return w, ess, ess / observed


def crossfit_propensity(X, groups, r, folds, prop_cfg, clip, prepared=None):
    r = np.asarray(r, dtype=int)
    if len(np.unique(r)) != 2:
        raise ValueError("observation indicator must contain both classes")
    p_raw = np.empty(len(r), dtype=float)
    prepared = prepared if prepared is not None else _folds_and_scaled(X, groups, folds)
    for train, valid, x_train, x_valid in prepared:
        if len(np.unique(r[train])) != 2:
            raise ValueError("a propensity training fold contains one class")
        model = _logistic(prop_cfg).fit(x_train, r[train])
        p_raw[valid] = model.predict_proba(x_valid)[:, 1]
    lo, hi = map(float, clip)
    p = np.clip(p_raw, lo, hi)
    _, ess, ess_fraction = _weights(r, p)
    return {
        "raw": p_raw,
        "clipped": p,
        "auc": float(roc_auc_score(r, p_raw)),
        "clip_fraction": float(np.mean((p_raw < lo) | (p_raw > hi))),
        "ess": ess,
        "ess_fraction": ess_fraction,
    }


def _cluster_se(phi, groups):
    phi = np.asarray(phi, dtype=float)
    codes = pd.factorize(np.asarray(groups, dtype=object), sort=False)[0]
    sums = np.bincount(codes, weights=phi)
    g = len(sums)
    if g < 2:
        raise ValueError("cluster influence-function SE needs at least two groups")
    return float(np.sqrt((g / (g - 1)) * np.square(sums).sum()) / len(phi))


def crossfit_estimators(X, groups, r, losses, folds, prop_cfg, ridge_cfg, clip, prepared=None):
    r = np.asarray(r, dtype=int)
    losses = np.asarray(losses, dtype=float)
    if losses.ndim == 1:
        losses = losses[:, None]
    if losses.shape[0] != len(r) or not np.isfinite(losses).all():
        raise ValueError("losses must be finite with shape (n, k)")
    p_raw = np.empty(len(r), dtype=float)
    mhat = np.empty_like(losses, dtype=float)
    prepared = prepared if prepared is not None else _folds_and_scaled(X, groups, folds)
    for train, valid, x_train, x_valid in prepared:
        if len(np.unique(r[train])) != 2:
            raise ValueError("a nuisance training fold contains one observation class")
        p_raw[valid] = _logistic(prop_cfg).fit(
            x_train, r[train]
        ).predict_proba(x_valid)[:, 1]
        observed = r[train].astype(bool)
        if observed.sum() == 0:
            raise ValueError("a nuisance training fold has no observed losses")
        ridge = Ridge(
            alpha=float(ridge_cfg["alpha"]),
            solver=ridge_cfg["solver"],
            tol=float(ridge_cfg["tol"]),
        ).fit(x_train[observed], losses[train][observed])
        pred = np.asarray(ridge.predict(x_valid), dtype=float)
        mhat[valid] = pred[:, None] if pred.ndim == 1 else pred
    lo, hi = map(float, clip)
    p = np.clip(p_raw, lo, hi)
    w, ess, ess_fraction = _weights(r, p)
    clip_fraction = float(np.mean((p_raw < lo) | (p_raw > hi)))
    estimates = []
    observed = r.astype(bool)
    for j in range(losses.shape[1]):
        loss = losses[:, j]
        cc = float(loss[observed].mean())
        ipw = float(np.sum(w * loss) / np.sum(w))
        pseudo = mhat[:, j] + r * (loss - mhat[:, j]) / p
        aipw = float(pseudo.mean())
        phi = pseudo - aipw
        se = _cluster_se(phi, groups)
        estimates.append(
            {
                "column": j,
                "cc": cc,
                "ipw": ipw,
                "aipw": aipw,
                "se": se,
                "lcl": aipw - 1.96 * se,
                "ucl": aipw + 1.96 * se,
                "ess": ess,
                "ess_fraction": ess_fraction,
                "clip_fraction": clip_fraction,
            }
        )
    return {"p_raw": p_raw, "p": p, "mhat": mhat, "estimates": estimates}


def bernoulli_mask(probability, seed):
    p = np.asarray(probability, dtype=float)
    if np.any((p < 0) | (p > 1) | ~np.isfinite(p)):
        raise ValueError("Bernoulli probabilities must be finite in [0, 1]")
    return (np.random.default_rng(seed).random(p.shape) < p).astype(int)


def mcar_mask(n, coverage, seed):
    return bernoulli_mask(np.full(int(n), float(coverage)), seed)


def mar_mask(base_p, seed):
    return bernoulli_mask(base_p, seed)


def mnar_probability(base_p, y, gamma):
    base_p = np.asarray(base_p, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(base_p) != len(y) or np.any((base_p <= 0) | (base_p >= 1)):
        raise ValueError("base_p must be strictly inside (0, 1) and align with y")
    sd = float(y.std(ddof=0))
    if not np.isfinite(sd) or sd == 0:
        raise ValueError("MNAR outcome must have nonzero finite population SD")
    z = (y - y.mean()) / sd
    target = float(base_p.mean())
    base_logit = logit(base_p) + float(gamma) * z
    low, high = -50.0, 50.0
    for _ in range(100):
        mid = (low + high) / 2
        if expit(base_logit + mid).mean() < target:
            low = mid
        else:
            high = mid
    shift = (low + high) / 2
    return expit(base_logit + shift), float(shift)


def mnar_mask(base_p, y, gamma, seed):
    probability, shift = mnar_probability(base_p, y, gamma)
    return bernoulli_mask(probability, seed), probability, shift


def gate4(diag, cfg):
    table = pd.DataFrame(diag).copy()
    required = ["mask", "auc", "ess_fraction"]
    expected = set(cfg["mask_endpoints"])
    if not all(c in table for c in required) or len(table) != len(expected) or set(table["mask"]) != expected or table["mask"].duplicated().any():
        return {"pass": False, "integrity_failure": True, "reason": "invalid gate-4 mask grid"}
    if not np.isfinite(table[["auc", "ess_fraction"]].to_numpy(float)).all():
        return {"pass": False, "integrity_failure": True, "reason": "nonfinite gate-4 diagnostic"}
    ok = (table["auc"] >= cfg["minimum_propensity_auc"]) & (
        table["ess_fraction"] >= cfg["minimum_ess_fraction"]
    )
    return {
        "pass": bool(ok.sum() >= cfg["minimum_predictable_masks"]),
        "integrity_failure": False,
        "qualifying_masks": int(ok.sum()),
        "required_masks": int(cfg["minimum_predictable_masks"]),
        "details": table.assign(qualifies=ok).to_dict("records"),
    }


def _winner(frame, value_col):
    values = frame.set_index("model")[value_col]
    minimum = values.min()
    winners = values.index[values == minimum].tolist()
    return winners[0] if len(winners) == 1 else None


def _strict_seeds(series):
    values = pd.to_numeric(series, errors="coerce").to_numpy(float)
    if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all():
        return None
    return values.astype(np.int64)


def gate5(rows, full, cfg):
    rows = pd.DataFrame(rows).copy()
    full = pd.DataFrame(full).copy()
    required_row_columns = {"mechanism", "truth_endpoint", "mask", "anchor", "seed", "gamma", "model", "cc"}
    required_full_columns = {"truth_endpoint", "model", "full_truth"}
    if not required_row_columns.issubset(rows.columns) or not required_full_columns.issubset(full.columns):
        return {"pass": False, "integrity_failure": True, "reason": "missing gate-5 columns"}
    models = list(cfg["models"])
    endpoints = list(cfg["truth_endpoints"])
    masks = list(cfg["mask_endpoints"])
    anchors = list(cfg["mar_anchors"])
    seeds = [int(x) for x in cfg["mask_seeds"]]
    expected_full = {(e, m) for e in endpoints for m in models}
    observed_full = set(map(tuple, full[["truth_endpoint", "model"]].itertuples(index=False, name=None)))
    if len(full) != len(expected_full) or observed_full != expected_full or full.duplicated(["truth_endpoint", "model"]).any():
        return {"pass": False, "integrity_failure": True, "reason": "invalid gate-5 full-truth grid"}
    if not np.isfinite(full["full_truth"].to_numpy(float)).all() or (full["full_truth"] <= 0).any():
        return {"pass": False, "integrity_failure": True, "reason": "nonpositive or nonfinite full truth"}
    actual = rows[rows["mechanism"] == "actual"].copy()
    mar = rows[rows["mechanism"] == "mar"].copy()
    expected_actual = {(e, mask, m) for e in endpoints for mask in masks for m in models}
    observed_actual = set(map(tuple, actual[["truth_endpoint", "mask", "model"]].itertuples(index=False, name=None)))
    if len(actual) != len(expected_actual) or observed_actual != expected_actual or actual.duplicated(["truth_endpoint", "mask", "model"]).any():
        return {"pass": False, "integrity_failure": True, "reason": "invalid gate-5 actual grid"}
    if not actual[["anchor", "seed", "gamma"]].isna().all().all() or not mar[["mask", "gamma"]].isna().all().all():
        return {"pass": False, "integrity_failure": True, "reason": "unexpected gate-5 mechanism key"}
    strict_mar_seeds = _strict_seeds(mar["seed"])
    if strict_mar_seeds is None:
        return {"pass": False, "integrity_failure": True, "reason": "invalid gate-5 seed"}
    mar["seed"] = strict_mar_seeds
    expected_mar = {(e, a, s, m) for e in endpoints for a in anchors for s in seeds for m in models}
    observed_mar = set((e, a, s, m) for e, a, s, m in mar[["truth_endpoint", "anchor", "seed", "model"]].itertuples(index=False, name=None))
    if len(mar) != len(expected_mar) or observed_mar != expected_mar or mar.duplicated(["truth_endpoint", "anchor", "seed", "model"]).any():
        return {"pass": False, "integrity_failure": True, "reason": "invalid gate-5 MAR grid"}
    if not np.isfinite(pd.concat([actual["cc"], mar["cc"]]).to_numpy(float)).all():
        return {"pass": False, "integrity_failure": True, "reason": "nonfinite gate-5 complete-case risk"}
    threshold = float(cfg["minimum_cc_relative_bias"])
    actual_details = []
    for keys, part in actual.groupby(["truth_endpoint", "mask"], sort=False):
        endpoint, mask = keys
        truth = full[full["truth_endpoint"] == endpoint]
        merged = part.merge(
            truth[["model", "full_truth"]].rename(
                columns={"full_truth": "canonical_full_truth"}
            ),
            on="model",
        )
        max_bias = float(
            (merged["cc"] - merged["canonical_full_truth"])
            .abs()
            .div(merged["canonical_full_truth"])
            .max()
        )
        truth_winner = _winner(truth, "full_truth")
        cc_winner = _winner(part, "cc")
        rank_change = truth_winner is not None and cc_winner is not None and cc_winner != truth_winner
        actual_details.append({"truth_endpoint": endpoint, "mask": mask, "max_relative_bias": max_bias, "rank_change": rank_change, "triggers": max_bias >= threshold or rank_change})
    mar_details = []
    for keys, part in mar.groupby(["truth_endpoint", "anchor"], sort=False):
        endpoint, anchor = keys
        truth = full[full["truth_endpoint"] == endpoint]
        merged = part.merge(
            truth[["model", "full_truth"]].rename(
                columns={"full_truth": "canonical_full_truth"}
            ),
            on="model",
        )
        med = merged.assign(
            relative_bias=(merged["cc"] - merged["canonical_full_truth"]).abs()
            / merged["canonical_full_truth"]
        ).groupby("model")["relative_bias"].median()
        max_bias = float(med.max())
        truth_winner = _winner(truth, "full_truth")
        changes = 0
        if truth_winner is not None:
            for _, seed_part in part.groupby("seed"):
                cc_winner = _winner(seed_part, "cc")
                changes += int(cc_winner is not None and cc_winner != truth_winner)
        trigger = max_bias >= threshold or changes >= 21
        mar_details.append({"truth_endpoint": endpoint, "anchor": anchor, "max_model_median_relative_bias": max_bias, "winner_changes": changes, "triggers": trigger})
    passed = any(x["triggers"] for x in actual_details + mar_details)
    return {"pass": bool(passed), "integrity_failure": False, "actual": actual_details, "mar": mar_details}


def gate6(rows, cfg):
    table = pd.DataFrame(rows).copy()
    required = ["cc", "aipw", "lcl", "ucl", "full_truth", "ess_fraction"]
    keys = ["mechanism", "truth_endpoint", "anchor", "seed", "model"]
    if not all(c in table for c in required + keys + ["mask", "gamma"]):
        return {"pass": False, "row_count": int(len(table)), "integrity_failure": True, "reason": "missing gate-6 columns"}
    models = list(cfg["models"])
    endpoints = list(cfg["truth_endpoints"])
    anchors = list(cfg["mar_anchors"])
    seeds = [int(x) for x in cfg["mask_seeds"]]
    strict_table_seeds = _strict_seeds(table["seed"])
    if strict_table_seeds is None:
        return {"pass": False, "row_count": int(len(table)), "integrity_failure": True, "reason": "invalid gate-6 seed"}
    table["seed"] = strict_table_seeds
    expected = {(e, a, s, m) for e in endpoints for a in anchors for s in seeds for m in models}
    observed = set((e, a, s, m) for mechanism, e, a, s, m in table[keys].itertuples(index=False, name=None) if mechanism == "mar")
    valid_grid = len(table) == len(expected) == 480 and observed == expected and not table.duplicated(["truth_endpoint", "anchor", "seed", "model"]).any() and (table["mechanism"] == "mar").all()
    finite = np.isfinite(table[required].to_numpy(dtype=float)).all()
    truth_constant = bool((table.groupby(["truth_endpoint", "model"])["full_truth"].nunique() == 1).all())
    mechanism_keys_empty = bool(table[["mask", "gamma"]].isna().all().all())
    if not valid_grid or not finite or not truth_constant or not mechanism_keys_empty:
        return {"pass": False, "row_count": int(len(table)), "integrity_failure": True, "reason": "invalid or nonfinite gate-6 MAR grid"}
    cc_bias = np.abs(table["cc"] - table["full_truth"])
    aipw_bias = np.abs(table["aipw"] - table["full_truth"])
    denominator = float(np.median(cc_bias))
    if denominator == 0:
        return {"pass": False, "row_count": 480, "integrity_failure": True, "reason": "zero median CC bias"}
    reduction = float(1 - np.median(aipw_bias) / denominator)
    covered = (table["lcl"] <= table["full_truth"]) & (table["ucl"] >= table["full_truth"])
    coverage_count = int(covered.sum())
    coverage = float(coverage_count / 480)
    width = table["ucl"] - table["lcl"]
    if (width < 0).any():
        return {"pass": False, "row_count": 480, "integrity_failure": True, "reason": "negative interval width"}
    ess_ok = bool((table["ess_fraction"] >= cfg["minimum_ess_fraction"]).all())
    passed = reduction >= cfg["minimum_aipw_bias_reduction"] and coverage >= cfg["minimum_aipw_ci_coverage"] and ess_ok
    return {"pass": bool(passed), "row_count": 480, "integrity_failure": False, "bias_reduction": reduction, "coverage": coverage, "coverage_count": coverage_count, "all_ess_acceptable": ess_ok, "minimum_ess_fraction": float(table["ess_fraction"].min()), "interval_width_median": float(width.median()), "interval_width_q05": float(width.quantile(0.05)), "interval_width_q95": float(width.quantile(0.95))}
