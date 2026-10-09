import argparse
import hashlib
import json
import platform
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import rdkit
import sklearn
import xgboost
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from xgboost import XGBRegressor

from pilot_estimators import (
    crossfit_estimators,
    crossfit_propensity,
    gate4,
    gate5,
    gate6,
    mar_mask,
    mcar_mask,
    mnar_mask,
    prepare_folds,
    scaffold_groups,
)


ENDPOINT_COLUMNS = {
    "LogD": "LogD",
    "KSOL": "KSOL",
    "HLM": "HLM CLint",
    "MLM": "MLM CLint",
    "Caco_Papp": "Caco-2 Permeability Papp A>B",
    "Caco_Efflux": "Caco-2 Permeability Efflux",
    "MPPB": "MPPB",
    "MBPB": "MBPB",
    "MGMB": "MGMB",
}
PRIMARY_MODELS = ["ecfp_rf", "ecfp_desc_xgb"]
NON_DUMMY_MODELS = ["ecfp_ridge", "ecfp_rf", "ecfp_desc_xgb"]


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding="utf-8")


def raw_state(value):
    text = str(value).strip()
    if not text or text.lower() in {"nan", "na", "n/a", "none"}:
        return "blank"
    try:
        number = float(text)
        return "numeric" if np.isfinite(number) else "blank"
    except ValueError:
        return "censored" if re.match(r"^[<>≤≥]\s*[-+]?\d", text) else "other"


def transform_y(values, endpoint):
    values = np.asarray(values, dtype=float)
    if endpoint == "KSOL":
        if np.any(values <= -1):
            raise ValueError("KSOL contains values outside log10(y+1) domain")
        return np.log10(values + 1.0)
    return values


def load_data(root, cfg):
    raw_dir = root / "data" / "raw"
    observed_hashes = {}
    for name, expected in cfg["dataset"]["files"].items():
        path = raw_dir / name
        observed = file_hash(path)
        if observed != expected:
            raise RuntimeError(f"hash mismatch for {name}: {observed}")
        observed_hashes[name] = observed
    train = pd.read_csv(raw_dir / "expansion_data_train.csv")
    test = pd.read_csv(raw_dir / "expansion_data_test.csv")
    raw = pd.read_csv(
        raw_dir / "expansion_data_raw.csv", dtype=str, keep_default_na=False
    )
    return train, test, raw, observed_hashes


def reconcile_sources(train, test, raw):
    id_col = "Molecule Name"
    if train[id_col].duplicated().any() or test[id_col].duplicated().any() or raw[id_col].duplicated().any():
        raise RuntimeError("duplicate molecule identifiers")
    if set(train[id_col]) & set(test[id_col]):
        raise RuntimeError("official train/test identifier overlap")
    ml_ids = set(train[id_col]) | set(test[id_col])
    raw_ids = set(raw[id_col])
    if not ml_ids.issubset(raw_ids):
        raise RuntimeError("ML-ready molecule missing from raw source")
    raw_by_id = raw.set_index(id_col)
    rows = []
    exact_numeric_match = True
    for phase, frame in [("train", train), ("test", test)]:
        for endpoint, column in ENDPOINT_COLUMNS.items():
            clean_numeric = pd.to_numeric(frame[column], errors="coerce").notna().to_numpy()
            states = frame[id_col].map(raw_by_id[column]).map(raw_state)
            exact_numeric_match &= bool(np.array_equal(clean_numeric, states.eq("numeric").to_numpy()))
            counts = states.value_counts().to_dict()
            rows.append(
                {
                    "phase": phase,
                    "endpoint": endpoint,
                    "rows": int(len(frame)),
                    "numeric": int(counts.get("numeric", 0)),
                    "censored": int(counts.get("censored", 0)),
                    "blank": int(counts.get("blank", 0)),
                    "other": int(counts.get("other", 0)),
                    "clean_numeric": int(clean_numeric.sum()),
                    "coverage": float(clean_numeric.mean()),
                }
            )
    audit = {
        "train_rows": int(len(train)),
        "test_rows": int(len(test)),
        "ml_ready_rows": int(len(train) + len(test)),
        "raw_rows": int(len(raw)),
        "raw_only_rows": int(len(raw_ids - ml_ids)),
        "identifier_overlap": 0,
        "raw_numeric_equals_clean_numeric": bool(exact_numeric_match),
    }
    return pd.DataFrame(rows), audit


def build_features(smiles, fp_cfg):
    descriptor_registry = sorted(Descriptors._descList, key=lambda item: item[0])
    descriptor_names = [name for name, _ in descriptor_registry]
    if len(descriptor_names) != len(set(descriptor_names)):
        raise RuntimeError("duplicate RDKit descriptor names")
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=int(fp_cfg["radius"]),
        fpSize=int(fp_cfg["n_bits"]),
        includeChirality=bool(fp_cfg["use_chirality"]),
    )
    fps = np.empty((len(smiles), int(fp_cfg["n_bits"])), dtype=np.float32)
    desc = np.empty((len(smiles), len(descriptor_registry)), dtype=np.float64)
    bitvectors = []
    canonical_iso = []
    canonical_noiso = []
    valid = []
    salts = []
    stereo = []
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(str(smi))
        valid.append(mol is not None)
        if mol is None:
            raise ValueError(f"invalid SMILES at row {i}")
        fp = generator.GetFingerprint(mol)
        arr = np.zeros(int(fp_cfg["n_bits"]), dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fps[i] = arr
        bitvectors.append(fp)
        desc[i] = [func(mol) for _, func in descriptor_registry]
        canonical_iso.append(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True))
        canonical_noiso.append(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False))
        salts.append("." in str(smi))
        stereo.append(bool(Chem.FindMolChiralCenters(mol, includeUnassigned=True)))
    if not np.isfinite(desc).all():
        raise RuntimeError("non-finite RDKit descriptor")
    return {
        "fps": fps,
        "desc": desc,
        "bitvectors": bitvectors,
        "descriptor_names": descriptor_names,
        "canonical_iso": np.asarray(canonical_iso, dtype=object),
        "canonical_noiso": np.asarray(canonical_noiso, dtype=object),
        "valid": np.asarray(valid, dtype=bool),
        "salts": np.asarray(salts, dtype=bool),
        "stereo": np.asarray(stereo, dtype=bool),
    }


def max_train_similarity(test_fps, train_fps):
    return np.asarray(
        [max(DataStructs.BulkTanimotoSimilarity(fp, train_fps)) for fp in test_fps],
        dtype=np.float64,
    )


def normalized_identifier_order(names):
    values = pd.Series(names, dtype=str).str.extract(r"(\d+)$", expand=False)
    if values.isna().any():
        raise RuntimeError("unparseable numeric identifier order")
    values = values.astype(float).to_numpy()
    span = values.max() - values.min()
    if span == 0:
        raise RuntimeError("constant numeric identifier order")
    return (values - values.min()) / span


def structure_audit(train_features, test_features):
    train_iso = set(train_features["canonical_iso"])
    test_iso = set(test_features["canonical_iso"])
    train_noiso = set(train_features["canonical_noiso"])
    test_noiso = set(test_features["canonical_noiso"])
    return {
        "valid_train": int(train_features["valid"].sum()),
        "valid_test": int(test_features["valid"].sum()),
        "train_exact_duplicates": int(len(train_features["canonical_iso"]) - len(train_iso)),
        "test_exact_duplicates": int(len(test_features["canonical_iso"]) - len(test_iso)),
        "train_test_isomeric_overlap": int(len(train_iso & test_iso)),
        "train_test_nonisomeric_overlap": int(len(train_noiso & test_noiso)),
        "train_salts": int(train_features["salts"].sum()),
        "test_salts": int(test_features["salts"].sum()),
        "train_stereo": int(train_features["stereo"].sum()),
        "test_stereo": int(test_features["stereo"].sum()),
    }


def fit_models(train, test, train_features, test_features, endpoint, cfg, smoke=False):
    column = ENDPOINT_COLUMNS[endpoint]
    train_numeric = pd.to_numeric(train[column], errors="coerce").to_numpy(float)
    test_numeric = pd.to_numeric(test[column], errors="coerce").to_numpy(float)
    train_idx = np.flatnonzero(np.isfinite(train_numeric))
    if smoke:
        train_idx = train_idx[: min(1200, len(train_idx))]
    y_train = transform_y(train_numeric[train_idx], endpoint)
    predictions = {}
    predictions["median"] = np.full(len(test), np.median(y_train), dtype=float)
    ridge = Ridge(alpha=float(cfg["models"]["ecfp_ridge"]["alpha"]), solver="lsqr", tol=1e-6)
    ridge.fit(train_features["fps"][train_idx], y_train)
    predictions["ecfp_ridge"] = ridge.predict(test_features["fps"])
    rf_cfg = cfg["models"]["ecfp_rf"]
    rf = RandomForestRegressor(
        n_estimators=32 if smoke else int(rf_cfg["n_estimators"]),
        min_samples_leaf=int(rf_cfg["min_samples_leaf"]),
        max_features=rf_cfg["max_features"],
        random_state=int(cfg["seed"]),
        n_jobs=-1,
    )
    rf.fit(train_features["fps"][train_idx], y_train)
    predictions["ecfp_rf"] = rf.predict(test_features["fps"])
    x_train = np.hstack([train_features["fps"][train_idx], train_features["desc"][train_idx]]).astype(np.float32)
    x_test = np.hstack([test_features["fps"], test_features["desc"]]).astype(np.float32)
    xgb_cfg = cfg["models"]["ecfp_desc_xgb"]
    xgb = XGBRegressor(
        objective="reg:squarederror",
        n_estimators=32 if smoke else int(xgb_cfg["n_estimators"]),
        max_depth=int(xgb_cfg["max_depth"]),
        learning_rate=float(xgb_cfg["learning_rate"]),
        subsample=float(xgb_cfg["subsample"]),
        colsample_bytree=float(xgb_cfg["colsample_bytree"]),
        reg_lambda=float(xgb_cfg["reg_lambda"]),
        random_state=int(cfg["seed"]),
        tree_method="hist",
        n_jobs=-1,
        verbosity=0,
    )
    xgb.fit(x_train, y_train)
    predictions["ecfp_desc_xgb"] = xgb.predict(x_test)
    if not all(np.isfinite(x).all() for x in predictions.values()):
        raise RuntimeError(f"non-finite prediction for {endpoint}")
    return predictions, train_numeric, test_numeric, int(len(train_idx))


def model_metrics(endpoint, predictions, train_numeric, test_numeric, fitted_train_n):
    observed = np.isfinite(test_numeric)
    y = transform_y(test_numeric[observed], endpoint)
    records = []
    for model, all_prediction in predictions.items():
        prediction = np.asarray(all_prediction)[observed]
        loss = np.abs(y - prediction)
        rho = spearmanr(y, prediction).statistic
        records.append(
            {
                "truth_endpoint": endpoint,
                "model": model,
                "train_n": int(fitted_train_n),
                "test_n": int(observed.sum()),
                "test_coverage": float(observed.mean()),
                "full_truth": float(loss.mean()),
                "r2": float(r2_score(y, prediction)),
                "spearman": float(rho),
            }
        )
    frame = pd.DataFrame(records)
    denominator = float(frame.loc[frame["model"] == "median", "full_truth"].iloc[0])
    frame["normalized_mae"] = frame["full_truth"] / denominator
    return frame


def cluster_se(values, groups):
    centered = np.asarray(values, dtype=float) - np.mean(values)
    codes = pd.factorize(np.asarray(groups, dtype=object), sort=False)[0]
    sums = np.bincount(codes, weights=centered)
    count = len(sums)
    if count < 2:
        raise RuntimeError("fewer than two scaffold clusters")
    return float(np.sqrt((count / (count - 1)) * np.square(sums).sum()) / len(values))


def add_scenario(
    risk_rows,
    pair_rows,
    endpoint,
    mechanism,
    mask_name,
    anchor,
    seed,
    gamma,
    r,
    losses,
    full_truth,
    x_obs,
    groups,
    prepared,
    cfg,
    generator_mean,
    mnar_shift,
):
    result = crossfit_estimators(
        x_obs,
        groups,
        r,
        losses,
        int(cfg["crossfit_folds"]),
        cfg["propensity_model"],
        cfg["loss_nuisance"],
        cfg["propensity_clip"],
        prepared=prepared,
    )
    base = {
        "truth_endpoint": endpoint,
        "mechanism": mechanism,
        "mask": mask_name,
        "anchor": anchor,
        "seed": seed,
        "gamma": gamma,
        "n_total": int(len(r)),
        "n_observed": int(np.sum(r)),
        "observed_fraction": float(np.mean(r)),
        "generator_mean_probability": generator_mean,
        "mnar_intercept_shift": mnar_shift,
    }
    for j, model in enumerate(PRIMARY_MODELS):
        estimate = result["estimates"][j]
        risk_rows.append(
            {
                **base,
                "model": model,
                "full_truth": float(full_truth[model]),
                **{key: value for key, value in estimate.items() if key != "column"},
                "cc_abs_bias": abs(estimate["cc"] - full_truth[model]),
                "aipw_abs_bias": abs(estimate["aipw"] - full_truth[model]),
                "ci_width": estimate["ucl"] - estimate["lcl"],
            }
        )
    p = result["p"]
    mhat = result["mhat"]
    pseudo = mhat + np.asarray(r)[:, None] * (losses - mhat) / p[:, None]
    delta_loss = losses[:, 0] - losses[:, 1]
    weights = np.asarray(r, dtype=float) / p
    delta_pseudo = pseudo[:, 0] - pseudo[:, 1]
    delta_aipw = float(delta_pseudo.mean())
    delta_se = cluster_se(delta_pseudo, groups)
    observed = np.asarray(r, dtype=bool)
    pair_rows.append(
        {
            **base,
            "model_a": PRIMARY_MODELS[0],
            "model_b": PRIMARY_MODELS[1],
            "full_truth_delta": float(full_truth[PRIMARY_MODELS[0]] - full_truth[PRIMARY_MODELS[1]]),
            "cc_delta": float(delta_loss[observed].mean()),
            "ipw_delta": float(np.sum(weights * delta_loss) / weights.sum()),
            "aipw_delta": delta_aipw,
            "aipw_delta_se": delta_se,
            "aipw_delta_lcl": delta_aipw - 1.96 * delta_se,
            "aipw_delta_ucl": delta_aipw + 1.96 * delta_se,
        }
    )
    return result


def gate_configuration(cfg):
    return {
        **cfg["gates"],
        "models": list(cfg["primary_models"]),
        "truth_endpoints": list(cfg["truth_endpoints"]),
        "mask_endpoints": list(cfg["mask_endpoints"]),
        "mar_anchors": list(cfg["mar_anchors"]),
        "mask_seeds": list(cfg["mask_seeds"]),
    }


def validate_scenario_grid(table, cfg, paired=False):
    table = pd.DataFrame(table).copy()
    key_columns = ["mechanism", "truth_endpoint", "mask", "anchor", "seed", "gamma"]
    if not set(key_columns).issubset(table.columns):
        return False, "missing scenario key column"
    if not paired and "model" not in table:
        return False, "missing scenario model column"
    observed = []
    for row in table.itertuples(index=False):
        values = row._asdict()
        seed = values["seed"]
        if pd.isna(seed):
            seed = None
        elif not np.isfinite(float(seed)) or float(seed) != np.floor(float(seed)):
            return False, "invalid scenario seed"
        else:
            seed = int(seed)
        gamma = values["gamma"]
        gamma = None if pd.isna(gamma) else float(gamma)
        key = (
            values["mechanism"],
            values["truth_endpoint"],
            None if pd.isna(values["mask"]) else values["mask"],
            None if pd.isna(values["anchor"]) else values["anchor"],
            seed,
            gamma,
        )
        observed.append(key if paired else (*key, values["model"]))
    expected = []
    models = [None] if paired else list(cfg["primary_models"])
    for endpoint in cfg["truth_endpoints"]:
        for model in models:
            suffix = () if paired else (model,)
            for mask_name in cfg["mask_endpoints"]:
                expected.append(("actual", endpoint, mask_name, None, None, None, *suffix))
                for seed in cfg["mask_seeds"]:
                    expected.append(("mcar", endpoint, mask_name, None, int(seed), None, *suffix))
            for anchor in cfg["mar_anchors"]:
                for seed in cfg["mask_seeds"]:
                    expected.append(("mar", endpoint, None, anchor, int(seed), None, *suffix))
                    for gamma in cfg["mnar_gamma"]:
                        expected.append(("mnar", endpoint, None, anchor, int(seed), float(gamma), *suffix))
    if len(observed) != len(expected) or len(set(observed)) != len(observed) or set(observed) != set(expected):
        return False, "scenario Cartesian grid mismatch"
    if paired and (
        not (table["model_a"] == PRIMARY_MODELS[0]).all()
        or not (table["model_b"] == PRIMARY_MODELS[1]).all()
    ):
        return False, "paired model identity mismatch"
    return True, None


def prepare_study(root, cfg):
    if file_hash(root / "protocol.md") != cfg["protocol_sha256"]:
        raise RuntimeError("protocol/config hash mismatch")
    if rdkit.__version__ != cfg["rdkit_descriptors"]["rdkit_version"]:
        raise RuntimeError(f"RDKit version mismatch: {rdkit.__version__}")
    train, test, raw, data_hashes = load_data(root, cfg)
    cohort_table, cohort_summary = reconcile_sources(train, test, raw)
    all_smiles = pd.concat([train["SMILES"], test["SMILES"]], ignore_index=True)
    all_features = build_features(all_smiles.tolist(), cfg["fingerprint"])
    split = len(train)
    train_features = {key: value[:split] if isinstance(value, np.ndarray) else value[:split] for key, value in all_features.items() if key not in {"descriptor_names"}}
    test_features = {key: value[split:] if isinstance(value, np.ndarray) else value[split:] for key, value in all_features.items() if key not in {"descriptor_names"}}
    train_features["descriptor_names"] = all_features["descriptor_names"]
    test_features["descriptor_names"] = all_features["descriptor_names"]
    if len(all_features["descriptor_names"]) != int(cfg["rdkit_descriptors"]["expected_count"]):
        raise RuntimeError("RDKit descriptor count mismatch")
    similarity = max_train_similarity(test_features["bitvectors"], train_features["bitvectors"])
    identifier_order = normalized_identifier_order(test["Molecule Name"])
    x_obs = np.hstack(
        [
            test_features["fps"],
            test_features["desc"],
            similarity[:, None],
            identifier_order[:, None],
        ]
    ).astype(np.float32)
    expected_observation_shape = (
        len(test),
        int(cfg["fingerprint"]["n_bits"])
        + int(cfg["rdkit_descriptors"]["expected_count"])
        + 2,
    )
    if x_obs.shape != expected_observation_shape or not np.isfinite(x_obs).all():
        raise RuntimeError("non-finite observation feature")
    if list(cfg["primary_models"]) != PRIMARY_MODELS:
        raise RuntimeError("primary model order mismatch")
    groups = scaffold_groups(test["SMILES"].tolist())
    structures = structure_audit(train_features, test_features)
    availability = {
        endpoint: pd.to_numeric(test[ENDPOINT_COLUMNS[endpoint]], errors="coerce").notna().to_numpy(dtype=int)
        for endpoint in cfg["mask_endpoints"]
    }
    return {
        "train": train,
        "test": test,
        "raw": raw,
        "data_hashes": data_hashes,
        "cohort_table": cohort_table,
        "cohort_summary": cohort_summary,
        "train_features": train_features,
        "test_features": test_features,
        "descriptor_names": all_features["descriptor_names"],
        "similarity": similarity,
        "identifier_order": identifier_order,
        "x_obs": x_obs,
        "groups": groups,
        "structures": structures,
        "availability": availability,
        "expected_observation_shape": list(expected_observation_shape),
    }


def save_common_outputs(output, study):
    output.mkdir(parents=True, exist_ok=True)
    study["cohort_table"].to_csv(output / "cohort_audit.csv", index=False)
    write_json(output / "cohort_summary.json", study["cohort_summary"])
    write_json(output / "structure_audit.json", study["structures"])


def run_smoke(root, cfg, study, output):
    predictions, train_y, test_y, fitted_n = fit_models(
        study["train"],
        study["test"],
        study["train_features"],
        study["test_features"],
        "LogD",
        cfg,
        smoke=True,
    )
    metrics = model_metrics("LogD", predictions, train_y, test_y, fitted_n)
    truth_all = np.flatnonzero(np.isfinite(test_y))
    truth_idx = truth_all[: min(600, len(truth_all))]
    y = transform_y(test_y[truth_idx], "LogD")
    losses = np.column_stack(
        [np.abs(y - predictions[model][truth_idx]) for model in PRIMARY_MODELS]
    )
    groups = study["groups"][truth_idx]
    x_obs = study["x_obs"][truth_idx]
    prepared = prepare_folds(x_obs, groups, int(cfg["crossfit_folds"]))
    r = study["availability"]["HLM"][truth_idx]
    full_truth = {model: float(losses[:, j].mean()) for j, model in enumerate(PRIMARY_MODELS)}
    risk_rows = []
    pair_rows = []
    result = add_scenario(
        risk_rows,
        pair_rows,
        "LogD",
        "actual",
        "HLM",
        None,
        None,
        None,
        r,
        losses,
        full_truth,
        x_obs,
        groups,
        prepared,
        cfg,
        None,
        None,
    )
    risk = pd.DataFrame(risk_rows)
    pair = pd.DataFrame(pair_rows)
    smoke_rows = pd.DataFrame(
        {
            "row_index": truth_idx,
            "molecule_id": study["test"].iloc[truth_idx]["Molecule Name"].to_numpy(),
            "scaffold_sha256": [
                hashlib.sha256(str(value).encode("utf-8")).hexdigest()
                for value in groups
            ],
            "y": y,
            "r": r,
            "p_raw": result["p_raw"],
            "p": result["p"],
        }
    )
    for j, model in enumerate(PRIMARY_MODELS):
        smoke_rows[f"loss_{model}"] = losses[:, j]
        smoke_rows[f"mhat_{model}"] = result["mhat"][:, j]
        smoke_rows[f"pseudo_{model}"] = (
            result["mhat"][:, j]
            + r * (losses[:, j] - result["mhat"][:, j]) / result["p"]
        )
    metrics.to_csv(output / "smoke_model_metrics.csv", index=False)
    risk.to_csv(output / "smoke_scenario.csv", index=False)
    pair.to_csv(output / "smoke_pairwise.csv", index=False)
    smoke_rows.to_csv(output / "smoke_row_audit.csv", index=False)
    formula_differences = []
    weights = r / result["p"]
    for j, model in enumerate(PRIMARY_MODELS):
        recorded = risk[risk["model"] == model].iloc[0]
        formula_differences.extend(
            [
                abs(recorded["cc"] - losses[r.astype(bool), j].mean()),
                abs(recorded["ipw"] - np.sum(weights * losses[:, j]) / weights.sum()),
                abs(recorded["aipw"] - smoke_rows[f"pseudo_{model}"].mean()),
            ]
        )
    recorded_pair = pair.iloc[0]
    loss_delta = losses[:, 0] - losses[:, 1]
    pseudo_delta = (
        smoke_rows[f"pseudo_{PRIMARY_MODELS[0]}"]
        - smoke_rows[f"pseudo_{PRIMARY_MODELS[1]}"]
    )
    formula_differences.extend(
        [
            abs(recorded_pair["full_truth_delta"] - loss_delta.mean()),
            abs(recorded_pair["cc_delta"] - loss_delta[r.astype(bool)].mean()),
            abs(recorded_pair["ipw_delta"] - np.sum(weights * loss_delta) / weights.sum()),
            abs(recorded_pair["aipw_delta"] - pseudo_delta.mean()),
        ]
    )
    pair_numeric = [
        "full_truth_delta",
        "cc_delta",
        "ipw_delta",
        "aipw_delta",
        "aipw_delta_se",
        "aipw_delta_lcl",
        "aipw_delta_ucl",
    ]
    smoke_ok = bool(
        len(risk) == 2
        and len(pair) == 1
        and np.isfinite(
            risk[["cc", "ipw", "aipw", "se", "ess_fraction"]].to_numpy()
        ).all()
        and np.isfinite(pair[pair_numeric].to_numpy()).all()
        and max(formula_differences) < 1e-12
    )
    checks = {
        "status": "PASS" if smoke_ok else "FAIL",
        "non_adjudicative": True,
        "truth_endpoint": "LogD",
        "mask": "HLM",
        "train_rows_used": int(fitted_n),
        "truth_rows_used": int(len(truth_idx)),
        "risk_rows": int(len(risk)),
        "pairwise_rows": int(len(pair)),
        "p_raw_finite": bool(np.isfinite(result["p_raw"]).all()),
        "p_clipped_range": [float(result["p"].min()), float(result["p"].max())],
        "mhat_finite": bool(np.isfinite(result["mhat"]).all()),
        "formula_recalculation_max_abs_difference": float(max(formula_differences)),
        "smoke_overrides": {"rf_n_estimators": 32, "xgb_n_estimators": 32, "maximum_train_rows": 1200, "maximum_truth_rows": 600},
    }
    write_json(output / "smoke_audit.json", checks)
    return checks


def run_full(root, cfg, study, output):
    model_frames = []
    predictions_by_endpoint = {}
    test_y_by_endpoint = {}
    fitted_train_n = {}
    for endpoint in cfg["truth_endpoints"]:
        print(f"MODEL {endpoint}", flush=True)
        predictions, train_y, test_y, fitted_n = fit_models(
            study["train"],
            study["test"],
            study["train_features"],
            study["test_features"],
            endpoint,
            cfg,
            smoke=False,
        )
        predictions_by_endpoint[endpoint] = predictions
        test_y_by_endpoint[endpoint] = test_y
        fitted_train_n[endpoint] = fitted_n
        model_frames.append(model_metrics(endpoint, predictions, train_y, test_y, fitted_n))
    model_table = pd.concat(model_frames, ignore_index=True)
    model_table.to_csv(output / "model_metrics.csv", index=False)

    prepared5 = prepare_folds(study["x_obs"], study["groups"], int(cfg["availability_auc_folds"]))
    diagnostics = []
    propensity = {}
    for mask_name in cfg["mask_endpoints"]:
        print(f"PROPENSITY {mask_name}", flush=True)
        result = crossfit_propensity(
            study["x_obs"],
            study["groups"],
            study["availability"][mask_name],
            int(cfg["availability_auc_folds"]),
            cfg["propensity_model"],
            cfg["propensity_clip"],
            prepared=prepared5,
        )
        propensity[mask_name] = result
        diagnostics.append(
            {
                "mask": mask_name,
                "coverage": float(study["availability"][mask_name].mean()),
                "observed_n": int(study["availability"][mask_name].sum()),
                "auc": result["auc"],
                "clip_fraction": result["clip_fraction"],
                "ess": result["ess"],
                "ess_fraction": result["ess_fraction"],
            }
        )
    diagnostics_table = pd.DataFrame(diagnostics)
    diagnostics_table.to_csv(output / "actual_mask_diagnostics.csv", index=False)

    risk_rows = []
    pair_rows = []
    for endpoint in cfg["truth_endpoints"]:
        test_y = test_y_by_endpoint[endpoint]
        truth_idx = np.flatnonzero(np.isfinite(test_y))
        y = transform_y(test_y[truth_idx], endpoint)
        losses = np.column_stack(
            [
                np.abs(y - predictions_by_endpoint[endpoint][model][truth_idx])
                for model in PRIMARY_MODELS
            ]
        )
        full_truth = {
            model: float(losses[:, j].mean()) for j, model in enumerate(PRIMARY_MODELS)
        }
        x_obs = study["x_obs"][truth_idx]
        groups = study["groups"][truth_idx]
        prepared3 = prepare_folds(x_obs, groups, int(cfg["crossfit_folds"]))
        for mask_name in cfg["mask_endpoints"]:
            r = study["availability"][mask_name][truth_idx]
            add_scenario(
                risk_rows, pair_rows, endpoint, "actual", mask_name, None, None, None,
                r, losses, full_truth, x_obs, groups, prepared3, cfg, None, None,
            )
        print(f"SCENARIOS {endpoint} actual", flush=True)
        for mask_name in cfg["mask_endpoints"]:
            coverage = float(study["availability"][mask_name][truth_idx].mean())
            for seed in cfg["mask_seeds"]:
                r = mcar_mask(len(truth_idx), coverage, int(seed))
                add_scenario(
                    risk_rows, pair_rows, endpoint, "mcar", mask_name, None, int(seed), None,
                    r, losses, full_truth, x_obs, groups, prepared3, cfg, coverage, None,
                )
            print(f"SCENARIOS {endpoint} mcar {mask_name}", flush=True)
        for anchor in cfg["mar_anchors"]:
            base_p = propensity[anchor]["clipped"][truth_idx]
            for seed in cfg["mask_seeds"]:
                r = mar_mask(base_p, int(seed))
                add_scenario(
                    risk_rows, pair_rows, endpoint, "mar", None, anchor, int(seed), None,
                    r, losses, full_truth, x_obs, groups, prepared3, cfg, float(base_p.mean()), None,
                )
            print(f"SCENARIOS {endpoint} mar {anchor}", flush=True)
            for gamma in cfg["mnar_gamma"]:
                for seed in cfg["mask_seeds"]:
                    r, probability, shift = mnar_mask(base_p, y, float(gamma), int(seed))
                    add_scenario(
                        risk_rows, pair_rows, endpoint, "mnar", None, anchor, int(seed), float(gamma),
                        r, losses, full_truth, x_obs, groups, prepared3, cfg, float(probability.mean()), float(shift),
                    )
                print(f"SCENARIOS {endpoint} mnar {anchor} gamma={gamma}", flush=True)

    risk_table = pd.DataFrame(risk_rows)
    pair_table = pd.DataFrame(pair_rows)
    risk_table.to_csv(output / "scenario_results.csv", index=False)
    pair_table.to_csv(output / "pairwise_results.csv", index=False)
    full_primary = model_table[model_table["model"].isin(PRIMARY_MODELS)][
        ["truth_endpoint", "model", "full_truth"]
    ].copy()
    risk_with_canonical = risk_table.merge(
        full_primary.rename(columns={"full_truth": "canonical_full_truth"}),
        on=["truth_endpoint", "model"],
        how="left",
        validate="many_to_one",
    )
    full_truth_equal = bool(
        len(risk_with_canonical) == len(risk_table)
        and risk_with_canonical["canonical_full_truth"].notna().all()
        and np.allclose(
            risk_with_canonical["full_truth"],
            risk_with_canonical["canonical_full_truth"],
            rtol=0,
            atol=1e-12,
        )
    )
    gate_cfg = gate_configuration(cfg)
    gate4_result = gate4(diagnostics, gate_cfg)
    gate5_result = gate5(risk_rows, full_primary, gate_cfg)
    mar_primary = risk_table[risk_table["mechanism"] == "mar"].copy()
    gate6_result = gate6(mar_primary, gate_cfg)

    coverage_by_endpoint = model_table.groupby("truth_endpoint")["test_coverage"].first()
    expected_cohort_keys = {
        (phase, endpoint)
        for phase in ["train", "test"]
        for endpoint in ENDPOINT_COLUMNS
    }
    observed_cohort_keys = set(
        map(
            tuple,
            study["cohort_table"][["phase", "endpoint"]].itertuples(
                index=False, name=None
            ),
        )
    )
    state_partition_ok = bool(
        (
            study["cohort_table"]["numeric"]
            + study["cohort_table"]["censored"]
            + study["cohort_table"]["blank"]
            == study["cohort_table"]["rows"]
        ).all()
        and (study["cohort_table"]["other"] == 0).all()
        and (
            study["cohort_table"]["clean_numeric"]
            == study["cohort_table"]["numeric"]
        ).all()
    )
    gate1_pass = bool(
        study["cohort_summary"]["train_rows"] == 5326
        and study["cohort_summary"]["test_rows"] == 2282
        and study["cohort_summary"]["ml_ready_rows"] == 7608
        and study["cohort_summary"]["raw_rows"] == 7618
        and study["cohort_summary"]["raw_only_rows"] == 10
        and study["cohort_summary"]["ml_ready_rows"] >= cfg["gates"]["minimum_parent_rows"]
        and study["cohort_summary"]["raw_numeric_equals_clean_numeric"]
        and len(study["cohort_table"]) == len(expected_cohort_keys)
        and observed_cohort_keys == expected_cohort_keys
        and state_partition_ok
    )
    gate2_pass = bool((coverage_by_endpoint >= cfg["gates"]["minimum_truth_coverage"]).all())
    performance = model_table[model_table["model"].isin(NON_DUMMY_MODELS)].pivot(
        index="model", columns="truth_endpoint", values="normalized_mae"
    )
    qualifying_models = performance.index[(performance < cfg["gates"]["maximum_normalized_mae"]).all(axis=1)].tolist()
    gate3_pass = len(qualifying_models) >= 2
    expected_risk_rows = len(cfg["truth_endpoints"]) * len(PRIMARY_MODELS) * (
        len(cfg["mask_endpoints"])
        + len(cfg["mask_endpoints"]) * len(cfg["mask_seeds"])
        + len(cfg["mar_anchors"]) * len(cfg["mask_seeds"])
        + len(cfg["mar_anchors"]) * len(cfg["mnar_gamma"]) * len(cfg["mask_seeds"])
    )
    required_risk_numeric = ["cc", "ipw", "aipw", "se", "lcl", "ucl", "ess", "ess_fraction", "full_truth"]
    required_pair_numeric = ["full_truth_delta", "cc_delta", "ipw_delta", "aipw_delta", "aipw_delta_se", "aipw_delta_lcl", "aipw_delta_ucl"]
    risk_grid_ok, risk_grid_reason = validate_scenario_grid(risk_table, cfg, paired=False)
    pair_grid_ok, pair_grid_reason = validate_scenario_grid(pair_table, cfg, paired=True)
    structural_ok = (
        study["structures"]["train_test_isomeric_overlap"] == 0
        and study["structures"]["train_test_nonisomeric_overlap"] == 0
        and study["structures"]["train_exact_duplicates"] == 0
        and study["structures"]["test_exact_duplicates"] == 0
        and study["structures"]["valid_train"] == len(study["train"])
        and study["structures"]["valid_test"] == len(study["test"])
    )
    finite_ok = bool(np.isfinite(risk_table[required_risk_numeric].to_numpy(float)).all())
    pair_finite_ok = bool(np.isfinite(pair_table[required_pair_numeric].to_numpy(float)).all())
    expected_deltas = {
        endpoint: float(
            full_primary[(full_primary["truth_endpoint"] == endpoint) & (full_primary["model"] == PRIMARY_MODELS[0])]["full_truth"].iloc[0]
            - full_primary[(full_primary["truth_endpoint"] == endpoint) & (full_primary["model"] == PRIMARY_MODELS[1])]["full_truth"].iloc[0]
        )
        for endpoint in cfg["truth_endpoints"]
    }
    pair_truth_equal = all(
        np.allclose(part["full_truth_delta"], expected_deltas[endpoint], rtol=0, atol=1e-12)
        for endpoint, part in pair_table.groupby("truth_endpoint")
    )
    gate_integrity_ok = not any(
        result.get("integrity_failure", False)
        for result in [gate4_result, gate5_result, gate6_result]
    )
    gate7_pass = bool(
        structural_ok
        and finite_ok
        and pair_finite_ok
        and len(risk_table) == expected_risk_rows
        and risk_grid_ok
        and pair_grid_ok
        and full_truth_equal
        and pair_truth_equal
        and gate_integrity_ok
        and len(study["descriptor_names"]) == cfg["rdkit_descriptors"]["expected_count"]
        and list(study["x_obs"].shape) == study["expected_observation_shape"]
        and list(cfg["primary_models"]) == PRIMARY_MODELS
        and rdkit.__version__ == cfg["rdkit_descriptors"]["rdkit_version"]
        and file_hash(root / "protocol.md") == cfg["protocol_sha256"]
    )
    gates = {
        "gate1": {
            "pass": gate1_pass,
            **study["cohort_summary"],
            "state_partition_ok": state_partition_ok,
            "cohort_grid_ok": observed_cohort_keys == expected_cohort_keys,
        },
        "gate2": {"pass": gate2_pass, "coverage": coverage_by_endpoint.to_dict()},
        "gate3": {"pass": gate3_pass, "qualifying_non_dummy_models": qualifying_models, "required": 2},
        "gate4": gate4_result,
        "gate5": gate5_result,
        "gate6": gate6_result,
        "gate7": {
            "pass": gate7_pass,
            "structural_ok": structural_ok,
            "finite_required_metrics": finite_ok,
            "finite_pairwise_metrics": pair_finite_ok,
            "expected_risk_rows": int(expected_risk_rows),
            "observed_risk_rows": int(len(risk_table)),
            "observed_pairwise_rows": int(len(pair_table)),
            "risk_grid_ok": risk_grid_ok,
            "risk_grid_reason": risk_grid_reason,
            "pair_grid_ok": pair_grid_ok,
            "pair_grid_reason": pair_grid_reason,
            "full_truth_matches_canonical_table": full_truth_equal,
            "pair_truth_matches_canonical_table": pair_truth_equal,
            "gate_input_integrity": gate_integrity_ok,
            "outcome_features_used": False,
            "observation_feature_shape": list(study["x_obs"].shape),
            "expected_observation_feature_shape": study["expected_observation_shape"],
            "primary_model_order_matches_config": list(cfg["primary_models"]) == PRIMARY_MODELS,
        },
    }
    critical = all(gates[name]["pass"] for name in ["gate1", "gate2", "gate3", "gate7"])
    decision = "GO_FULL" if critical and all(gates[name]["pass"] for name in ["gate4", "gate5", "gate6"]) else ("STOP_PILOT_NEGATIVE" if critical else "STOP_CRITICAL")
    decision_payload = {"protocol_id": cfg["protocol_id"], "decision": decision, "gates": gates}
    write_json(output / "gate_decision.json", decision_payload)

    row_audit = pd.DataFrame(
        {
            "row_index": np.arange(len(study["test"])),
            "molecule_id": study["test"]["Molecule Name"],
            "canonical_sha256": [hashlib.sha256(x.encode("utf-8")).hexdigest() for x in study["test_features"]["canonical_iso"]],
            "scaffold_sha256": [hashlib.sha256(str(x).encode("utf-8")).hexdigest() for x in study["groups"]],
            "max_train_similarity": study["similarity"],
            "normalized_identifier_order": study["identifier_order"],
        }
    )
    for mask_name in cfg["mask_endpoints"]:
        row_audit[f"available_{mask_name}"] = study["availability"][mask_name]
        row_audit[f"propensity_raw_{mask_name}"] = propensity[mask_name]["raw"]
        row_audit[f"propensity_clipped_{mask_name}"] = propensity[mask_name]["clipped"]
    predictions_dir = output / "predictions"
    predictions_dir.mkdir(exist_ok=True)
    row_audit.to_csv(predictions_dir / "row_audit.csv", index=False)
    prediction_rows = []
    for endpoint, model_predictions in predictions_by_endpoint.items():
        for model, values in model_predictions.items():
            prediction_rows.append(
                pd.DataFrame(
                    {
                        "row_index": np.arange(len(study["test"])),
                        "molecule_id": study["test"]["Molecule Name"],
                        "truth_endpoint": endpoint,
                        "model": model,
                        "prediction": values,
                    }
                )
            )
    pd.concat(prediction_rows, ignore_index=True).to_parquet(predictions_dir / "predictions.parquet", index=False)
    return decision_payload, {
        "model_table": model_table,
        "diagnostics_table": diagnostics_table,
        "risk_table": risk_table,
        "pair_table": pair_table,
        "fitted_train_n": fitted_train_n,
    }


def output_hashes(output):
    return {
        str(path.relative_to(output)).replace("\\", "/"): file_hash(path)
        for path in sorted(output.rglob("*"))
        if path.is_file() and path.name != "run_manifest.json"
    }


def write_summary(output, mode, result):
    if mode == "smoke":
        text = (
            "# Pilot smoke test\n\n"
            f"Status: `{result['status']}`. This was a non-adjudicative wiring test with reduced model sizes and row counts.\n"
        )
    else:
        lines = ["# Pilot gate decision", "", f"Decision: `{result['decision']}`", ""]
        for name, value in result["gates"].items():
            lines.append(f"- {name}: {'PASS' if value['pass'] else 'FAIL'}")
        lines.extend(
            [
                "",
                "This decision follows the frozen protocol without threshold changes. Real sparse-endpoint conclusions remain conditional on MAR and positivity; MNAR results are sensitivity analyses.",
            ]
        )
        text = "\n".join(lines) + "\n"
    (output / "summary.md").write_text(text, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full"], required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    cfg_path = root / "config" / "pilot.json"
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    default_output = root / "artifacts" / "experiment" / cfg["protocol_id"] / args.mode
    output = args.output.resolve() if args.output else default_output
    if output.exists() and any(output.iterdir()):
        raise RuntimeError(f"output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    started = time.time()
    started_utc = datetime.now(timezone.utc).isoformat()
    print(f"START mode={args.mode} protocol={cfg['protocol_sha256']}", flush=True)
    study = prepare_study(root, cfg)
    save_common_outputs(output, study)
    if args.mode == "smoke":
        result = run_smoke(root, cfg, study, output)
        details = {}
    else:
        result, details = run_full(root, cfg, study, output)
    write_summary(output, args.mode, result)
    manifest = {
        "protocol_id": cfg["protocol_id"],
        "mode": args.mode,
        "started_utc": started_utc,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "duration_seconds": float(time.time() - started),
        "command": [sys.executable, *sys.argv],
        "hashes": {
            "protocol.md": file_hash(root / "protocol.md"),
            "config/pilot.json": file_hash(cfg_path),
            "src/run_pilot.py": file_hash(Path(__file__)),
            "src/pilot_estimators.py": file_hash(Path(__file__).with_name("pilot_estimators.py")),
            "baseline/metrics.json": file_hash(root / "baseline" / "metrics.json"),
            **{f"data/raw/{key}": value for key, value in study["data_hashes"].items()},
        },
        "baseline_protocol_sha256": json.loads((root / "baseline" / "metrics.json").read_text(encoding="utf-8"))["protocol_sha256"],
        "current_protocol_sha256": cfg["protocol_sha256"],
        "dataset_revision": cfg["dataset"]["revision"],
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
            "rdkit": rdkit.__version__,
            "xgboost": xgboost.__version__,
            "platform": platform.platform(),
        },
        "features": {
            "fingerprint": cfg["fingerprint"],
            "rdkit_descriptor_source": cfg["rdkit_descriptors"]["source"],
            "rdkit_descriptor_count": len(study["descriptor_names"]),
            "rdkit_descriptor_names": study["descriptor_names"],
            "observation_feature_components": [
                "ECFP radius 2 2048-bit with chirality",
                "all sorted finite RDKit 2D descriptors",
                "maximum test-to-all-official-train ECFP Tanimoto similarity",
                "test-cohort min-max normalized numeric identifier suffix",
            ],
            "observation_feature_shape": list(study["x_obs"].shape),
            "observation_feature_dtype": str(study["x_obs"].dtype),
            "standard_scaler": "dense input; with_mean=True; fit separately within every nuisance-training fold",
            "scaffold_grouping": "chiral Bemis-Murcko; acyclic keys are singleton by fixed full-test row index",
        },
        "result": result,
        "details": {
            "fitted_train_n": details.get("fitted_train_n", {}),
            "risk_rows": int(len(details.get("risk_table", []))),
            "pairwise_rows": int(len(details.get("pair_table", []))),
        },
    }
    manifest["output_hashes_before_manifest"] = output_hashes(output)
    write_json(output / "run_manifest.json", manifest)
    print(f"DONE mode={args.mode} result={result.get('decision', result.get('status'))} seconds={manifest['duration_seconds']:.1f}", flush=True)


if __name__ == "__main__":
    main()
