import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupKFold


EXPECTED_CONFIG_SHA256 = "72dfed2fd134fff7b5e6ccc6d6f0300659def6bd0d6d6f510aec4372a929c203"
TRUSTED_FORENSIC_RUNNER_SHA256 = "d465ba4b162d64a52d807fa520efffd3a944a8dcfcb43d9289b6ce9bcc9044d8"
ENDPOINT_COLUMNS = {"LogD": "LogD", "KSOL": "KSOL"}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def vector_hash(values, dtype):
    return hashlib.sha256(np.asarray(values, dtype=dtype).tobytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def close(left, right, tolerance, label):
    left = float(left)
    right = float(right)
    require(np.isfinite(left) and np.isfinite(right), f"nonfinite {label}")
    require(abs(left - right) <= tolerance, f"{label}: {left} != {right}")


def compare_tree(expected, observed, tolerance, path="root"):
    if isinstance(expected, dict):
        require(isinstance(observed, dict), f"{path} is not an object")
        require(set(expected) == set(observed), f"{path} keys differ")
        for key in expected:
            compare_tree(expected[key], observed[key], tolerance, f"{path}.{key}")
        return
    if isinstance(expected, list):
        require(isinstance(observed, list), f"{path} is not a list")
        require(len(expected) == len(observed), f"{path} length differs")
        for index, (left, right) in enumerate(zip(expected, observed)):
            compare_tree(left, right, tolerance, f"{path}[{index}]")
        return
    if expected is None or isinstance(expected, (str, bool)):
        require(expected == observed, f"{path} differs")
        return
    if isinstance(expected, (int, float, np.integer, np.floating)):
        close(expected, observed, tolerance, path)
        return
    require(expected == observed, f"{path} differs")


def transformed(values, endpoint):
    values = np.asarray(values, dtype=float)
    return np.log10(values + 1.0) if endpoint == "KSOL" else values


def cluster_se(influence, groups):
    codes = pd.factorize(np.asarray(groups, dtype=object), sort=False)[0]
    sums = np.bincount(codes, weights=np.asarray(influence, dtype=float))
    count = len(sums)
    require(count > 1, "fewer than two scaffold clusters")
    return float(np.sqrt((count / (count - 1)) * np.square(sums).sum()) / len(influence))


def scaffold_contract(test):
    groups = []
    for index, smiles in enumerate(test["SMILES"]):
        molecule = Chem.MolFromSmiles(str(smiles))
        require(molecule is not None, f"invalid test SMILES at row {index}")
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(
            mol=molecule, includeChirality=True
        )
        groups.append(scaffold or f"ACYCLIC_SINGLETON_{index}")
    groups = np.asarray(groups, dtype=object)
    hashes = np.asarray(
        [hashlib.sha256(str(value).encode("utf-8")).hexdigest() for value in groups]
    )
    return groups, hashes


def fold_contract(groups, folds):
    split = GroupKFold(n_splits=int(folds)).split(
        np.zeros(len(groups)), groups=groups
    )
    fold_id = np.full(len(groups), -1, dtype=np.int64)
    for value, (_, valid) in enumerate(split):
        fold_id[valid] = value
    require(np.all(fold_id >= 0), "incomplete independent fold vector")
    return fold_id, vector_hash(fold_id, "<i8")


def expected_input_hashes(root, cfg, cfg_path):
    expected = {
        "artifacts/diagnostics/forensic_001/config.json": sha256(cfg_path),
        "artifacts/diagnostics/forensic_001/protocol.md": cfg["protocol_sha256"],
        "protocol.md": cfg["pilot_protocol_sha256"],
        "config/pilot.json": cfg["pilot_config_sha256"],
        "src/run_pilot.py": cfg["trusted_runner_sha256"],
        "src/pilot_estimators.py": cfg["trusted_estimators_sha256"],
        **cfg["source_data_sha256"],
        **cfg["pilot_inputs"],
    }
    for relative, wanted in expected.items():
        path = root / relative
        require(path.is_file(), f"missing frozen input: {relative}")
        require(sha256(path) == wanted, f"frozen input hash mismatch: {relative}")
    return expected


def load_contexts(root, cfg):
    test = pd.read_csv(root / "data/raw/expansion_data_test.csv")
    row_audit = pd.read_csv(
        root / "artifacts/experiment/pilot_001/full/predictions/row_audit.csv"
    )
    predictions = pd.read_parquet(
        root / "artifacts/experiment/pilot_001/full/predictions/predictions.parquet"
    )
    require(len(row_audit) == len(test), "row-audit length mismatch")
    require(
        np.array_equal(row_audit["row_index"].to_numpy(int), np.arange(len(test))),
        "row-audit order mismatch",
    )
    require(
        np.array_equal(
            row_audit["molecule_id"].astype(str).to_numpy(),
            test["Molecule Name"].astype(str).to_numpy(),
        ),
        "row-audit molecule alignment mismatch",
    )
    groups, scaffold_hashes = scaffold_contract(test)
    require(
        np.array_equal(
            row_audit["scaffold_sha256"].astype(str).to_numpy(), scaffold_hashes
        ),
        "row-audit scaffold hash mismatch",
    )
    selected = predictions[
        predictions["truth_endpoint"].isin(cfg["truth_endpoints"])
        & predictions["model"].isin(cfg["models"])
    ].copy()
    prediction_columns = {
        "truth_endpoint",
        "model",
        "row_index",
        "molecule_id",
        "prediction",
    }
    require(prediction_columns.issubset(selected.columns), "prediction schema mismatch")
    require(
        len(selected)
        == len(test) * len(cfg["truth_endpoints"]) * len(cfg["models"]),
        "prediction grid length mismatch",
    )
    require(
        not selected.duplicated(["truth_endpoint", "model", "row_index"]).any(),
        "duplicate prediction key",
    )
    require(
        np.isfinite(selected["prediction"].to_numpy(float)).all(),
        "nonfinite prediction",
    )
    contexts = {}
    for endpoint in cfg["truth_endpoints"]:
        numeric = pd.to_numeric(
            test[ENDPOINT_COLUMNS[endpoint]], errors="coerce"
        ).to_numpy(float)
        truth_idx = np.flatnonzero(np.isfinite(numeric))
        truth = transformed(numeric[truth_idx], endpoint)
        losses = {}
        for model in cfg["models"]:
            part = selected[
                (selected["truth_endpoint"] == endpoint)
                & (selected["model"] == model)
            ].sort_values("row_index")
            require(
                np.array_equal(part["row_index"].to_numpy(int), np.arange(len(test))),
                f"prediction row order mismatch: {endpoint}/{model}",
            )
            require(
                np.array_equal(
                    part["molecule_id"].astype(str).to_numpy(),
                    test["Molecule Name"].astype(str).to_numpy(),
                ),
                f"prediction molecule alignment mismatch: {endpoint}/{model}",
            )
            losses[model] = np.abs(
                truth - part["prediction"].to_numpy(float)[truth_idx]
            )
        endpoint_groups = groups[truth_idx]
        fold_id, fold_hash = fold_contract(endpoint_groups, cfg["crossfit_folds"])
        contexts[endpoint] = {
            "row_index": truth_idx,
            "molecule_id": test.iloc[truth_idx]["Molecule Name"].astype(str).to_numpy(),
            "groups": endpoint_groups,
            "scaffold_sha256": scaffold_hashes[truth_idx],
            "fold_id": fold_id,
            "fold_sha256": fold_hash,
            "losses": losses,
            "available": row_audit["available_HLM"].to_numpy(int)[truth_idx],
            "mar_probability": row_audit["propensity_clipped_HLM"].to_numpy(float)[truth_idx],
        }
    require(len(contexts["LogD"]["row_index"]) == 2270, "LogD cohort mismatch")
    require(len(contexts["KSOL"]["row_index"]) == 2170, "KSOL cohort mismatch")
    return contexts


def make_masks(contexts, cfg):
    masks = {}
    for endpoint, context in contexts.items():
        mcar_probability = float(context["available"].mean())
        mar_probability = context["mar_probability"]
        require(
            np.isfinite(mar_probability).all()
            and ((mar_probability >= 0.05) & (mar_probability <= 0.95)).all(),
            f"invalid frozen MAR probability: {endpoint}",
        )
        for seed in cfg["mask_seeds"]:
            seed = int(seed)
            mcar_p = np.full(len(context["row_index"]), mcar_probability)
            mcar_r = (
                np.random.default_rng(seed).random(len(mcar_p)) < mcar_p
            ).astype(int)
            mar_r = (
                np.random.default_rng(seed).random(len(mar_probability))
                < mar_probability
            ).astype(int)
            masks[(endpoint, "mcar_hlm", seed)] = (mcar_r, mcar_p)
            masks[(endpoint, "mar_hlm", seed)] = (mar_r, mar_probability.copy())
    return masks


def expected_scenario_keys(cfg):
    return {
        (endpoint, mechanism, int(seed), feature_set, model)
        for endpoint in cfg["truth_endpoints"]
        for mechanism in cfg["mechanisms"]
        for seed in cfg["mask_seeds"]
        for feature_set in cfg["feature_sets"]
        for model in cfg["models"]
    }


def verify_scenario_grid(scenarios, contexts, masks, cfg):
    keys = ["endpoint", "mechanism", "seed", "feature_set", "model"]
    required = set(
        keys
        + [
            "n_total",
            "n_observed",
            "observed_fraction",
            "generator_mean_probability",
            "fold_sha256",
            "mask_sha256",
            "oracle_probability_sha256",
            "full_truth",
            "cc",
            "or_mean",
            "fitted_ipw",
            "fitted_aipw",
            "fitted_aipw_se",
            "fitted_aipw_lcl",
            "fitted_aipw_ucl",
            "fitted_ess",
            "fitted_ess_fraction",
            "fitted_clip_fraction",
            "oracle_ipw",
            "oracle_aipw",
            "oracle_aipw_se",
            "oracle_aipw_lcl",
            "oracle_aipw_ucl",
            "oracle_ess",
            "oracle_ess_fraction",
        ]
    )
    require(required.issubset(scenarios.columns), "scenario schema mismatch")
    require(len(scenarios) == 640, "scenario row count mismatch")
    require(not scenarios.duplicated(keys).any(), "duplicate scenario key")
    observed = set(scenarios[keys].itertuples(index=False, name=None))
    require(observed == expected_scenario_keys(cfg), "scenario key grid mismatch")
    numeric = list(required - set(keys) - {
        "fold_sha256",
        "mask_sha256",
        "oracle_probability_sha256",
    })
    require(
        np.isfinite(scenarios[numeric].to_numpy(float)).all(),
        "nonfinite scenario value",
    )
    tolerance = float(cfg["reproduction_tolerance"])
    for column in [
        "observed_fraction",
        "fitted_ess_fraction",
        "fitted_clip_fraction",
        "oracle_ess_fraction",
    ]:
        require(
            (
                (scenarios[column] >= -tolerance)
                & (scenarios[column] <= 1 + tolerance)
            ).all(),
            f"scenario fraction outside [0, 1]: {column}",
        )
    shared_model_keys = ["endpoint", "mechanism", "seed", "feature_set"]
    for column in [
        "n_total",
        "n_observed",
        "observed_fraction",
        "generator_mean_probability",
        "fold_sha256",
        "mask_sha256",
        "oracle_probability_sha256",
        "fitted_ess",
        "fitted_ess_fraction",
        "fitted_clip_fraction",
        "oracle_ess",
        "oracle_ess_fraction",
    ]:
        require(
            scenarios.groupby(shared_model_keys)[column].nunique().max() == 1,
            f"scenario value is not shared across models: {column}",
        )
    for mask_key, (r, probability) in masks.items():
        endpoint, mechanism, seed = mask_key
        part = scenarios[
            (scenarios["endpoint"] == endpoint)
            & (scenarios["mechanism"] == mechanism)
            & (scenarios["seed"] == seed)
        ]
        require(
            len(part) == len(cfg["feature_sets"]) * len(cfg["models"]),
            f"scenario cell multiplicity mismatch: {mask_key}",
        )
        context = contexts[endpoint]
        expected_mask_hash = vector_hash(r, "u1")
        expected_probability_hash = vector_hash(probability, "<f8")
        for row in part.to_dict("records"):
            require(int(row["n_total"]) == len(r), "scenario n_total mismatch")
            require(int(row["n_observed"]) == int(r.sum()), "scenario n_observed mismatch")
            close(row["observed_fraction"], r.mean(), tolerance, "observed_fraction")
            close(
                row["generator_mean_probability"],
                probability.mean(),
                tolerance,
                "generator_mean_probability",
            )
            require(row["mask_sha256"] == expected_mask_hash, "mask hash mismatch")
            require(
                row["oracle_probability_sha256"] == expected_probability_hash,
                "oracle probability hash mismatch",
            )
            require(row["fold_sha256"] == context["fold_sha256"], "fold hash mismatch")
            close(
                row["full_truth"],
                context["losses"][row["model"]].mean(),
                tolerance,
                "full_truth",
            )
            close(
                row["fitted_aipw_lcl"],
                row["fitted_aipw"] - 1.96 * row["fitted_aipw_se"],
                tolerance,
                "fitted interval lower",
            )
            close(
                row["fitted_aipw_ucl"],
                row["fitted_aipw"] + 1.96 * row["fitted_aipw_se"],
                tolerance,
                "fitted interval upper",
            )
            close(
                row["oracle_aipw_lcl"],
                row["oracle_aipw"] - 1.96 * row["oracle_aipw_se"],
                tolerance,
                "oracle interval lower",
            )
            close(
                row["oracle_aipw_ucl"],
                row["oracle_aipw"] + 1.96 * row["oracle_aipw_se"],
                tolerance,
                "oracle interval upper",
            )
    return {
        "one_fold_vector_per_endpoint": bool(
            (scenarios.groupby("endpoint")["fold_sha256"].nunique() == 1).all()
        ),
        "one_mask_per_endpoint_mechanism_seed": bool(
            (
                scenarios.groupby(["endpoint", "mechanism", "seed"])["mask_sha256"]
                .nunique()
                == 1
            ).all()
        ),
        "one_oracle_probability_per_endpoint_mechanism_seed": bool(
            (
                scenarios.groupby(["endpoint", "mechanism", "seed"])[
                    "oracle_probability_sha256"
                ].nunique()
                == 1
            ).all()
        ),
    }


def f1_reproduction(root, scenarios, cfg):
    pilot = pd.read_csv(
        root / "artifacts/experiment/pilot_001/full/scenario_results.csv"
    )
    mcar = pilot[
        (pilot["mechanism"] == "mcar")
        & (pilot["mask"] == "HLM")
        & pilot["anchor"].isna()
    ].copy()
    mar = pilot[
        (pilot["mechanism"] == "mar")
        & pilot["mask"].isna()
        & (pilot["anchor"] == "HLM")
    ].copy()
    mcar["diagnostic_mechanism"] = "mcar_hlm"
    mar["diagnostic_mechanism"] = "mar_hlm"
    source = pd.concat([mcar, mar], ignore_index=True).rename(
        columns={"truth_endpoint": "endpoint"}
    )
    source["seed"] = source["seed"].astype(int)
    diagnostic = scenarios[scenarios["feature_set"] == "full_2267"].copy()
    diagnostic = diagnostic.rename(columns={"mechanism": "diagnostic_mechanism"})
    keys = ["endpoint", "diagnostic_mechanism", "seed", "model"]
    expected = {
        (endpoint, mechanism, int(seed), model)
        for endpoint in cfg["truth_endpoints"]
        for mechanism in cfg["mechanisms"]
        for seed in cfg["mask_seeds"]
        for model in cfg["models"]
    }
    source_keys = set(source[keys].itertuples(index=False, name=None))
    diagnostic_keys = set(diagnostic[keys].itertuples(index=False, name=None))
    grid_ok = bool(
        len(mcar) == cfg["f1_contract"]["mcar_source"]["expected_rows"]
        and len(mar) == cfg["f1_contract"]["mar_source"]["expected_rows"]
        and len(source) == len(diagnostic) == cfg["f1_contract"]["expected_rows"]
        and not source.duplicated(keys).any()
        and not diagnostic.duplicated(keys).any()
        and source_keys == diagnostic_keys == expected
    )
    if not grid_ok:
        return {"pass": False, "grid_ok": False, "reason": "F1 key grid mismatch"}
    mapping = {
        "n_total": "n_total",
        "n_observed": "n_observed",
        "generator_mean_probability": "generator_mean_probability",
        "full_truth": "full_truth",
        "cc": "cc",
        "ipw": "fitted_ipw",
        "aipw": "fitted_aipw",
        "se": "fitted_aipw_se",
        "ess_fraction": "fitted_ess_fraction",
        "clip_fraction": "fitted_clip_fraction",
    }
    source_columns = source[keys + list(mapping)].rename(
        columns={column: f"pilot_{column}" for column in mapping}
    )
    diagnostic_columns = diagnostic[keys + list(mapping.values())].rename(
        columns={column: f"diagnostic_{column}" for column in mapping.values()}
    )
    merged = source_columns.merge(
        diagnostic_columns, on=keys, validate="one_to_one"
    )
    maximum = {}
    for pilot_column, diagnostic_column in mapping.items():
        difference = np.abs(
            merged[f"pilot_{pilot_column}"].to_numpy(float)
            - merged[f"diagnostic_{diagnostic_column}"].to_numpy(float)
        )
        require(np.isfinite(difference).all(), f"nonfinite F1 difference: {pilot_column}")
        maximum[pilot_column] = float(difference.max())
    tolerance = float(cfg["reproduction_tolerance"])
    return {
        "pass": bool(all(value <= tolerance for value in maximum.values())),
        "grid_ok": True,
        "rows": int(len(merged)),
        "tolerance": tolerance,
        "maximum_absolute_difference": maximum,
    }


def f2_replay(replay, scenarios, contexts, masks, cfg):
    contract = cfg["f2_contract"]
    key = contract["unique_key"]
    cells = ["endpoint", "mechanism", "seed", "feature_set", "model"]
    required = set(key + contract["required_row_columns"])
    require(required.issubset(replay.columns), "F2 row schema mismatch")
    expected_cells = {
        (endpoint, mechanism, int(cfg["row_replay_seed"]), feature_set, model)
        for endpoint in cfg["truth_endpoints"]
        for mechanism in cfg["mechanisms"]
        for feature_set in cfg["feature_sets"]
        for model in cfg["models"]
    }
    observed_cells = set(replay[cells].drop_duplicates().itertuples(index=False, name=None))
    grid_ok = bool(
        len(replay) == contract["expected_rows"]
        and not replay.duplicated(key).any()
        and len(observed_cells) == contract["expected_cells"]
        and observed_cells == expected_cells
    )
    if not grid_ok:
        return {"pass": False, "grid_ok": False, "reason": "F2 row grid mismatch"}
    numeric = [
        "seed",
        "row_index",
        "fold_id",
        "R",
        "loss",
        "mhat",
        "p_fitted",
        "p_oracle",
        "augmentation_fitted",
        "pseudo_outcome_fitted",
        "influence_fitted",
        "augmentation_oracle",
        "pseudo_outcome_oracle",
        "influence_oracle",
    ]
    require(np.isfinite(replay[numeric].to_numpy(float)).all(), "nonfinite F2 value")
    require(np.isin(replay["R"].to_numpy(float), [0.0, 1.0]).all(), "nonbinary F2 R")
    lo, hi = map(float, cfg["propensity_clip"])
    require(
        ((replay["p_fitted"] >= lo) & (replay["p_fitted"] <= hi)).all(),
        "fitted propensity outside frozen clipping interval",
    )
    require(
        ((replay["p_oracle"] > 0) & (replay["p_oracle"] <= 1)).all(),
        "invalid oracle propensity",
    )
    tolerance = float(cfg["reproduction_tolerance"])
    maximum = {
        "augmentation_fitted": 0.0,
        "pseudo_outcome_fitted": 0.0,
        "influence_fitted": 0.0,
        "augmentation_oracle": 0.0,
        "pseudo_outcome_oracle": 0.0,
        "influence_oracle": 0.0,
        "fitted_ipw": 0.0,
        "fitted_aipw": 0.0,
        "fitted_aipw_se": 0.0,
        "oracle_ipw": 0.0,
        "oracle_aipw": 0.0,
        "oracle_aipw_se": 0.0,
    }
    scenario_index = scenarios.set_index(cells)
    vectors_ok = True
    frozen_generator_ok = True
    scenario_hashes_ok = True
    for cell, part in replay.groupby(cells, sort=False):
        part = part.sort_values("row_index")
        endpoint, mechanism, seed, _, model = cell
        context = contexts[endpoint]
        expected_r, expected_p = masks[(endpoint, mechanism, int(seed))]
        vectors_ok &= bool(
            np.array_equal(part["row_index"].to_numpy(int), context["row_index"])
            and np.array_equal(part["molecule_id"].astype(str).to_numpy(), context["molecule_id"])
            and np.array_equal(part["scaffold_group"].astype(str).to_numpy(), context["groups"].astype(str))
            and np.array_equal(part["scaffold_sha256"].astype(str).to_numpy(), context["scaffold_sha256"])
            and np.array_equal(part["fold_id"].to_numpy(int), context["fold_id"])
        )
        frozen_generator_ok &= bool(
            np.array_equal(part["R"].to_numpy(int), expected_r)
            and np.array_equal(part["p_oracle"].to_numpy(float), expected_p)
        )
        require(
            np.allclose(
                part["loss"].to_numpy(float),
                context["losses"][model],
                rtol=0,
                atol=tolerance,
            ),
            f"F2 loss vector mismatch: {cell}",
        )
        stored = scenario_index.loc[cell]
        scenario_hashes_ok &= bool(
            stored["mask_sha256"] == vector_hash(expected_r, "u1")
            and stored["oracle_probability_sha256"] == vector_hash(expected_p, "<f8")
            and stored["fold_sha256"] == context["fold_sha256"]
        )
        r = part["R"].to_numpy(float)
        loss = part["loss"].to_numpy(float)
        mhat = part["mhat"].to_numpy(float)
        groups = context["groups"]
        observed = r.astype(bool)
        close(stored["cc"], loss[observed].mean(), tolerance, "F2 complete case")
        close(stored["or_mean"], mhat.mean(), tolerance, "F2 outcome regression")
        for suffix, probability_column in [("fitted", "p_fitted"), ("oracle", "p_oracle")]:
            probability = part[probability_column].to_numpy(float)
            augmentation = r * (loss - mhat) / probability
            pseudo = mhat + augmentation
            influence = pseudo - pseudo.mean()
            weights = r / probability
            require(weights.sum() > 0, "F2 empty weighted sample")
            ipw = float(np.sum(weights * loss) / weights.sum())
            aipw = float(pseudo.mean())
            se = cluster_se(influence, groups)
            require(
                np.isfinite(
                    np.concatenate([augmentation, pseudo, influence, [ipw, aipw, se]])
                ).all(),
                "nonfinite recomputed F2 value",
            )
            for name, expected_values in [
                (f"augmentation_{suffix}", augmentation),
                (f"pseudo_outcome_{suffix}", pseudo),
                (f"influence_{suffix}", influence),
            ]:
                difference = np.abs(
                    expected_values - part[name].to_numpy(float)
                )
                maximum[name] = max(maximum[name], float(difference.max()))
            maximum[f"{suffix}_ipw"] = max(
                maximum[f"{suffix}_ipw"], abs(ipw - float(stored[f"{suffix}_ipw"]))
            )
            maximum[f"{suffix}_aipw"] = max(
                maximum[f"{suffix}_aipw"], abs(aipw - float(stored[f"{suffix}_aipw"]))
            )
            maximum[f"{suffix}_aipw_se"] = max(
                maximum[f"{suffix}_aipw_se"], abs(se - float(stored[f"{suffix}_aipw_se"]))
            )
            ess = float(weights.sum() ** 2 / np.square(weights).sum())
            close(stored[f"{suffix}_ess"], ess, tolerance, f"F2 {suffix} ESS")
            close(
                stored[f"{suffix}_ess_fraction"],
                ess / observed.sum(),
                tolerance,
                f"F2 {suffix} ESS fraction",
            )
    shared_ok = bool(
        replay.groupby(["endpoint", "mechanism", "row_index"])["R"].nunique().max() == 1
        and replay.groupby(["endpoint", "mechanism", "row_index"])["p_oracle"].nunique().max() == 1
        and replay.groupby(["endpoint", "row_index"])["fold_id"].nunique().max() == 1
        and replay.groupby(["endpoint", "row_index"])["scaffold_sha256"].nunique().max() == 1
        and replay.groupby(["endpoint", "mechanism", "feature_set", "row_index"])["p_fitted"].nunique().max() == 1
    )
    hash_ok = all(
        vector_hash(contexts[endpoint]["fold_id"], "<i8")
        == contexts[endpoint]["fold_sha256"]
        for endpoint in cfg["truth_endpoints"]
    )
    return {
        "pass": bool(
            vectors_ok
            and shared_ok
            and hash_ok
            and frozen_generator_ok
            and scenario_hashes_ok
            and all(value <= tolerance for value in maximum.values())
        ),
        "grid_ok": True,
        "rows": int(len(replay)),
        "cells": int(len(observed_cells)),
        "tolerance": tolerance,
        "frozen_vectors_match": bool(vectors_ok),
        "shared_mask_probability_fold_keys": shared_ok,
        "fold_hashes_match": hash_ok,
        "frozen_generator_vectors_match": bool(frozen_generator_ok),
        "scenario_mask_oracle_fold_hashes_match": bool(scenario_hashes_ok),
        "maximum_absolute_difference": maximum,
    }


def safe_ratio(numerator, denominator):
    return None if denominator == 0 else float(numerator / denominator)


def summarize(scenarios, f1, f2, contexts, shared_contracts, cfg):
    columns = [
        "cc",
        "or_mean",
        "fitted_ipw",
        "fitted_aipw",
        "oracle_ipw",
        "oracle_aipw",
    ]
    pooled = []
    for (mechanism, feature_set), part in scenarios.groupby(
        ["mechanism", "feature_set"], sort=False
    ):
        require(len(part) == 160, "pooled summary row count mismatch")
        pooled.append(
            {
                "mechanism": mechanism,
                "feature_set": feature_set,
                "rows": int(len(part)),
                "median_absolute_bias": {
                    column: float(np.median(np.abs(part[column] - part["full_truth"])))
                    for column in columns
                },
                "median_fitted_clip_fraction": float(part["fitted_clip_fraction"].median()),
                "median_fitted_ess_fraction": float(part["fitted_ess_fraction"].median()),
                "median_oracle_ess_fraction": float(part["oracle_ess_fraction"].median()),
            }
        )
    lookup = {(row["mechanism"], row["feature_set"]): row for row in pooled}
    thresholds = cfg["signature_thresholds"]
    propensity_details = []
    augmentation_details = []
    dimension_details = []
    for mechanism in cfg["mechanisms"]:
        full = lookup[(mechanism, "full_2267")]
        low = lookup[(mechanism, "lowdim_12")]
        full_bias = full["median_absolute_bias"]
        low_bias = low["median_absolute_bias"]
        clip_reduction = (
            full["median_fitted_clip_fraction"] - low["median_fitted_clip_fraction"]
        )
        propensity_details.append(
            {
                "mechanism": mechanism,
                "ratio": safe_ratio(full_bias["oracle_aipw"], full_bias["fitted_aipw"]),
                "meets": full_bias["oracle_aipw"]
                <= thresholds["oracle_aipw_fraction_of_fitted_aipw"]
                * full_bias["fitted_aipw"],
            }
        )
        augmentation_details.append(
            {
                "mechanism": mechanism,
                "ratio": safe_ratio(full_bias["oracle_aipw"], full_bias["oracle_ipw"]),
                "meets": full_bias["oracle_aipw"]
                > thresholds["oracle_aipw_multiple_of_oracle_ipw"]
                * full_bias["oracle_ipw"],
            }
        )
        dimension_details.append(
            {
                "mechanism": mechanism,
                "bias_ratio": safe_ratio(low_bias["fitted_aipw"], full_bias["fitted_aipw"]),
                "clip_fraction_reduction": clip_reduction,
                "meets": low_bias["fitted_aipw"]
                <= thresholds["lowdim_aipw_fraction_of_full"]
                * full_bias["fitted_aipw"]
                and clip_reduction >= thresholds["minimum_clip_fraction_reduction"],
            }
        )
    signatures = {
        "PROPENSITY_INSTABILITY_SUPPORTED": {
            "supported": bool(all(row["meets"] for row in propensity_details)),
            "details": propensity_details,
        },
        "AUGMENTATION_INSTABILITY_SUPPORTED": {
            "supported": bool(all(row["meets"] for row in augmentation_details)),
            "details": augmentation_details,
        },
        "HIGH_DIMENSION_INSTABILITY_SUPPORTED": {
            "supported": bool(all(row["meets"] for row in dimension_details)),
            "details": dimension_details,
        },
    }
    if not f1["pass"] or not f2["pass"]:
        decision = "IMPLEMENTATION_MISMATCH"
    elif any(value["supported"] for value in signatures.values()):
        decision = "MECHANISM_SIGNAL_DETECTED"
    else:
        decision = "MECHANISM_INDETERMINATE"
    return {
        "diagnostic_id": cfg["diagnostic_id"],
        "pilot_decision_unchanged": "STOP_PILOT_NEGATIVE",
        "decision": decision,
        "fold_id_sha256": {
            endpoint: contexts[endpoint]["fold_sha256"]
            for endpoint in cfg["truth_endpoints"]
        },
        "fold_id_hash_encoding": "ordered little-endian signed int64 raw bytes",
        "shared_contracts": shared_contracts,
        "F1_reproduction": f1,
        "F2_row_replay": f2,
        "pooled_mechanism_summaries": pooled,
        "signature_thresholds": thresholds,
        "signatures": signatures,
    }


def verify(output):
    root = Path(__file__).resolve().parents[1]
    output = Path(output).resolve()
    cfg_path = root / "artifacts/diagnostics/forensic_001/config.json"
    require(sha256(cfg_path) == EXPECTED_CONFIG_SHA256, "forensic config hash mismatch")
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    require(
        cfg["protocol_amendment"] == 1
        and cfg["protocol_amended_before_results"] is True,
        "forensic amendment state mismatch",
    )
    require(rdkit.__version__ == cfg["rdkit_version"], "RDKit version mismatch")
    inputs = expected_input_hashes(root, cfg, cfg_path)
    require(
        sha256(root / "src/run_forensic_diagnostic.py")
        == TRUSTED_FORENSIC_RUNNER_SHA256,
        "forensic runner differs from audited trust anchor",
    )
    required_outputs = {
        "scenario_decomposition.csv",
        "row_replay_seed31001.parquet",
        "summary.json",
        "verification.md",
    }
    require(output.is_dir(), f"forensic output directory is missing: {output}")
    manifest_path = output / "run_manifest.json"
    require(manifest_path.is_file(), "run manifest is missing")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    require(
        set(manifest["output_hashes_before_manifest"]) == required_outputs,
        "manifest output grid mismatch",
    )
    actual_outputs = {
        str(path.relative_to(output)).replace("\\", "/")
        for path in output.rglob("*")
        if path.is_file() and path != manifest_path
    }
    require(actual_outputs == required_outputs, "unmanifested forensic output file")
    for relative, wanted in manifest["output_hashes_before_manifest"].items():
        require(sha256(output / relative) == wanted, f"output hash mismatch: {relative}")
    require(manifest["input_hashes"] == inputs, "manifest input hashes mismatch")
    require(
        manifest["runner_sha256"] == TRUSTED_FORENSIC_RUNNER_SHA256,
        "manifest runner hash mismatch",
    )
    require(
        manifest["diagnostic_id"] == cfg["diagnostic_id"]
        and manifest["mode"] == "full"
        and manifest["protocol_amendment"] == 1,
        "manifest identity mismatch",
    )
    scenarios = pd.read_csv(output / "scenario_decomposition.csv")
    replay = pd.read_parquet(output / "row_replay_seed31001.parquet")
    contexts = load_contexts(root, cfg)
    masks = make_masks(contexts, cfg)
    shared_contracts = verify_scenario_grid(scenarios, contexts, masks, cfg)
    require(all(shared_contracts.values()), "shared scenario contract mismatch")
    f1 = f1_reproduction(root, scenarios, cfg)
    f2 = f2_replay(replay, scenarios, contexts, masks, cfg)
    recomputed = summarize(scenarios, f1, f2, contexts, shared_contracts, cfg)
    tolerance = float(cfg["reproduction_tolerance"])
    compare_tree(recomputed, summary, tolerance, "summary")
    require(
        manifest["result_rows"]
        == {
            "scenario_decomposition.csv": 640,
            "row_replay_seed31001.parquet": 35520,
        },
        "manifest result row counts mismatch",
    )
    require(manifest["decision"] == recomputed["decision"], "manifest decision mismatch")
    require(
        manifest["fold_id_sha256"] == recomputed["fold_id_sha256"]
        and manifest["fold_id_hash_encoding"]
        == "ordered little-endian signed int64 raw bytes",
        "manifest fold contract mismatch",
    )
    require(
        manifest["shared_contracts"] == shared_contracts,
        "manifest shared-contract mismatch",
    )
    expected_shapes = {
        endpoint: {
            "full_2267": {
                "shape": [len(contexts[endpoint]["row_index"]), 2267],
                "dtype": "float32",
            },
            "lowdim_12": {
                "shape": [len(contexts[endpoint]["row_index"]), 12],
                "dtype": "float32",
            },
        }
        for endpoint in cfg["truth_endpoints"]
    }
    require(manifest["feature_sets"] == expected_shapes, "manifest feature-shape mismatch")
    require(manifest["versions"]["rdkit"] == cfg["rdkit_version"], "manifest RDKit mismatch")
    return {
        "status": "PASS",
        "status_scope": "artifact replay agrees with the stored diagnostic decision",
        "adjudicative_integrity_pass": bool(f1["pass"] and f2["pass"]),
        "verification_scope": "frozen input and output hashes; independent source-row, prediction-loss, scaffold, fold, scenario-grid, all-seed mask, seed-31001 row-formula, F1/F2, signature, decision, and manifest replay; excludes an independent refit of the nuisance models",
        "forensic_runner_sha256": TRUSTED_FORENSIC_RUNNER_SHA256,
        "manifest_sha256": sha256(manifest_path),
        "scenario_rows": len(scenarios),
        "replay_rows": len(replay),
        "F1_pass": bool(f1["pass"]),
        "F2_pass": bool(f2["pass"]),
        "decision": recomputed["decision"],
        "maximum_F1_difference": max(f1.get("maximum_absolute_difference", {"x": np.inf}).values()),
        "maximum_F2_difference": max(f2.get("maximum_absolute_difference", {"x": np.inf}).values()),
        "output_hashes_verified": len(required_outputs),
        "input_hashes_verified": len(inputs),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "artifacts/diagnostics/forensic_001/full",
    )
    args = parser.parse_args()
    print(json.dumps(verify(args.output), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
