import argparse
import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import rdkit
import sklearn
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

from pilot_estimators import crossfit_estimators
from run_pilot import ENDPOINT_COLUMNS, PRIMARY_MODELS, file_hash, prepare_study, transform_y


EXPECTED_CONFIG_SHA256 = "72dfed2fd134fff7b5e6ccc6d6f0300659def6bd0d6d6f510aec4372a929c203"


def write_json(path, value):
    Path(path).write_text(
        json.dumps(value, indent=2, ensure_ascii=False), encoding="utf-8"
    )


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def vector_hash(values, dtype):
    return hashlib.sha256(np.asarray(values, dtype=dtype).tobytes()).hexdigest()


def verify_inputs(root, cfg_path):
    require(file_hash(cfg_path) == EXPECTED_CONFIG_SHA256, "forensic config hash mismatch")
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    require(
        cfg["protocol_amendment"] == 1 and cfg["protocol_amended_before_results"] is True,
        "unexpected forensic protocol amendment state",
    )
    expected = {
        "artifacts/diagnostics/forensic_001/protocol.md": cfg["protocol_sha256"],
        "protocol.md": cfg["pilot_protocol_sha256"],
        "config/pilot.json": cfg["pilot_config_sha256"],
        "src/run_pilot.py": cfg["trusted_runner_sha256"],
        "src/pilot_estimators.py": cfg["trusted_estimators_sha256"],
        **cfg["source_data_sha256"],
        **cfg["pilot_inputs"],
    }
    observed = {
        "artifacts/diagnostics/forensic_001/config.json": file_hash(cfg_path)
    }
    for relative, wanted in expected.items():
        path = root / relative
        require(path.is_file(), f"missing frozen input: {relative}")
        got = file_hash(path)
        require(got == wanted, f"hash mismatch for {relative}: {got}")
        observed[relative] = got
    require(rdkit.__version__ == cfg["rdkit_version"], "RDKit version mismatch")
    pilot_cfg = json.loads((root / "config" / "pilot.json").read_text(encoding="utf-8"))
    require(pilot_cfg["protocol_sha256"] == cfg["pilot_protocol_sha256"], "pilot protocol mismatch")
    require(pilot_cfg["truth_endpoints"] == cfg["truth_endpoints"], "truth endpoint mismatch")
    require(pilot_cfg["primary_models"] == cfg["models"] == PRIMARY_MODELS, "model order mismatch")
    require(cfg["anchor"] == "HLM", "unexpected forensic anchor")
    require(pilot_cfg["mask_seeds"] == cfg["mask_seeds"], "mask seed mismatch")
    require(pilot_cfg["crossfit_folds"] == cfg["crossfit_folds"], "fold count mismatch")
    require(pilot_cfg["propensity_model"] == cfg["propensity_model"], "propensity contract mismatch")
    require(pilot_cfg["propensity_clip"] == cfg["propensity_clip"], "propensity clip mismatch")
    for key in ["alpha", "solver", "tol", "truncate_predictions"]:
        require(pilot_cfg["loss_nuisance"][key] == cfg["loss_nuisance"][key], f"loss nuisance mismatch: {key}")
    pilot_data_hashes = {
        f"data/raw/{name}": value for name, value in pilot_cfg["dataset"]["files"].items()
    }
    require(pilot_data_hashes == cfg["source_data_sha256"], "pilot/source hash contract mismatch")
    return cfg, pilot_cfg, observed


def fold_contract(groups, folds):
    split = list(GroupKFold(n_splits=int(folds)).split(np.zeros(len(groups)), groups=groups))
    fold_id = np.full(len(groups), -1, dtype=np.int64)
    for value, (_, valid) in enumerate(split):
        fold_id[valid] = value
    require(np.all(fold_id >= 0), "incomplete fold assignment")
    return split, fold_id, vector_hash(fold_id, "<i8")


def scaled_folds(X, split):
    prepared = []
    for train, valid in split:
        scaler = StandardScaler()
        prepared.append(
            (train, valid, scaler.fit_transform(X[train]), scaler.transform(X[valid]))
        )
    return prepared


def cluster_se(influence, groups):
    codes = pd.factorize(np.asarray(groups, dtype=object), sort=False)[0]
    sums = np.bincount(codes, weights=np.asarray(influence, dtype=float))
    count = len(sums)
    require(count > 1, "cluster SE needs at least two scaffolds")
    return float(
        np.sqrt((count / (count - 1)) * np.square(sums).sum()) / len(influence)
    )


def estimates_for_probability(r, loss, mhat, probability, groups):
    r = np.asarray(r, dtype=float)
    p = np.asarray(probability, dtype=float)
    weights = r / p
    denominator = weights.sum()
    require(denominator > 0, "mask has no observed rows")
    ipw = float(np.sum(weights * loss) / denominator)
    augmentation = r * (loss - mhat) / p
    pseudo = mhat + augmentation
    aipw = float(pseudo.mean())
    influence = pseudo - aipw
    se = cluster_se(influence, groups)
    ess = float(denominator**2 / np.square(weights).sum())
    observed = int(r.sum())
    return {
        "ipw": ipw,
        "aipw": aipw,
        "se": se,
        "lcl": aipw - 1.96 * se,
        "ucl": aipw + 1.96 * se,
        "ess": ess,
        "ess_fraction": ess / observed,
        "augmentation": augmentation,
        "pseudo": pseudo,
        "influence": influence,
    }


def load_contexts(root, cfg, pilot_cfg):
    study = prepare_study(root, pilot_cfg)
    audit = pd.read_csv(root / "artifacts/experiment/pilot_001/full/predictions/row_audit.csv")
    predictions = pd.read_parquet(
        root / "artifacts/experiment/pilot_001/full/predictions/predictions.parquet"
    )
    n_test = len(study["test"])
    require(len(audit) == n_test, "row audit length mismatch")
    require(np.array_equal(audit["row_index"].to_numpy(), np.arange(n_test)), "row audit order mismatch")
    require(
        audit["molecule_id"].astype(str).tolist()
        == study["test"]["Molecule Name"].astype(str).tolist(),
        "row audit molecule alignment mismatch",
    )
    expected_scaffold_hash = np.asarray(
        [hashlib.sha256(str(x).encode("utf-8")).hexdigest() for x in study["groups"]]
    )
    require(
        np.array_equal(expected_scaffold_hash, audit["scaffold_sha256"].astype(str).to_numpy()),
        "row audit scaffold mismatch",
    )
    require(
        np.array_equal(
            audit["available_HLM"].to_numpy(int), study["availability"]["HLM"]
        ),
        "row audit HLM availability mismatch",
    )
    require(
        np.allclose(audit["max_train_similarity"], study["similarity"], rtol=0, atol=1e-15)
        and np.allclose(
            audit["normalized_identifier_order"], study["identifier_order"], rtol=0, atol=1e-15
        ),
        "row audit covariate mismatch",
    )
    selected = predictions[
        predictions["truth_endpoint"].isin(cfg["truth_endpoints"])
        & predictions["model"].isin(cfg["models"])
    ].copy()
    require(
        len(selected) == n_test * len(cfg["truth_endpoints"]) * len(cfg["models"]),
        "prediction grid length mismatch",
    )
    require(
        not selected.duplicated(["truth_endpoint", "model", "row_index"]).any(),
        "duplicate prediction key",
    )
    require(np.isfinite(selected["prediction"].to_numpy(float)).all(), "nonfinite prediction")
    descriptor_index = {name: i for i, name in enumerate(study["descriptor_names"])}
    low_names = cfg["feature_sets"]["lowdim_12"]["descriptors"]
    require(all(name in descriptor_index for name in low_names), "missing low-dimensional descriptor")
    lowdim = np.hstack(
        [
            study["test_features"]["desc"][:, [descriptor_index[name] for name in low_names]],
            study["similarity"][:, None],
            study["identifier_order"][:, None],
        ]
    ).astype(np.float32)
    require(lowdim.shape == (n_test, 12) and np.isfinite(lowdim).all(), "invalid lowdim_12")
    contexts = {}
    for endpoint in cfg["truth_endpoints"]:
        numeric = pd.to_numeric(
            study["test"][ENDPOINT_COLUMNS[endpoint]], errors="coerce"
        ).to_numpy(float)
        truth_idx = np.flatnonzero(np.isfinite(numeric))
        y = transform_y(numeric[truth_idx], endpoint)
        loss_columns = []
        for model in cfg["models"]:
            part = selected[
                (selected["truth_endpoint"] == endpoint) & (selected["model"] == model)
            ].sort_values("row_index")
            require(
                np.array_equal(part["row_index"].to_numpy(int), np.arange(n_test)),
                f"prediction order mismatch: {endpoint}/{model}",
            )
            require(
                part["molecule_id"].astype(str).tolist()
                == study["test"]["Molecule Name"].astype(str).tolist(),
                f"prediction molecule mismatch: {endpoint}/{model}",
            )
            loss_columns.append(
                np.abs(y - part["prediction"].to_numpy(float)[truth_idx])
            )
        groups = study["groups"][truth_idx]
        split, fold_id, fold_hash = fold_contract(groups, cfg["crossfit_folds"])
        contexts[endpoint] = {
            "truth_idx": truth_idx,
            "row_index": truth_idx,
            "molecule_id": study["test"].iloc[truth_idx]["Molecule Name"].astype(str).to_numpy(),
            "groups": groups,
            "scaffold_sha256": expected_scaffold_hash[truth_idx],
            "fold_id": fold_id,
            "fold_sha256": fold_hash,
            "split": split,
            "losses": np.column_stack(loss_columns),
            "features": {
                "full_2267": study["x_obs"][truth_idx],
                "lowdim_12": lowdim[truth_idx],
            },
            "available_hlm": audit["available_HLM"].to_numpy(int)[truth_idx],
            "mar_probability": audit["propensity_clipped_HLM"].to_numpy(float)[truth_idx],
        }
    require(len(contexts["LogD"]["truth_idx"]) == 2270, "unexpected LogD truth cohort")
    require(len(contexts["KSOL"]["truth_idx"]) == 2170, "unexpected KSOL truth cohort")
    return study, contexts


def make_masks(context, cfg):
    masks = {}
    mcar_probability = float(context["available_hlm"].mean())
    for seed in cfg["mask_seeds"]:
        random_values = np.random.default_rng(int(seed)).random(len(context["truth_idx"]))
        masks[("mcar_hlm", int(seed))] = (
            (random_values < mcar_probability).astype(int),
            np.full(len(random_values), mcar_probability, dtype=float),
        )
        random_values = np.random.default_rng(int(seed)).random(len(context["truth_idx"]))
        masks[("mar_hlm", int(seed))] = (
            (random_values < context["mar_probability"]).astype(int),
            context["mar_probability"].copy(),
        )
    return masks


def replay_rows(context, endpoint, mechanism, seed, feature_set, models, result, r, oracle_p):
    frames = []
    fitted_p = result["p"]
    for j, model in enumerate(models):
        fitted = estimates_for_probability(
            r, context["losses"][:, j], result["mhat"][:, j], fitted_p, context["groups"]
        )
        oracle = estimates_for_probability(
            r, context["losses"][:, j], result["mhat"][:, j], oracle_p, context["groups"]
        )
        frames.append(
            pd.DataFrame(
                {
                    "endpoint": endpoint,
                    "mechanism": mechanism,
                    "seed": int(seed),
                    "feature_set": feature_set,
                    "model": model,
                    "row_index": context["row_index"],
                    "molecule_id": context["molecule_id"],
                    "scaffold_group": context["groups"],
                    "scaffold_sha256": context["scaffold_sha256"],
                    "fold_id": context["fold_id"],
                    "R": r,
                    "loss": context["losses"][:, j],
                    "mhat": result["mhat"][:, j],
                    "p_fitted": fitted_p,
                    "p_oracle": oracle_p,
                    "augmentation_fitted": fitted["augmentation"],
                    "pseudo_outcome_fitted": fitted["pseudo"],
                    "influence_fitted": fitted["influence"],
                    "augmentation_oracle": oracle["augmentation"],
                    "pseudo_outcome_oracle": oracle["pseudo"],
                    "influence_oracle": oracle["influence"],
                }
            )
        )
    return frames


def run_grid(contexts, cfg):
    rows = []
    replay = []
    fold_hashes = {}
    for endpoint in cfg["truth_endpoints"]:
        context = contexts[endpoint]
        fold_hashes[endpoint] = context["fold_sha256"]
        masks = make_masks(context, cfg)
        for feature_set in cfg["feature_sets"]:
            X = context["features"][feature_set]
            prepared = scaled_folds(X, context["split"])
            for mechanism in cfg["mechanisms"]:
                for seed in cfg["mask_seeds"]:
                    r, oracle_p = masks[(mechanism, int(seed))]
                    result = crossfit_estimators(
                        X,
                        context["groups"],
                        r,
                        context["losses"],
                        int(cfg["crossfit_folds"]),
                        cfg["propensity_model"],
                        cfg["loss_nuisance"],
                        cfg["propensity_clip"],
                        prepared=prepared,
                    )
                    mask_sha256 = vector_hash(r, "u1")
                    oracle_sha256 = vector_hash(oracle_p, "<f8")
                    for j, model in enumerate(cfg["models"]):
                        fitted = result["estimates"][j]
                        oracle = estimates_for_probability(
                            r,
                            context["losses"][:, j],
                            result["mhat"][:, j],
                            oracle_p,
                            context["groups"],
                        )
                        truth = float(context["losses"][:, j].mean())
                        rows.append(
                            {
                                "endpoint": endpoint,
                                "mechanism": mechanism,
                                "seed": int(seed),
                                "feature_set": feature_set,
                                "model": model,
                                "n_total": int(len(r)),
                                "n_observed": int(r.sum()),
                                "observed_fraction": float(r.mean()),
                                "generator_mean_probability": float(oracle_p.mean()),
                                "fold_sha256": context["fold_sha256"],
                                "mask_sha256": mask_sha256,
                                "oracle_probability_sha256": oracle_sha256,
                                "full_truth": truth,
                                "cc": fitted["cc"],
                                "or_mean": float(result["mhat"][:, j].mean()),
                                "fitted_ipw": fitted["ipw"],
                                "fitted_aipw": fitted["aipw"],
                                "fitted_aipw_se": fitted["se"],
                                "fitted_aipw_lcl": fitted["lcl"],
                                "fitted_aipw_ucl": fitted["ucl"],
                                "fitted_ess": fitted["ess"],
                                "fitted_ess_fraction": fitted["ess_fraction"],
                                "fitted_clip_fraction": fitted["clip_fraction"],
                                "oracle_ipw": oracle["ipw"],
                                "oracle_aipw": oracle["aipw"],
                                "oracle_aipw_se": oracle["se"],
                                "oracle_aipw_lcl": oracle["lcl"],
                                "oracle_aipw_ucl": oracle["ucl"],
                                "oracle_ess": oracle["ess"],
                                "oracle_ess_fraction": oracle["ess_fraction"],
                            }
                        )
                    if int(seed) == int(cfg["row_replay_seed"]):
                        replay.extend(
                            replay_rows(
                                context,
                                endpoint,
                                mechanism,
                                seed,
                                feature_set,
                                cfg["models"],
                                result,
                                r,
                                oracle_p,
                            )
                        )
                print(f"DONE {endpoint} {feature_set} {mechanism}", flush=True)
    scenarios = pd.DataFrame(rows)
    replay_table = pd.concat(replay, ignore_index=True)
    expected_keys = {
        (endpoint, mechanism, int(seed), feature_set, model)
        for endpoint in cfg["truth_endpoints"]
        for mechanism in cfg["mechanisms"]
        for seed in cfg["mask_seeds"]
        for feature_set in cfg["feature_sets"]
        for model in cfg["models"]
    }
    actual_keys = set(
        scenarios[["endpoint", "mechanism", "seed", "feature_set", "model"]].itertuples(
            index=False, name=None
        )
    )
    require(
        len(scenarios) == 640
        and not scenarios.duplicated(
            ["endpoint", "mechanism", "seed", "feature_set", "model"]
        ).any()
        and actual_keys == expected_keys,
        "scenario grid mismatch",
    )
    shared_contracts = {
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
    require(all(shared_contracts.values()), "shared fold/mask contract mismatch")
    return scenarios, replay_table, fold_hashes, shared_contracts


def f1_reproduction(root, scenarios, cfg):
    pilot = pd.read_csv(root / "artifacts/experiment/pilot_001/full/scenario_results.csv")
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
    source = pd.concat([mcar, mar], ignore_index=True)
    source["seed"] = source["seed"].astype(int)
    source = source.rename(columns={"truth_endpoint": "endpoint"})
    keys = ["endpoint", "diagnostic_mechanism", "seed", "model"]
    expected = {
        (endpoint, mechanism, int(seed), model)
        for endpoint in cfg["truth_endpoints"]
        for mechanism in cfg["mechanisms"]
        for seed in cfg["mask_seeds"]
        for model in cfg["models"]
    }
    source_keys = set(source[keys].itertuples(index=False, name=None))
    diagnostic = scenarios[scenarios["feature_set"] == "full_2267"].copy()
    diagnostic = diagnostic.rename(columns={"mechanism": "diagnostic_mechanism"})
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
        maximum[pilot_column] = float(
            np.max(
                np.abs(
                    merged[f"pilot_{pilot_column}"].to_numpy(float)
                    - merged[f"diagnostic_{diagnostic_column}"].to_numpy(float)
                )
            )
        )
    tolerance = float(cfg["reproduction_tolerance"])
    return {
        "pass": bool(all(value <= tolerance for value in maximum.values())),
        "grid_ok": True,
        "rows": int(len(merged)),
        "tolerance": tolerance,
        "maximum_absolute_difference": maximum,
    }


def f2_replay(replay_path, scenarios, contexts, cfg):
    replay = pd.read_parquet(replay_path)
    contract = cfg["f2_contract"]
    key = contract["unique_key"]
    required_columns = set(key + contract["required_row_columns"])
    cells = ["endpoint", "mechanism", "seed", "feature_set", "model"]
    scenario_numeric = [
        "fitted_ipw",
        "fitted_aipw",
        "fitted_aipw_se",
        "oracle_ipw",
        "oracle_aipw",
        "oracle_aipw_se",
    ]
    scenario_required = set(cells + scenario_numeric + [
        "fold_sha256",
        "mask_sha256",
        "oracle_probability_sha256",
    ])
    if not required_columns.issubset(replay.columns):
        return {"pass": False, "grid_ok": False, "reason": "F2 missing required row column"}
    if not scenario_required.issubset(scenarios.columns):
        return {"pass": False, "grid_ok": False, "reason": "F2 missing required scenario column"}
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
    replay_numeric = [
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
    try:
        replay_values = replay[replay_numeric].to_numpy(float)
        scenario_values = scenarios.loc[
            scenarios["seed"] == int(cfg["row_replay_seed"]), scenario_numeric
        ].to_numpy(float)
    except (TypeError, ValueError):
        return {"pass": False, "grid_ok": True, "reason": "F2 nonnumeric stored value"}
    if not np.isfinite(replay_values).all() or not np.isfinite(scenario_values).all():
        return {"pass": False, "grid_ok": True, "reason": "F2 nonfinite stored value"}
    if not np.isin(replay["R"].to_numpy(float), [0.0, 1.0]).all():
        return {"pass": False, "grid_ok": True, "reason": "F2 R is not binary"}
    for column in ["p_fitted", "p_oracle"]:
        probability = replay[column].to_numpy(float)
        if not ((probability > 0) & (probability <= 1)).all():
            return {"pass": False, "grid_ok": True, "reason": f"F2 invalid {column}"}
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
    vectors_ok = True
    frozen_generator_ok = True
    scenario_hashes_ok = True
    frozen_masks = {}
    for endpoint in cfg["truth_endpoints"]:
        generated = make_masks(contexts[endpoint], cfg)
        for mechanism in cfg["mechanisms"]:
            frozen_masks[(endpoint, mechanism)] = generated[
                (mechanism, int(cfg["row_replay_seed"]))
            ]
    scenario_index = scenarios.set_index(cells)
    for cell, part in replay.groupby(cells, sort=False):
        part = part.sort_values("row_index")
        context = contexts[cell[0]]
        vectors_ok &= bool(
            np.array_equal(part["row_index"].to_numpy(int), context["row_index"])
            and np.array_equal(part["molecule_id"].astype(str).to_numpy(), context["molecule_id"])
            and np.array_equal(part["scaffold_group"].astype(str).to_numpy(), context["groups"].astype(str))
            and np.array_equal(part["scaffold_sha256"].astype(str).to_numpy(), context["scaffold_sha256"])
            and np.array_equal(part["fold_id"].to_numpy(int), context["fold_id"])
        )
        stored = scenario_index.loc[cell]
        r = part["R"].to_numpy(float)
        loss = part["loss"].to_numpy(float)
        mhat = part["mhat"].to_numpy(float)
        groups = part["scaffold_sha256"].astype(str).to_numpy()
        expected_r, expected_oracle = frozen_masks[(cell[0], cell[1])]
        frozen_generator_ok &= bool(
            np.array_equal(r, expected_r.astype(float))
            and np.array_equal(part["p_oracle"].to_numpy(float), expected_oracle)
        )
        scenario_hashes_ok &= bool(
            stored["mask_sha256"] == vector_hash(expected_r, "u1")
            and stored["oracle_probability_sha256"]
            == vector_hash(expected_oracle, "<f8")
            and stored["fold_sha256"] == context["fold_sha256"]
        )
        for suffix, p_column in [("fitted", "p_fitted"), ("oracle", "p_oracle")]:
            p = part[p_column].to_numpy(float)
            augmentation = r * (loss - mhat) / p
            pseudo = mhat + augmentation
            influence = pseudo - pseudo.mean()
            weights = r / p
            ipw = float(np.sum(weights * loss) / weights.sum())
            aipw = float(pseudo.mean())
            se = cluster_se(influence, groups)
            if not np.isfinite(
                np.concatenate(
                    [augmentation, pseudo, influence, np.asarray([ipw, aipw, se])]
                )
            ).all():
                return {
                    "pass": False,
                    "grid_ok": True,
                    "reason": "F2 nonfinite recomputed value",
                }
            maximum[f"augmentation_{suffix}"] = max(
                maximum[f"augmentation_{suffix}"],
                float(np.max(np.abs(augmentation - part[f"augmentation_{suffix}"].to_numpy(float)))),
            )
            maximum[f"pseudo_outcome_{suffix}"] = max(
                maximum[f"pseudo_outcome_{suffix}"],
                float(np.max(np.abs(pseudo - part[f"pseudo_outcome_{suffix}"].to_numpy(float)))),
            )
            maximum[f"influence_{suffix}"] = max(
                maximum[f"influence_{suffix}"],
                float(np.max(np.abs(influence - part[f"influence_{suffix}"].to_numpy(float)))),
            )
            maximum[f"{suffix}_ipw"] = max(
                maximum[f"{suffix}_ipw"], abs(ipw - float(stored[f"{suffix}_ipw"]))
            )
            maximum[f"{suffix}_aipw"] = max(
                maximum[f"{suffix}_aipw"], abs(aipw - float(stored[f"{suffix}_aipw"]))
            )
            maximum[f"{suffix}_aipw_se"] = max(
                maximum[f"{suffix}_aipw_se"], abs(se - float(stored[f"{suffix}_aipw_se"]))
            )
    shared_ok = bool(
        replay.groupby(["endpoint", "mechanism", "row_index"])["R"].nunique().max() == 1
        and replay.groupby(["endpoint", "mechanism", "row_index"])["p_oracle"].nunique().max() == 1
        and replay.groupby(["endpoint", "row_index"])["fold_id"].nunique().max() == 1
        and replay.groupby(["endpoint", "row_index"])["scaffold_sha256"].nunique().max() == 1
    )
    hash_ok = all(
        vector_hash(contexts[endpoint]["fold_id"], "<i8") == contexts[endpoint]["fold_sha256"]
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


def median_bias(frame, column):
    return float(np.median(np.abs(frame[column] - frame["full_truth"])))


def safe_ratio(numerator, denominator):
    return None if denominator == 0 else float(numerator / denominator)


def summarize(scenarios, f1, f2, fold_hashes, shared_contracts, cfg):
    pooled = []
    for (mechanism, feature_set), part in scenarios.groupby(
        ["mechanism", "feature_set"], sort=False
    ):
        pooled.append(
            {
                "mechanism": mechanism,
                "feature_set": feature_set,
                "rows": int(len(part)),
                "median_absolute_bias": {
                    column: median_bias(part, column)
                    for column in [
                        "cc",
                        "or_mean",
                        "fitted_ipw",
                        "fitted_aipw",
                        "oracle_ipw",
                        "oracle_aipw",
                    ]
                },
                "median_fitted_clip_fraction": float(part["fitted_clip_fraction"].median()),
                "median_fitted_ess_fraction": float(part["fitted_ess_fraction"].median()),
                "median_oracle_ess_fraction": float(part["oracle_ess_fraction"].median()),
            }
        )
    lookup = {(x["mechanism"], x["feature_set"]): x for x in pooled}
    thresholds = cfg["signature_thresholds"]
    propensity_details = []
    augmentation_details = []
    dimension_details = []
    for mechanism in cfg["mechanisms"]:
        full = lookup[(mechanism, "full_2267")]
        low = lookup[(mechanism, "lowdim_12")]
        full_bias = full["median_absolute_bias"]
        low_bias = low["median_absolute_bias"]
        propensity_ratio = safe_ratio(
            full_bias["oracle_aipw"], full_bias["fitted_aipw"]
        )
        augmentation_ratio = safe_ratio(
            full_bias["oracle_aipw"], full_bias["oracle_ipw"]
        )
        dimension_ratio = safe_ratio(
            low_bias["fitted_aipw"], full_bias["fitted_aipw"]
        )
        clip_reduction = (
            full["median_fitted_clip_fraction"] - low["median_fitted_clip_fraction"]
        )
        propensity_details.append(
            {
                "mechanism": mechanism,
                "ratio": propensity_ratio,
                "meets": full_bias["oracle_aipw"]
                <= thresholds["oracle_aipw_fraction_of_fitted_aipw"]
                * full_bias["fitted_aipw"],
            }
        )
        augmentation_details.append(
            {
                "mechanism": mechanism,
                "ratio": augmentation_ratio,
                "meets": full_bias["oracle_aipw"]
                > thresholds["oracle_aipw_multiple_of_oracle_ipw"]
                * full_bias["oracle_ipw"],
            }
        )
        dimension_details.append(
            {
                "mechanism": mechanism,
                "bias_ratio": dimension_ratio,
                "clip_fraction_reduction": clip_reduction,
                "meets": low_bias["fitted_aipw"]
                <= thresholds["lowdim_aipw_fraction_of_full"]
                * full_bias["fitted_aipw"]
                and clip_reduction >= thresholds["minimum_clip_fraction_reduction"],
            }
        )
    signatures = {
        "PROPENSITY_INSTABILITY_SUPPORTED": {
            "supported": bool(all(x["meets"] for x in propensity_details)),
            "details": propensity_details,
        },
        "AUGMENTATION_INSTABILITY_SUPPORTED": {
            "supported": bool(all(x["meets"] for x in augmentation_details)),
            "details": augmentation_details,
        },
        "HIGH_DIMENSION_INSTABILITY_SUPPORTED": {
            "supported": bool(all(x["meets"] for x in dimension_details)),
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
        "fold_id_sha256": fold_hashes,
        "fold_id_hash_encoding": "ordered little-endian signed int64 raw bytes",
        "shared_contracts": shared_contracts,
        "F1_reproduction": f1,
        "F2_row_replay": f2,
        "pooled_mechanism_summaries": pooled,
        "signature_thresholds": thresholds,
        "signatures": signatures,
    }


def output_hashes(output):
    return {
        str(path.relative_to(output)).replace("\\", "/"): file_hash(path)
        for path in sorted(output.rglob("*"))
        if path.is_file() and path.name != "run_manifest.json"
    }


def run_full(root, output, cfg, pilot_cfg, input_hashes, started, started_utc):
    _, contexts = load_contexts(root, cfg, pilot_cfg)
    scenarios, replay, fold_hashes, shared_contracts = run_grid(contexts, cfg)
    scenario_path = output / "scenario_decomposition.csv"
    replay_path = output / "row_replay_seed31001.parquet"
    scenarios.to_csv(scenario_path, index=False)
    replay.to_parquet(replay_path, index=False)
    f1 = f1_reproduction(root, scenarios, cfg)
    f2 = f2_replay(replay_path, scenarios, contexts, cfg)
    summary = summarize(scenarios, f1, f2, fold_hashes, shared_contracts, cfg)
    write_json(output / "summary.json", summary)
    verification = (
        "# Forensic diagnostic verification\n\n"
        f"- F1 frozen-pilot reproduction: `{'PASS' if f1['pass'] else 'FAIL'}`.\n"
        f"- F2 reloaded row-level formula replay: `{'PASS' if f2['pass'] else 'FAIL'}`.\n"
        f"- Diagnostic decision: `{summary['decision']}`.\n\n"
        "Scope: F2 reloads the Parquet row table and independently recomputes the frozen "
        "Hájek, AIPW, and scaffold-cluster SE formulas. It is not a second nuisance-model "
        "implementation or an independent scientific validation. The pilot decision remains "
        "`STOP_PILOT_NEGATIVE`.\n"
    )
    (output / "verification.md").write_text(verification, encoding="utf-8")
    manifest = {
        "diagnostic_id": cfg["diagnostic_id"],
        "mode": "full",
        "protocol_amendment": cfg["protocol_amendment"],
        "started_utc": started_utc,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "duration_seconds": float(time.time() - started),
        "command": [sys.executable, *sys.argv],
        "input_hashes": input_hashes,
        "runner_sha256": file_hash(Path(__file__)),
        "fold_id_sha256": fold_hashes,
        "fold_id_hash_encoding": "ordered little-endian signed int64 raw bytes",
        "shared_contracts": shared_contracts,
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
            "rdkit": rdkit.__version__,
            "platform": platform.platform(),
        },
        "feature_sets": {
            endpoint: {
                name: {"shape": list(context["features"][name].shape), "dtype": str(context["features"][name].dtype)}
                for name in cfg["feature_sets"]
            }
            for endpoint, context in contexts.items()
        },
        "result_rows": {
            "scenario_decomposition.csv": int(len(scenarios)),
            "row_replay_seed31001.parquet": int(len(replay)),
        },
        "decision": summary["decision"],
        "output_hashes_before_manifest": output_hashes(output),
    }
    write_json(output / "run_manifest.json", manifest)


def run_smoke(output, cfg, input_hashes):
    rng = np.random.default_rng(20260911)
    X = rng.normal(size=(36, 12))
    groups = np.repeat(np.arange(12), 3).astype(str)
    split, fold_id, fold_sha256 = fold_contract(groups, 3)
    r = np.tile([0, 1, 1], 12)
    losses = np.abs(rng.normal(size=(36, 2)))
    result = crossfit_estimators(
        X,
        groups,
        r,
        losses,
        3,
        cfg["propensity_model"],
        cfg["loss_nuisance"],
        cfg["propensity_clip"],
        prepared=scaled_folds(X, split),
    )
    oracle = estimates_for_probability(
        r, losses[:, 0], result["mhat"][:, 0], np.full(36, 2 / 3), groups
    )
    checks = {
        "status": "PASS",
        "non_adjudicative": True,
        "rows": 36,
        "fold_id_sha256": fold_sha256,
        "folds_complete": bool(np.array_equal(np.sort(np.unique(fold_id)), np.arange(3))),
        "fitted_outputs_finite": bool(
            np.isfinite(result["p_raw"]).all()
            and np.isfinite(result["p"]).all()
            and np.isfinite(result["mhat"]).all()
        ),
        "oracle_outputs_finite": bool(
            all(np.isfinite(oracle[key]) for key in ["ipw", "aipw", "se", "ess_fraction"])
        ),
        "input_hashes": input_hashes,
        "runner_sha256": file_hash(Path(__file__)),
    }
    require(
        checks["folds_complete"]
        and checks["fitted_outputs_finite"]
        and checks["oracle_outputs_finite"],
        "smoke check failed",
    )
    write_json(output / "smoke.json", checks)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full"], required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    cfg_path = root / "artifacts" / "diagnostics" / "forensic_001" / "config.json"
    cfg, pilot_cfg, input_hashes = verify_inputs(root, cfg_path)
    formal_output = (root / "artifacts" / "diagnostics" / "forensic_001" / "full").resolve()
    if args.mode == "smoke":
        require(args.output is not None, "smoke requires an explicit non-formal --output")
        output = args.output.resolve()
        require(output != formal_output, "smoke cannot write the formal result directory")
    else:
        output = args.output.resolve() if args.output else formal_output
    require(not output.exists() or not any(output.iterdir()), f"output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    started = time.time()
    started_utc = datetime.now(timezone.utc).isoformat()
    print(f"START mode={args.mode} protocol={cfg['protocol_sha256']}", flush=True)
    if args.mode == "smoke":
        run_smoke(output, cfg, input_hashes)
    else:
        run_full(root, output, cfg, pilot_cfg, input_hashes, started, started_utc)
    print(f"DONE mode={args.mode} seconds={time.time() - started:.3f}", flush=True)


if __name__ == "__main__":
    main()
