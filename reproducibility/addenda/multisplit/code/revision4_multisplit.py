import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import time

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import run_confirmatory as core


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "artifacts/revision_004_multisplit"
SCREEN = BASE / "screen"
SMOKE = BASE / "smoke"
FULL = BASE / "full"
SOURCE = ROOT / "evidence/upstream/expansion_ml_comparison/biogen_adme_3521.csv"
CFG = ROOT / "artifacts/confirmatory/confirm_001/config.json"
FROZEN = ROOT / "artifacts/confirmatory/confirm_001/full_recovery_001"
REVISION2 = ROOT / "artifacts/revision_002/original/estimates.csv"
ASSIGNMENTS = SCREEN / "split_assignments.csv"
SEEDS = list(range(61453, 61463))
MASK_SEEDS = list(range(41001, 41101))
CASE = 2945
ENDPOINTS = ["LOG_HLM", "LOG_RLM"]
ANCHORS = ["LOG_SOL", "LOG_MDR1_ER"]
MECHANISMS = ["mcar", "mar"]
MODELS = ["ecfp_rf", "ecfp_desc_xgb"]
EXPECTED_SOURCE_SHA256 = "cd310dd3911a9f9763d0581e6202bc20aaf8a810523d01e9ab5e403681a386b4"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs():
    assert sha(SOURCE) == EXPECTED_SOURCE_SHA256
    split_manifest = json.loads((SCREEN / "split_manifest.json").read_text(encoding="utf-8"))
    assert split_manifest["status"] == "PASS_OUTCOME_BLIND_SPLIT_RECONSTRUCTION"
    assert sha(ASSIGNMENTS) == split_manifest["output_sha256"]["split_assignments.csv"]
    cfg = json.loads(CFG.read_text(encoding="utf-8"))
    environment = core.verify_environment(cfg)
    frame = pd.read_csv(SOURCE)
    assignments = pd.read_csv(ASSIGNMENTS)
    assert len(assignments) == len(frame) * len(SEEDS)
    assert sorted(assignments.split_seed.unique().tolist()) == SEEDS
    original = assignments[assignments.split_seed.eq(61453)].sort_values("source_row_index")
    assert original.ds.tolist() == frame.ds.tolist()
    registry = core.descriptor_registry()
    features = core.build_features(frame.SMILES.tolist(), registry)
    names = [name for name, _ in registry]
    assert names.count("Ipc") == 1 and names.count("AvgIpc") == 1
    return frame, assignments, cfg, features, names, environment


def labels_for(assignments, seed):
    part = assignments[assignments.split_seed.eq(seed)].sort_values("source_row_index")
    assert part.source_row_index.tolist() == list(range(3521))
    return part.ds.to_numpy()


def run_screen():
    path = SCREEN / "descriptor_screen.csv"
    assert not path.exists()
    frame, assignments, _, features, names, environment = load_inputs()
    ipc = features["desc"][:, names.index("Ipc")]
    rows = []
    for seed in SEEDS:
        labels = labels_for(assignments, seed)
        for endpoint in ENDPOINTS:
            values = pd.to_numeric(frame[endpoint], errors="coerce").to_numpy(float)
            train = np.flatnonzero((labels == "train") & np.isfinite(values))
            target = np.flatnonzero((labels == "test") & np.isfinite(values))
            target_max = target[np.argmax(ipc[target])]
            train_max = train[np.argmax(ipc[train])]
            rows.append({
                "split_seed": seed,
                "endpoint": endpoint,
                "train_finite": len(train),
                "target_finite": len(target),
                "case_2945_in_target": bool(CASE in target),
                "train_ipc_max": ipc[train_max],
                "train_ipc_max_source_row": int(train_max),
                "target_ipc_max": ipc[target_max],
                "target_ipc_max_source_row": int(target_max),
                "target_to_train_ipc_max_ratio": ipc[target_max] / ipc[train_max],
            })
    result = pd.DataFrame(rows)
    result.to_csv(path, index=False, float_format="%.17g")
    manifest = {
        "status": "PASS_SCREEN_CLEAR_TO_FULL",
        "screen_role": "outcome_blind_partition_and_descriptor_support_only",
        "seeds": SEEDS,
        "case_in_target_splits": sorted(result.loc[result.case_2945_in_target, "split_seed"].unique().tolist()),
        "minimum_target_endpoint_rows": int(result.target_finite.min()),
        "all_target_endpoint_rows_at_least_900": bool((result.target_finite >= 900).all()),
        "environment": environment,
        "source_hashes": {
            "protocol.md": sha(BASE / "protocol.md"),
            "split_assignments.csv": sha(ASSIGNMENTS),
            "source_derivative": sha(SOURCE),
            "script": sha(Path(__file__)),
        },
        "output_sha256": sha(path),
    }
    (SCREEN / "screen_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: manifest[key] for key in ["status", "case_in_target_splits", "minimum_target_endpoint_rows"]}, indent=2))


def build_context(frame, labels, features, names, cfg, endpoint):
    values = pd.to_numeric(frame[endpoint], errors="coerce").to_numpy(float)
    train = np.flatnonzero((labels == "train") & np.isfinite(values))
    target = np.flatnonzero((labels == "test") & np.isfinite(values))
    rf_cfg = cfg["prediction_models"]["ecfp_rf"]
    rf = core.RandomForestRegressor(
        n_estimators=rf_cfg["n_estimators"], criterion=rf_cfg["criterion"],
        min_samples_leaf=rf_cfg["min_samples_leaf"], max_features=rf_cfg["max_features"],
        n_jobs=rf_cfg["n_jobs"], random_state=rf_cfg["random_state"],
    ).fit(features["fps"][train], values[train])
    full_train = np.hstack([features["fps"][train], features["desc"][train]]).astype(np.float32)
    full_target = np.hstack([features["fps"][target], features["desc"][target]]).astype(np.float32)
    xgb_cfg = cfg["prediction_models"]["ecfp_desc_xgb"]
    xgb = core.XGBRegressor(**{key: value for key, value in xgb_cfg.items() if key != "features"})
    xgb.fit(full_train, values[train])
    predictions = np.column_stack([
        rf.predict(features["fps"][target]),
        xgb.predict(full_target),
    ]).astype(float)
    similarity = core.max_similarity(
        [features["bitvectors"][i] for i in target],
        [features["bitvectors"][i] for i in train],
    )
    descriptor_index = {name: index for index, name in enumerate(names)}
    low = np.column_stack([
        features["desc"][target][:, [descriptor_index[name] for name in core.LOW_DESCRIPTOR_NAMES]],
        similarity,
    ]).astype(float)
    raw = np.hstack([features["fps"][target], features["desc"][target], similarity[:, None]]).astype(float)
    drop = np.delete(raw, 2048 + descriptor_index["Ipc"], axis=1)
    groups = core.scaffold_groups([features["molecules"][i] for i in target], target)
    split, fold = core.fold_contract(groups)
    losses = np.abs(values[target, None] - predictions)
    return {
        "endpoint": endpoint,
        "source_rows": target,
        "groups": groups,
        "fold": fold,
        "split": split,
        "low": low,
        "prepared": {
            "low": core.scaled_folds(low, split),
            "raw": core.scaled_folds(raw, split),
            "drop": core.scaled_folds(drop, split),
        },
        "losses": losses,
        "availability": {
            anchor: np.isfinite(pd.to_numeric(frame.iloc[target][anchor], errors="coerce").to_numpy(float)).astype(int)
            for anchor in ANCHORS
        },
        "representations": {"lowdim_11": low},
        "audit_features": low,
        "audit_names": [*core.LOW_DESCRIPTOR_NAMES, "max_train_similarity"],
        "predictions": predictions,
        "reference": values[target],
        "train_rows": train,
    }


def one_mask(context, split_seed, case_in_target, anchor, mechanism, seed, probability):
    r = (np.random.default_rng(seed).random(len(probability)) < probability).astype(int)
    support = core.fold_support(r, context["split"])["fold_class_support"]
    prop = core.fit_propensity(context["prepared"]["low"], r)
    fitted = {
        variant: core.fit_loss(context["prepared"][variant], r, context["losses"])
        for variant in ["low", "raw", "drop"]
    }
    case_position = np.flatnonzero(context["source_rows"] == CASE)
    case_position = int(case_position[0]) if len(case_position) else None
    rows = []
    for model_index, model in enumerate(MODELS):
        loss = context["losses"][:, model_index]
        truth = float(loss.mean())
        observed = r.astype(bool)
        complete_case = float(loss[observed].mean())
        weights = r / prop["p"]
        hajek = float(np.sum(weights * loss) / np.sum(weights))
        base = {
            "split_seed": split_seed, "endpoint": context["endpoint"], "anchor": anchor,
            "mechanism": mechanism, "mask_seed": seed, "model": model, "truth": truth,
            "case_2945_in_target": case_in_target, "case_R": int(r[case_position]) if case_position is not None else np.nan,
            "fold_class_support": support, "propensity_clip_fraction": prop["clip_fraction"],
        }
        for method, estimate in [("complete_case", complete_case), ("hajek", hajek)]:
            rows.append({**base, "method": method, "loss_variant": "none", "estimate": estimate,
                         "error": estimate - truth, "abs_error": abs(estimate - truth),
                         "largest_source_row": np.nan, "top1_abs_share": np.nan, "case_abs_rank": np.nan})
        for variant, mhat in fitted.items():
            m = mhat[:, model_index]
            pseudo = m + weights * (loss - m)
            estimate = float(pseudo.mean())
            shares = np.abs(pseudo) / np.abs(pseudo).sum()
            largest = int(np.argmax(shares))
            case_rank = np.nan if case_position is None else int(1 + np.sum(shares > shares[case_position]))
            rows.append({**base, "method": "aipw", "loss_variant": variant, "estimate": estimate,
                         "error": estimate - truth, "abs_error": abs(estimate - truth),
                         "largest_source_row": int(context["source_rows"][largest]),
                         "top1_abs_share": float(shares[largest]), "case_abs_rank": case_rank})
    return rows


def verify_original(predictions, estimates):
    frozen_predictions = pd.read_parquet(FROZEN / "predictions.parquet")
    expected = frozen_predictions[["endpoint", "model", "source_row_index", "prediction"]]
    observed = predictions[predictions.split_seed.eq(61453)][["endpoint", "model", "source_row_index", "prediction"]]
    joined = observed.merge(expected, on=["endpoint", "model", "source_row_index"], suffixes=("_new", "_old"), validate="one_to_one")
    prediction_difference = float(np.max(np.abs(joined.prediction_new - joined.prediction_old)))

    old = pd.read_csv(REVISION2, float_precision="round_trip")
    old = old[
        ((old.method == "aipw") & (old.propensity == "low") & old.loss_variant.isin(["low", "raw", "drop"])) |
        ((old.method == "hajek") & (old.propensity == "low") & (old.loss_variant == "none")) |
        ((old.method == "complete_case") & (old.propensity == "none") & (old.loss_variant == "none"))
    ].copy()
    old = old.rename(columns={"seed": "mask_seed"})
    new = estimates[estimates.split_seed.eq(61453)].copy()
    keys = ["endpoint", "anchor", "mechanism", "mask_seed", "model", "method", "loss_variant"]
    compared = new.merge(old[keys + ["estimate"]], on=keys, suffixes=("_new", "_old"), validate="one_to_one")
    risk_difference = float(np.max(np.abs(compared.estimate_new - compared.estimate_old)))
    return {
        "status": "PASS_EXACT_OR_TOLERANCE" if prediction_difference <= 1e-12 and risk_difference <= 1e-10 else "FAIL",
        "prediction_rows": len(joined),
        "risk_rows": len(compared),
        "maximum_prediction_difference": prediction_difference,
        "maximum_risk_estimate_difference": risk_difference,
        "prediction_tolerance": 1e-12,
        "risk_tolerance": 1e-10,
    }


def run_smoke():
    assert not SMOKE.exists()
    screen_manifest = json.loads((SCREEN / "screen_manifest.json").read_text(encoding="utf-8"))
    assert screen_manifest["status"] == "PASS_SCREEN_CLEAR_TO_FULL"
    frame, assignments, cfg, features, names, environment = load_inputs()
    labels = labels_for(assignments, 61453)
    estimates, predictions = [], []
    with threadpool_limits(limits=1):
        for endpoint in ENDPOINTS:
            context = build_context(frame, labels, features, names, cfg, endpoint)
            for model_index, model in enumerate(MODELS):
                for cohort_index, source_row in enumerate(context["source_rows"]):
                    predictions.append({
                        "split_seed": 61453, "endpoint": endpoint, "model": model,
                        "source_row_index": int(source_row), "reference": context["reference"][cohort_index],
                        "prediction": context["predictions"][cohort_index, model_index],
                        "loss": context["losses"][cohort_index, model_index],
                    })
            for anchor in ANCHORS:
                generators = core.generator_contract(context, anchor)
                for mechanism in MECHANISMS:
                    estimates.extend(one_mask(
                        context, 61453, True, anchor, mechanism, 41001, generators[mechanism]["p"]
                    ))
    estimates = pd.DataFrame(estimates)
    predictions = pd.DataFrame(predictions)
    verification = verify_original(predictions, estimates)
    SMOKE.mkdir()
    estimates.to_csv(SMOKE / "estimates.csv", index=False, float_format="%.17g")
    predictions.to_csv(SMOKE / "predictions.csv", index=False, float_format="%.17g")
    record = {
        "status": "PASS_SMOKE_CLEAR_TO_FULL" if verification["status"].startswith("PASS") else "FAIL",
        "verification": verification,
        "environment": environment,
        "source_hashes": {"protocol.md": sha(BASE / "protocol.md"), "script": sha(Path(__file__))},
        "output_sha256": {"estimates.csv": sha(SMOKE / "estimates.csv"), "predictions.csv": sha(SMOKE / "predictions.csv")},
    }
    (SMOKE / "manifest.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))
    assert record["status"] == "PASS_SMOKE_CLEAR_TO_FULL"


def summarize(estimates, predictions):
    group = ["split_seed", "case_2945_in_target", "endpoint", "anchor", "mechanism", "model", "method", "loss_variant"]
    summaries = estimates.groupby(group, dropna=False).agg(
        masks=("mask_seed", "size"), median_abs_error=("abs_error", "median"),
        mean_abs_error=("abs_error", "mean"), p10_abs_error=("abs_error", lambda x: np.quantile(x, .1)),
        p90_abs_error=("abs_error", lambda x: np.quantile(x, .9)),
        median_signed_error=("error", "median"), all_fold_support=("fold_class_support", "all"),
    ).reset_index()
    aipw = summaries[(summaries.method == "aipw")].copy()
    pivot = aipw.pivot(index=group[:-2], columns="loss_variant", values="median_abs_error").reset_index()
    pivot.columns.name = None
    pivot["raw_to_drop_error_ratio"] = pivot["raw"] / pivot["drop"]
    pivot["raw_to_low_error_ratio"] = pivot["raw"] / pivot["low"]

    mar = pivot[pivot.mechanism.eq("mar")]
    present = mar[mar.case_2945_in_target]
    absent = mar[~mar.case_2945_in_target]
    conditional = (present.raw_to_drop_error_ratio.ge(10) & present.raw_to_low_error_ratio.ge(10))
    raw_rows = estimates[
        estimates.method.eq("aipw") & estimates.loss_variant.eq("raw") & estimates.mechanism.eq("mar") &
        estimates.case_2945_in_target
    ]
    rank_one = raw_rows.case_abs_rank.eq(1)
    split_medians = estimates[
        estimates.method.eq("aipw") & estimates.loss_variant.eq("raw") & estimates.mechanism.eq("mar")
    ].groupby(["split_seed", "case_2945_in_target"]).abs_error.median().reset_index(name="median_raw_abs_error")
    state_medians = split_medians.groupby("case_2945_in_target").median_raw_abs_error.median()
    partition_ratio = float(state_medians.get(True, np.nan) / state_medians.get(False, np.nan))
    broader = absent.raw_to_drop_error_ratio.ge(10)

    risks = predictions.groupby(["split_seed", "endpoint", "model"], as_index=False).agg(
        n=("source_row_index", "size"), finite_cohort_mae=("loss", "mean")
    )
    risk_pivot = risks.pivot(index=["split_seed", "endpoint"], columns="model", values="finite_cohort_mae").reset_index()
    risk_pivot.columns.name = None
    risk_pivot["lower_mae_model"] = np.where(risk_pivot.ecfp_rf < risk_pivot.ecfp_desc_xgb, "ecfp_rf", "ecfp_desc_xgb")
    decision = {
        "conditional_mechanism_recurrence": {
            "supported": bool(conditional.mean() >= .80 and rank_one.mean() >= .95),
            "eligible_mar_cells": int(len(present)),
            "cells_meeting_both_10x_rules": int(conditional.sum()),
            "cell_fraction": float(conditional.mean()),
            "raw_estimates_with_case": int(len(raw_rows)),
            "case_rank_one_count": int(rank_one.sum()),
            "case_rank_one_fraction": float(rank_one.mean()),
        },
        "partition_sensitivity": {
            "supported": bool(partition_ratio >= 10),
            "median_raw_error_case_in_target": float(state_medians.get(True, np.nan)),
            "median_raw_error_case_not_in_target": float(state_medians.get(False, np.nan)),
            "ratio": partition_ratio,
        },
        "broader_raw_ipc_instability_without_case": {
            "supported": bool(broader.mean() >= .50),
            "eligible_mar_cells": int(len(absent)),
            "cells_raw_to_drop_at_least_10x": int(broader.sum()),
            "cell_fraction": float(broader.mean()),
        },
        "predictor_ordering": {
            "split_endpoint_comparisons": int(len(risk_pivot)),
            "ecfp_rf_lower_mae": int(risk_pivot.lower_mae_model.eq("ecfp_rf").sum()),
            "ecfp_desc_xgb_lower_mae": int(risk_pivot.lower_mae_model.eq("ecfp_desc_xgb").sum()),
        },
    }
    return summaries, pivot, split_medians, risks, risk_pivot, decision


def run_full():
    assert not FULL.exists()
    screen_manifest = json.loads((SCREEN / "screen_manifest.json").read_text(encoding="utf-8"))
    assert screen_manifest["status"] == "PASS_SCREEN_CLEAR_TO_FULL"
    smoke_manifest = json.loads((SMOKE / "manifest.json").read_text(encoding="utf-8"))
    assert smoke_manifest["status"] == "PASS_SMOKE_CLEAR_TO_FULL"
    FULL.mkdir(parents=True)
    (FULL / "by_split").mkdir()
    started = time.time()
    frame, assignments, cfg, features, names, environment = load_inputs()
    initial_hashes = {path.name: sha(path) for path in FROZEN.iterdir() if path.is_file()}
    all_estimates, all_predictions = [], []
    with threadpool_limits(limits=1):
        for split_index, split_seed in enumerate(SEEDS, 1):
            labels = labels_for(assignments, split_seed)
            split_estimates, split_predictions = [], []
            for endpoint in ENDPOINTS:
                context = build_context(frame, labels, features, names, cfg, endpoint)
                for model_index, model in enumerate(MODELS):
                    for cohort_index, source_row in enumerate(context["source_rows"]):
                        split_predictions.append({
                            "split_seed": split_seed, "endpoint": endpoint, "model": model,
                            "source_row_index": int(source_row), "reference": context["reference"][cohort_index],
                            "prediction": context["predictions"][cohort_index, model_index],
                            "loss": context["losses"][cohort_index, model_index],
                        })
                generators = {
                    anchor: core.generator_contract(context, anchor)
                    for anchor in ANCHORS
                }
                tasks = [
                    (context, split_seed, bool(labels[CASE] == "test"), anchor, mechanism, mask_seed,
                     generators[anchor][mechanism]["p"])
                    for anchor in ANCHORS for mechanism in MECHANISMS for mask_seed in MASK_SEEDS
                ]
                with ThreadPoolExecutor(max_workers=4) as pool:
                    futures = [pool.submit(one_mask, *task) for task in tasks]
                    for done, future in enumerate(as_completed(futures), 1):
                        split_estimates.extend(future.result())
                        if done % 100 == 0:
                            message = f"split {split_seed} {endpoint}: {done}/{len(tasks)} masks; {time.time()-started:.1f}s"
                            print(message, flush=True)
                            with (FULL / "progress.log").open("a", encoding="utf-8") as handle:
                                handle.write(message + "\n")
            split_estimates = pd.DataFrame(split_estimates)
            split_predictions = pd.DataFrame(split_predictions)
            split_estimates.to_csv(FULL / "by_split" / f"estimates_{split_seed}.csv", index=False, float_format="%.17g")
            split_predictions.to_csv(FULL / "by_split" / f"predictions_{split_seed}.csv", index=False, float_format="%.17g")
            all_estimates.append(split_estimates)
            all_predictions.append(split_predictions)
            print(f"completed split {split_index}/{len(SEEDS)}: {split_seed}", flush=True)
    estimates = pd.concat(all_estimates, ignore_index=True)
    predictions = pd.concat(all_predictions, ignore_index=True)
    estimates.to_csv(FULL / "estimates.csv", index=False, float_format="%.17g")
    predictions.to_csv(FULL / "predictions.csv", index=False, float_format="%.17g")
    verification = verify_original(predictions, estimates)
    summaries, comparisons, split_medians, risks, risk_pivot, decision = summarize(estimates, predictions)
    summaries.to_csv(FULL / "summary.csv", index=False, float_format="%.17g")
    comparisons.to_csv(FULL / "aipw_loss_comparisons.csv", index=False, float_format="%.17g")
    split_medians.to_csv(FULL / "split_raw_error_medians.csv", index=False, float_format="%.17g")
    risks.to_csv(FULL / "predictor_risk.csv", index=False, float_format="%.17g")
    risk_pivot.to_csv(FULL / "predictor_ordering.csv", index=False, float_format="%.17g")
    (FULL / "verification.json").write_text(json.dumps(verification, indent=2) + "\n", encoding="utf-8")
    (FULL / "decision.json").write_text(json.dumps(decision, indent=2) + "\n", encoding="utf-8")
    current_hashes = {path.name: sha(path) for path in FROZEN.iterdir() if path.is_file()}
    manifest = {
        "status": "PASS_FULL_VERIFIED" if verification["status"].startswith("PASS") and initial_hashes == current_hashes else "FAIL",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "duration_seconds": time.time() - started,
        "split_seeds": SEEDS,
        "mask_seeds": [MASK_SEEDS[0], MASK_SEEDS[-1]],
        "estimate_rows": len(estimates),
        "prediction_rows": len(predictions),
        "all_fold_support": bool(estimates.fold_class_support.all()),
        "environment": {**environment, "platform": platform.platform(), "xgboost": importlib.metadata.version("xgboost")},
        "frozen_parent_unchanged": initial_hashes == current_hashes,
        "source_hashes": {
            "protocol.md": sha(BASE / "protocol.md"), "script": sha(Path(__file__)),
            "source_derivative": sha(SOURCE), "split_assignments.csv": sha(ASSIGNMENTS),
            "confirmatory_config": sha(CFG), "frozen_predictions": sha(FROZEN / "predictions.parquet"),
            "revision2_estimates": sha(REVISION2),
        },
        "output_sha256": {
            path.name: sha(path) for path in FULL.iterdir() if path.is_file() and path.name != "manifest.json"
        },
    }
    (FULL / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"manifest": manifest["status"], "verification": verification, "decision": decision}, indent=2), flush=True)
    assert manifest["status"] == "PASS_FULL_VERIFIED"


def finalize_recovered_full():
    assert FULL.exists() and not (FULL / "manifest.json").exists()
    assert (FULL / "estimates.csv").exists() and (FULL / "predictions.csv").exists()
    estimates = pd.read_csv(FULL / "estimates.csv", float_precision="round_trip")
    predictions = pd.read_csv(FULL / "predictions.csv", float_precision="round_trip")
    assert len(estimates) == 80000 and estimates.split_seed.nunique() == 10
    verification = verify_original(predictions, estimates)
    summaries, comparisons, split_medians, risks, risk_pivot, decision = summarize(estimates, predictions)
    summaries.to_csv(FULL / "summary.csv", index=False, float_format="%.17g")
    comparisons.to_csv(FULL / "aipw_loss_comparisons.csv", index=False, float_format="%.17g")
    split_medians.to_csv(FULL / "split_raw_error_medians.csv", index=False, float_format="%.17g")
    risks.to_csv(FULL / "predictor_risk.csv", index=False, float_format="%.17g")
    risk_pivot.to_csv(FULL / "predictor_ordering.csv", index=False, float_format="%.17g")
    (FULL / "verification.json").write_text(json.dumps(verification, indent=2) + "\n", encoding="utf-8")
    (FULL / "decision.json").write_text(json.dumps(decision, indent=2) + "\n", encoding="utf-8")
    parent_hashes = {path.name: sha(path) for path in FROZEN.iterdir() if path.is_file()}
    manifest = {
        "status": "PASS_FULL_VERIFIED_AFTER_SUMMARY_RECOVERY" if verification["status"].startswith("PASS") and estimates.fold_class_support.all() else "FAIL",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "recovery_note": "All ten split computations and primary CSV writes completed before a postprocessing column-name collision; summary-only recovery reused those preserved CSVs without refitting.",
        "split_seeds": SEEDS,
        "mask_seeds": [MASK_SEEDS[0], MASK_SEEDS[-1]],
        "estimate_rows": len(estimates),
        "prediction_rows": len(predictions),
        "all_fold_support": bool(estimates.fold_class_support.all()),
        "parent_frozen_file_count": len(parent_hashes),
        "parent_frozen_sha256": parent_hashes,
        "source_hashes": {
            "protocol.md": sha(BASE / "protocol.md"), "script": sha(Path(__file__)),
            "source_derivative": sha(SOURCE), "split_assignments.csv": sha(ASSIGNMENTS),
            "confirmatory_config": sha(CFG), "frozen_predictions": sha(FROZEN / "predictions.parquet"),
            "revision2_estimates": sha(REVISION2),
        },
        "output_sha256": {
            path.name: sha(path) for path in FULL.iterdir() if path.is_file() and path.name != "manifest.json"
        },
    }
    (FULL / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"manifest": manifest["status"], "verification": verification, "decision": decision}, indent=2))
    assert manifest["status"] == "PASS_FULL_VERIFIED_AFTER_SUMMARY_RECOVERY"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["screen", "smoke", "full", "finalize"], required=True)
    args = parser.parse_args()
    {"screen": run_screen, "smoke": run_smoke, "full": run_full, "finalize": finalize_recovered_full}[args.phase]()
