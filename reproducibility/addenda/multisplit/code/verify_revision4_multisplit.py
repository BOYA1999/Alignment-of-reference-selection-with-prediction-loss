import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "artifacts/revision_004_multisplit"
FULL = BASE / "full"
FROZEN = ROOT / "artifacts/confirmatory/confirm_001/full_recovery_001"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    manifest = json.loads((FULL / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "PASS_FULL_VERIFIED_AFTER_SUMMARY_RECOVERY"
    for name, digest in manifest["output_sha256"].items():
        assert sha(FULL / name) == digest
    for name, digest in manifest["parent_frozen_sha256"].items():
        assert sha(FROZEN / name) == digest

    estimates = pd.read_csv(FULL / "estimates.csv", float_precision="round_trip")
    predictions = pd.read_csv(FULL / "predictions.csv", float_precision="round_trip")
    assert len(estimates) == 80000 and len(predictions) == manifest["prediction_rows"] == 37170
    assert sorted(estimates.split_seed.unique()) == list(range(61453, 61463))
    assert sorted(estimates.mask_seed.unique()) == list(range(41001, 41101))
    assert estimates.fold_class_support.all()
    assert np.isfinite(estimates[["truth", "estimate", "error", "abs_error"]].to_numpy()).all()

    keys = ["split_seed", "endpoint", "anchor", "mechanism", "mask_seed", "model", "method", "loss_variant"]
    assert not estimates.duplicated(keys).any()
    assert estimates.groupby(["split_seed", "endpoint", "anchor", "mechanism", "mask_seed"]).size().eq(10).all()
    assert predictions.groupby(["split_seed", "endpoint", "model"]).size().between(908, 953).all()

    parts = pd.concat(
        [pd.read_csv(FULL / "by_split" / f"estimates_{seed}.csv", float_precision="round_trip") for seed in range(61453, 61463)],
        ignore_index=True,
    )
    left = estimates.sort_values(keys).reset_index(drop=True)
    right = parts.sort_values(keys).reset_index(drop=True)
    assert left[keys].equals(right[keys])
    numeric = ["truth", "estimate", "error", "abs_error", "top1_abs_share"]
    assert np.allclose(left[numeric], right[numeric], rtol=0, atol=0, equal_nan=True)

    aipw = estimates[(estimates.method == "aipw") & (estimates.mechanism == "mar")]
    cell = aipw.groupby([
        "split_seed", "case_2945_in_target", "endpoint", "anchor", "model", "loss_variant"
    ]).abs_error.median().unstack("loss_variant")
    cell["raw_drop"] = cell["raw"] / cell["drop"]
    cell["raw_low"] = cell["raw"] / cell["low"]
    present = cell[cell.index.get_level_values("case_2945_in_target")]
    absent = cell[~cell.index.get_level_values("case_2945_in_target")]
    conditional = (present.raw_drop >= 10) & (present.raw_low >= 10)
    broader = absent.raw_drop >= 10
    raw_present = aipw[(aipw.loss_variant == "raw") & aipw.case_2945_in_target]
    rank_one = raw_present.case_abs_rank.eq(1)
    split_medians = aipw[aipw.loss_variant.eq("raw")].groupby(
        ["split_seed", "case_2945_in_target"]
    ).abs_error.median()
    state = split_medians.groupby(level="case_2945_in_target").median()
    partition_ratio = float(state[True] / state[False])

    ordering = pd.read_csv(FULL / "predictor_ordering.csv")
    decision = json.loads((FULL / "decision.json").read_text(encoding="utf-8"))
    recomputed = {
        "conditional_cells": int(len(present)),
        "conditional_pass": int(conditional.sum()),
        "conditional_fraction": float(conditional.mean()),
        "raw_case_estimates": int(len(raw_present)),
        "rank_one": int(rank_one.sum()),
        "rank_one_fraction": float(rank_one.mean()),
        "absent_cells": int(len(absent)),
        "absent_10x": int(broader.sum()),
        "absent_fraction": float(broader.mean()),
        "partition_ratio": partition_ratio,
        "xgb_lower_mae": int(ordering.lower_mae_model.eq("ecfp_desc_xgb").sum()),
    }
    assert recomputed["conditional_cells"] == decision["conditional_mechanism_recurrence"]["eligible_mar_cells"] == 32
    assert recomputed["conditional_pass"] == decision["conditional_mechanism_recurrence"]["cells_meeting_both_10x_rules"] == 32
    assert recomputed["rank_one"] == decision["conditional_mechanism_recurrence"]["case_rank_one_count"] == 3200
    assert recomputed["absent_cells"] == decision["broader_raw_ipc_instability_without_case"]["eligible_mar_cells"] == 48
    assert recomputed["absent_10x"] == decision["broader_raw_ipc_instability_without_case"]["cells_raw_to_drop_at_least_10x"] == 0
    assert abs(recomputed["partition_ratio"] - decision["partition_sensitivity"]["ratio"]) < 1e-12
    assert recomputed["xgb_lower_mae"] == decision["predictor_ordering"]["ecfp_desc_xgb_lower_mae"] == 20

    report = {
        "status": "PASS_INDEPENDENT_AGGREGATE_VERIFICATION",
        "manifest_and_parent_hashes_verified": True,
        "aggregate_matches_ten_split_files_exactly": True,
        "grid_complete_and_unique": True,
        "all_fold_support": True,
        "all_required_estimates_finite": True,
        "recomputed": recomputed,
    }
    path = FULL / "independent_verification.json"
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
