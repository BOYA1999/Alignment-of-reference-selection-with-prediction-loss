import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "artifacts/confirmatory/confirm_001/full_recovery_001"
OUT = ROOT / "artifacts/revision_001/pxr"
SOURCE = ROOT / "data/external/source_freeze_20260911/pxr"
DESCRIPTORS = ["MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors", "NumRotatableBonds", "RingCount", "FractionCSP3", "HeavyAtomCount", "MolMR"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    frozen_before = {p.name: sha(p) for p in FROZEN.iterdir() if p.is_file()}
    assert len(frozen_before) == 19 and not any(p.is_dir() for p in FROZEN.iterdir())
    inputs = [FROZEN / "pxr_linkage_audit.csv", FROZEN / "pxr_overlap_summary.csv", SOURCE / "pxr-challenge_TRAIN.csv", SOURCE / "pxr-challenge_single_concentration_TRAIN.csv", SOURCE / "README.md"]
    linkage, original, dose, screen = [pd.read_csv(p, float_precision="round_trip") for p in inputs[:4]]
    original = original.set_index("specification")
    parent = screen.loc[screen.compound_class.eq("Library") & screen.concentration_M.eq(8.251e-06)]
    assert len(linkage) == len(parent) == 10747 and len(dose) == 4139
    assert not dose.OCNT_ID.duplicated().any()
    assert linkage.OCNT_ID.tolist() == parent.OCNT_ID.tolist()
    observed = linkage.R_followup.to_numpy(bool)
    hit = linkage.formal_hit.to_numpy(bool)
    assert np.array_equal(observed, parent.OCNT_ID.isin(dose.OCNT_ID).to_numpy())
    assert np.array_equal(hit, ((parent.log2_fc_estimate > 1) & (parent.fdr_bh < 0.05)).to_numpy())
    selections = {"all_parent": np.ones(len(linkage), bool), "formal_hit": hit, "non_hit": ~hit, "table_member": observed, "not_table_member": ~observed}
    for h in (0, 1):
        for r in (0, 1):
            selections[f"formal_hit_{h}_table_member_{r}"] = (hit == h) & (observed == r)
    tail_rows, weight_rows, balance_rows = [], [], []
    for spec in ("structure_10", "structure_screen_12"):
        raw = linkage[f"p_raw_{spec}"].to_numpy()
        frozen_p = linkage[f"p_clipped_{spec}"].to_numpy()
        assert np.array_equal(np.clip(raw, 0.05, 0.95), frozen_p)
        assert np.isfinite(raw).all() and ((raw > 0) & (raw <= 1)).all()
        for name, selected in selections.items():
            p = raw[selected]
            low, high = p < 0.05, p > 0.95
            tail_rows.append({"specification": f"pxr_{spec}", "stratum": name, "n": int(selected.sum()), "n_table_member": int(observed[selected].sum()), "low_count": int(low.sum()), "low_fraction": float(low.mean()), "high_count": int(high.sum()), "high_fraction": float(high.mean()), "combined_fraction": float((low | high).mean()), "p_min": float(p.min()), "p_median": float(np.median(p)), "p_max": float(p.max())})
        names = DESCRIPTORS + (["log2_fc_estimate", "fdr_score"] if spec == "structure_screen_12" else [])
        x = linkage[names].to_numpy()
        for rule, p in (("frozen_two_sided_005_095", frozen_p), ("post_result_lower_only_005", np.maximum(raw, 0.05))):
            weights = observed / p
            weighted_mean = weights @ x / weights.sum()
            smds = np.abs(weighted_mean - x.mean(axis=0)) / x.std(axis=0, ddof=0)
            ess = float(weights.sum() ** 2 / (weights @ weights))
            row = {"specification": f"pxr_{spec}", "weight_rule": rule, "n": len(p), "n_table_member": int(observed.sum()), "ess": ess, "ess_fraction": ess / observed.sum(), "max_abs_smd": float(smds.max()), "max_normalized_weight": float(weights.max() / weights.sum()), "low_fraction": float((raw < 0.05).mean()), "high_fraction": float((raw > 0.95).mean())}
            weight_rows.append(row)
            balance_rows.extend({"specification": f"pxr_{spec}", "weight_rule": rule, "feature": name, "abs_smd": float(smd)} for name, smd in zip(names, smds))
            if rule.startswith("frozen"):
                for metric in ("ess", "ess_fraction", "max_abs_smd", "max_normalized_weight"):
                    assert np.isclose(row[metric], original.loc[f"pxr_{spec}", metric], rtol=0, atol=1e-11), metric
                assert np.isclose(row["low_fraction"] + row["high_fraction"], original.loc[f"pxr_{spec}", "clip_fraction"], rtol=0, atol=1e-15)
    joined = linkage[["OCNT_ID", "formal_hit", "R_followup"]].merge(dose, on="OCNT_ID", how="left", validate="one_to_one")
    assert int(joined.pEC50.notna().sum()) <= int(observed.sum())
    endpoints = {"pEC50": ("pEC50", "pEC50_std.error (-log10(molarity))", "pEC50_ci.lower (-log10(molarity))", "pEC50_ci.upper (-log10(molarity))"), "Emax": ("Emax_estimate (log2FC vs. baseline)", "Emax_std.error (log2FC vs. baseline)", "Emax_ci.lower (log2FC vs. baseline)", "Emax_ci.upper (log2FC vs. baseline)"), "relative_Emax": ("Emax.vs.pos.ctrl_estimate (dimensionless)", "Emax.vs.pos.ctrl_std.error (dimensionless)", "Emax.vs.pos.ctrl_ci.lower (dimensionless)", "Emax.vs.pos.ctrl_ci.upper (dimensionless)")}
    subsets = {"all_dose_table": dose, "parent_table_members": joined.loc[observed], "parent_hit_members": joined.loc[observed & hit], "parent_non_hit_members": joined.loc[observed & ~hit]}
    parameter_rows = []
    for subset, data in subsets.items():
        for endpoint, columns in endpoints.items():
            estimate, se, lower, upper = (pd.to_numeric(data[c], errors="coerce").to_numpy() for c in columns)
            finite = np.isfinite(estimate)
            intervals = np.isfinite(lower) & np.isfinite(upper)
            values = estimate[finite]
            parameter_rows.append({"subset": subset, "endpoint": endpoint, "n_rows": len(data), "finite_estimate_count": int(finite.sum()), "nonfinite_estimate_count": int((~finite).sum()), "finite_nonnegative_se_count": int((np.isfinite(se) & (se >= 0)).sum()), "finite_interval_count": int(intervals.sum()), "reversed_interval_count": int((intervals & (lower > upper)).sum()), "estimate_outside_interval_count": int((intervals & finite & ((estimate < lower) | (estimate > upper))).sum()), "negative_estimate_count": int((values < 0).sum()), "zero_estimate_count": int((values == 0).sum()), "minimum": float(values.min()), "q01": float(np.quantile(values, 0.01)), "median": float(np.median(values)), "q99": float(np.quantile(values, 0.99)), "maximum": float(values.max())})
    summary = {
        "analysis_id": "revision_001_pxr",
        "review_items": ["PXR_one_sided_positivity", "PXR_clipping_tails", "PXR_table_membership_vs_endpoint_validity"],
        "design": "Post-result descriptive reaggregation of frozen probabilities and pinned source tables; no nuisance refit or risk estimation.",
        "parent_n": len(parent), "dose_table_n": len(dose), "parent_table_member_n": int(observed.sum()),
        "formal_hit": {"n": int(hit.sum()), "table_member_n": int((hit & observed).sum()), "fraction": float(observed[hit].mean())},
        "non_hit": {"n": int((~hit).sum()), "table_member_n": int((~hit & observed).sum()), "fraction": float(observed[~hit].mean())},
        "stored_upper_extreme_trigger_is_not_missing_outcome_positivity_failure": True,
        "positivity_condition": "p(R=1 | X=x)>0 on the target support; p near 1 alone is not a missing-outcome identification obstacle.",
        "selection_assumption_status": "NOT_ESTABLISHED_BY_PUBLIC_RECORDS",
        "revised_interpretation": "FULL_LIBRARY_RISK_NOT_ESTIMATED_ASSUMPTIONS_AND_LOSS_CONTRACT_UNESTABLISHED",
        "dose_columns": dose.columns.tolist(),
        "curve_class_or_rollover_or_individual_qc_columns_present": [c for c in dose if any(s in c.lower() for s in ("curve", "rollover", "quality", "qc"))],
        "predictor_loss_eligibility_contract": "Not specified for PXR; finite released parameters are not an independently validated loss-eligibility indicator.",
        "input_sha256": {p.relative_to(ROOT).as_posix(): sha(p) for p in inputs},
        "frozen_19_sha256": frozen_before,
        "sources": [{"url": "https://openadmet.ghost.io/announcing-the-next-openadmet-blind-challenge-predicting-pxr-induction/", "retrieved": "2026-09-12", "supports": "additional analog/diversity follow-up, curve classification and constrained Bayesian fitting, post-processing QC and reduced final training size"}, {"url": "https://arxiv.org/html/2606.14506v1#S3", "retrieved": "2026-09-12", "supports": "Assumption 3.3 specifies label-observation positivity as a lower bound; preprint, June 2026"}],
        "script_sha256": sha(Path(__file__)),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    tables = {"propensity_tails.csv": tail_rows, "weighting_sensitivity.csv": weight_rows, "feature_balance.csv": balance_rows, "parameter_completeness.csv": parameter_rows}
    for name, rows in tables.items():
        pd.DataFrame(rows).to_csv(OUT / name, index=False, float_format="%.17g")
    frozen_after = {p.name: sha(p) for p in FROZEN.iterdir() if p.is_file()}
    assert frozen_before == frozen_after
    summary["frozen_19_unchanged"] = True
    summary["output_sha256"] = {name: sha(OUT / name) for name in tables}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "parent_n": len(parent), "frozen_19_unchanged": True, "tables": list(tables)}, indent=2))


if __name__ == "__main__":
    main()
