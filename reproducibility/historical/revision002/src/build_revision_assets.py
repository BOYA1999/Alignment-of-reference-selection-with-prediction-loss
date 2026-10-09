import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import build_paper_assets as base


ROOT = Path(__file__).resolve().parents[1]
REVISION = ROOT / "artifacts/revision_001"
TABLES = ROOT / "paper/tables"
SOURCE = ROOT / "paper/source_data"
QA = ROOT / "paper/qa"
KEYS = ["endpoint", "anchor", "model", "mechanism"]
LABELS = {"LOG_HLM": "HLM", "LOG_RLM": "RLM", "LOG_SOL": "SOL", "LOG_MDR1_ER": "MDR1", "ecfp_rf": "RF400", "ecfp_desc_xgb": "XGB800", "mar": "MAR", "mcar": "MCAR"}
MECHANISMS = {"original_mcar": "Original MCAR", "original_mar": "Original MAR", "loss_rank_beta_0.5": "Loss rank, β=0.5", "loss_rank_beta_2.0": "Loss rank, β=2.0"}


def readable(frame, columns, digits=4):
    result = frame[list(columns)].rename(columns=columns).copy()
    for name in result:
        if pd.api.types.is_numeric_dtype(result[name]) and not pd.api.types.is_integer_dtype(result[name]):
            result[name] = result[name].map(lambda value: f"{value:.{digits}f}" if pd.notna(value) else "NA")
        else:
            result[name] = result[name].replace(LABELS)
    return result


def value_range(series, digits=3):
    return f"[{series.min():.{digits}f}, {series.max():.{digits}f}]"


def figure4(calibration, forensic_summary):
    base.style()
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.1))
    fig.subplots_adjust(left=0.07, right=0.985, top=0.85, bottom=0.235, wspace=0.58)
    colors = base.COLORS
    offsets = np.linspace(-0.15, 0.15, 8)
    mechanisms = ["original_mar", "loss_rank_beta_0.5", "loss_rank_beta_2.0"]
    calibration = calibration.sort_values(KEYS)
    ax = axes[0]
    for index, mechanism in enumerate(mechanisms):
        values = calibration.loc[calibration.mechanism.eq(mechanism), "median_abs_error_ratio_vs_cc"].to_numpy()
        assert len(values) == 8
        ax.scatter(index + offsets, values, s=17, color=colors["teal"], edgecolor="white", linewidth=0.35, zorder=3)
    ax.axhline(0.75, color=colors["red"], linestyle="--", linewidth=0.9)
    ax.text(2.39, 0.77, "0.75", color=colors["red"], ha="right", va="bottom", fontsize=6.8)
    ax.set(xlim=(-0.45, 2.45), ylim=(0, 1.2), xticks=range(3), xticklabels=["Original\nMAR", "Loss rank\nβ = 0.5", "Loss rank\nβ = 2.0"], ylabel="Median absolute-error ratio vs CC")
    ax.set_title("a  Selection–loss alignment", loc="left", fontweight="bold", fontsize=8.2, pad=12)
    ax.text(0.5, -0.31, "True p; fixed replay loss model", transform=ax.transAxes, ha="center", fontsize=7, color=colors["grey"])
    ax = axes[1]
    for index, row in forensic_summary.iterrows():
        ax.plot([index - 0.10, index + 0.10], [row.stored_median, row.clipped_median], color="#BEC3C8", lw=1.1, zorder=1)
        ax.scatter(index - 0.10, row.stored_median, s=29, marker="o", color=colors["red"], zorder=3)
        ax.scatter(index + 0.10, row.clipped_median, s=29, marker="s", color=colors["teal"], zorder=3)
    ax.set_yscale("log")
    ax.set(xlim=(-0.45, 2.45), ylim=(0.002, 2500), xticks=range(3), xticklabels=["Low P", "Full P", "Oracle p"], ylabel="Median absolute risk error")
    ax.set_title("b  Descriptor extrapolation", loc="left", fontweight="bold", fontsize=8.2, pad=12)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=colors["red"], label="Stored", markersize=4), Line2D([], [], marker="s", ls="", color=colors["teal"], label="IPC clipped", markersize=4)], loc="center", bbox_to_anchor=(0.5, 0.52), frameon=False, fontsize=6.8, handletextpad=0.3, borderpad=0.1)
    ax.text(0.5, -0.31, "Seed 41,001; 16 estimates per point", transform=ax.transAxes, ha="center", fontsize=7, color=colors["grey"])
    ax = axes[2]
    original_mar = calibration.loc[calibration.mechanism.eq("original_mar")]
    for offset, (_, row) in zip(offsets, original_mar.iterrows()):
        values = [row.mean_se_over_mc_sd, row.legacy_cluster_mean_se_over_mc_sd]
        ax.plot([offset, 1 + offset], values, color="#D2D5D8", lw=0.6, zorder=1)
        ax.scatter(offset, values[0], color=colors["teal"], marker="s", s=17, edgecolor="white", linewidth=0.35, zorder=3)
        ax.scatter(1 + offset, values[1], color=colors["grey"], marker="o", s=17, edgecolor="white", linewidth=0.35, zorder=3)
    ax.axhline(1, color=colors["red"], linestyle="--", linewidth=0.9)
    ax.set(xlim=(-0.4, 1.4), ylim=(0.8, 2.65), xticks=[0, 1], xticklabels=["Mask-design\nvariance", "Legacy\ncluster SE"], ylabel="Mean estimated SE / Monte Carlo SD")
    ax.set_title("c  Uncertainty target", loc="left", fontweight="bold", fontsize=8.2, pad=12)
    ax.text(0.5, -0.31, "Original MAR; fixed loss model", transform=ax.transAxes, ha="center", fontsize=7, color=colors["grey"])
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#E6E8EA", linewidth=0.45, zorder=0)
        ax.tick_params(axis="x", length=0, pad=6)
        ax.tick_params(axis="y", width=0.6)
    return fig


def main():
    input_paths = {
        "estimators": REVISION / "risk/original_estimator_summary.csv",
        "candidate": REVISION / "risk/original_candidate_16_cells.csv",
        "calibration": REVISION / "risk/mechanism_calibration_summary.csv",
        "nonminimum": REVISION / "risk/nonminimum_sensitivity.csv",
        "forensic": REVISION / "forensic/stable_sensitivity_estimates.csv",
        "tails": REVISION / "pxr/propensity_tails.csv",
        "weights": REVISION / "pxr/weighting_sensitivity.csv",
        "parameters": REVISION / "pxr/parameter_completeness.csv",
    }
    data = {name: pd.read_csv(path, float_precision="round_trip") for name, path in input_paths.items()}
    expected = {"estimators": 384, "candidate": 16, "calibration": 192, "nonminimum": 64, "forensic": 256, "tails": 18, "weights": 4, "parameters": 12}
    assert {name: len(frame) for name, frame in data.items()} == expected
    tables = {}
    descriptions = {}

    def add(name, frame, description):
        assert frame.columns.is_unique
        assert not {"molecule_id", "OCNT_ID", "SMILES", "source_row_index", "canonical_smiles"} & set(frame.columns)
        tables[name] = frame
        descriptions[name] = description

    estimator, candidate = data["estimators"], data["candidate"].sort_values(KEYS)
    add("Table_S8_original_estimator_summary.csv", estimator, "All 384 original estimator/configuration summaries across 100 masks; repeated CC/oracle controls remain indexed by original configurations.")
    add("Table_S8_candidate_gate_chain.csv", candidate, "The 16 original low-P/low-M candidate cells with selection strength, CC bias, and the original gate components.")
    mar = candidate[candidate.mechanism.eq("mar")]
    add("Table_S8A_MAR_bias_display.csv", readable(mar, {"endpoint": "Endpoint", "anchor": "Anchor", "model": "Model", "p_loss_correlation": "Corr(p,L)", "probability_weighted_shift": "Tp − T", "cc_signed_bias": "CC bias", "cc_rmse": "CC RMSE"}), "Eight original MAR cells. Bias is Monte Carlo signed error; Tp is the probability-weighted finite-cohort loss, not exact expected complete-case risk.")
    interval = estimator[estimator.propensity_representation.eq("lowdim_11") & estimator.loss_representation.eq("lowdim_11") & estimator.method.eq("fitted_aipw")].sort_values(KEYS)
    assert len(interval) == 16 and not interval.duplicated(KEYS).any()
    add("Table_S8B_candidate_intervals_display.csv", readable(interval, {"endpoint": "Endpoint", "anchor": "Anchor", "model": "Model", "mechanism": "Mechanism", "mean_se": "Mean SE", "mc_sd": "MC SD", "mean_interval_width": "Mean width", "coverage": "Coverage"}), "Sixteen original fitted low-P/low-M AIPW cells; original cluster SE and nominal 95% interval performance across 100 masks.")
    calibration = data["calibration"]
    add("Table_S9_mechanism_calibration.csv", calibration, "All 192 post-result fixed-nuisance calibration summaries; 10,000 Bernoulli masks per cell. Full original gate was not evaluated.")
    rows = []
    methods = {"true_p_fixed_replay_m_aipw": "True p; replay m", "constant_p_fixed_mean_aipw": "Constant p; mean m"}
    for mechanism, label in MECHANISMS.items():
        for method, method_label in methods.items():
            selected = calibration[calibration.mechanism.eq(mechanism) & calibration.method.eq(method)]
            assert len(selected) == 8
            rows.append({"Mechanism": label, "Reference": method_label, "Bias range": value_range(selected.signed_bias, 4), "Error-ratio range": value_range(selected.median_abs_error_ratio_vs_cc), "Coverage range": value_range(selected.coverage, 4), "Ratio passes": f"{int(selected.original_error_component_pass.sum())}/8"})
    add("Table_S9_calibration_display.csv", pd.DataFrame(rows), "Ranges across eight endpoint/anchor/model cells. Ratio passes evaluate only the original error component (MCAR <=1.25; other mechanisms <=0.75), not the whole original gate. Loss-rank mechanisms intentionally use known full-cohort loss and are diagnostic designs.")
    forensic = data["forensic"]
    add("Table_S10_forensic_estimates.csv", forensic, "All 256 single-seed (41,001), high-dimensional loss-regression forensic estimates, including numerical replay and controlled IPC clipping; oracle rows are repeated by propensity representation in the electronic source.")
    forensic_rows = []
    for probability, representation, label in [("fitted", "lowdim_11", "Low P"), ("fitted", "full_2266", "Full P"), ("oracle", "lowdim_11", "Oracle p")]:
        subset = forensic[forensic.probability.eq(probability) & forensic.propensity_representation.eq(representation)]
        stored = subset[subset.variant.eq("stored")].sort_values(KEYS)
        clipped = subset[subset.variant.eq("ipc_train_range_clip")].sort_values(KEYS)
        assert len(stored) == len(clipped) == 16 and not stored.duplicated(KEYS).any()
        assert stored[KEYS].reset_index(drop=True).equals(clipped[KEYS].reset_index(drop=True))
        if probability == "oracle":
            redundant = forensic[forensic.probability.eq("oracle") & forensic.propensity_representation.eq("full_2266")]
            for variant, reference in [("stored", stored), ("ipc_train_range_clip", clipped)]:
                other = redundant[redundant.variant.eq(variant)].sort_values(KEYS)
                assert np.allclose(reference.absolute_error.to_numpy(), other.absolute_error.to_numpy(), rtol=0, atol=1e-12)
        forensic_rows.append({"probability_specification": label, "n_distinct": len(stored), "stored_median": stored.absolute_error.median(), "clipped_median": clipped.absolute_error.median(), "clipped_max": clipped.absolute_error.max()})
    forensic_summary = pd.DataFrame(forensic_rows)
    add("Table_S10_forensic_display.csv", readable(forensic_summary, {"probability_specification": "Probability", "n_distinct": "N", "stored_median": "Stored median error", "clipped_median": "Clipped median error", "clipped_max": "Clipped maximum error"}), "Three probability specifications, each summarizing 16 distinct seed-41,001 estimates. Duplicate oracle representation rows are excluded. Clipping targets only IPC at fold-training ranges.")
    add("Table_S11A_propensity_tails.csv", data["tails"], "Stored PXR fitted-probability tails in all nine strata for two specifications; lower and upper tails have different missing-outcome interpretations.")
    add("Table_S11B_weighting_sensitivity.csv", data["weights"], "Frozen two-sided clipping and post-result lower-only clipping summaries; all fitted vectors unchanged.")
    add("Table_S11C_parameter_completeness.csv", data["parameters"], "Numerical parameter completeness in released PXR tables; this is not a predictor/loss/eligibility contract.")
    tails = data["tails"][data["tails"].stratum.eq("all_parent")].copy()
    pxr_labels = {"pxr_structure_10": "Structure only", "pxr_structure_screen_12": "Structure + screen"}
    tails["specification"] = tails.specification.replace(pxr_labels)
    tails["low_percent"] = tails.low_fraction * 100
    tails["high_percent"] = tails.high_fraction * 100
    add("Table_S11A_tails_display.csv", readable(tails, {"specification": "Specification", "low_count": "p < 0.05: N", "low_percent": "Lower tail (%)", "high_count": "p > 0.95: N", "high_percent": "Upper tail (%)"}, digits=2), "Tail denominator is all 10,747 PXR parent compounds for each specification.")
    weights = data["weights"][data["weights"].weight_rule.eq("post_result_lower_only_005")].copy()
    weights["specification"] = weights.specification.replace(pxr_labels)
    add("Table_S11B_lower_only_display.csv", readable(weights, {"specification": "Specification", "ess": "ESS", "ess_fraction": "ESS / observed", "max_abs_smd": "Maximum |SMD|"}), "Lower-only clipping at 0.05; observed denominator 2,722. These diagnostic quantities do not prove or disprove identification.")
    nonminimum = data["nonminimum"]
    add("Table_S12_nonminimum_sensitivity.csv", nonminimum, "All 64 above-minimum/all-measured CC and oracle-Hajek sensitivity summaries across the original 100 masks.")
    subset = nonminimum[nonminimum.subset.eq("strictly_above_released_minimum") & nonminimum.mechanism.eq("mar") & nonminimum.method.eq("cc")].sort_values(KEYS)
    assert len(subset) == 8
    add("Table_S12_nonminimum_display.csv", readable(subset, {"endpoint": "Endpoint", "anchor": "Anchor", "model": "Model", "n_total": "N", "truth": "True risk", "p_loss_correlation": "Corr(p,L)", "probability_weighted_shift": "Tp − T", "signed_bias": "CC bias", "rmse": "CC RMSE"}), "Eight original MAR CC cells restricted to released values strictly above the cohort minimum; the target cohort and truth change, and the selection mechanism is not recalibrated.")
    for name, frame in tables.items():
        base.write_csv(frame, TABLES / name)
    plot_calibration = calibration[calibration.method.eq("true_p_fixed_replay_m_aipw")].copy()
    base.write_csv(plot_calibration, SOURCE / "figure4_fixed_nuisance_calibration.csv")
    base.write_csv(forensic_summary, SOURCE / "figure4_ipc_sensitivity.csv")
    figure_paths = base.save_figure(figure4(plot_calibration, forensic_summary), "Figure_4_diagnostic_mechanisms")
    report = {"build": "revision_001_tables_S8_S12_figure4", "source_sha256": {path.relative_to(ROOT).as_posix(): base.sha256(path) for path in input_paths.values()}, "script_sha256": base.sha256(Path(__file__)), "tables": {name: {"rows": len(frame), "columns": frame.columns.tolist(), "description": descriptions[name], "sha256": base.sha256(TABLES / name)} for name, frame in tables.items()}, "figure_sources": {name: base.sha256(SOURCE / name) for name in ["figure4_fixed_nuisance_calibration.csv", "figure4_ipc_sensitivity.csv"]}, "figure_exports": {path.relative_to(ROOT).as_posix(): base.sha256(path) for path in figure_paths}, "figure_surface": "paper_main", "visual_review": "Separate rendered-image inspection required after generation.", "preservation": "Original builders and accepted confirmatory outputs are not modified by this script."}
    (QA / "revision_assets_build_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (QA / "revision_assets_build_report.md").write_text("# Revision figure and table build\n\n" + "\n\n".join(f"- `{name}`: {len(frame)} rows. {descriptions[name]}" for name, frame in tables.items()) + "\n\nFigure 4 uses post-result diagnostics only. Panel a: eight fixed-predictor cells per selection mechanism. Panel b: medians of 16 distinct estimates for one replay seed. Panel c: eight paired SE/MC-SD values under original MAR with fixed nuisance predictions.\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "tables": len(tables), "figure_exports": len(figure_paths), "display_tables": [name for name in tables if "display" in name]}, indent=2))


if __name__ == "__main__":
    main()
