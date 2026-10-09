import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


BASE = Path(__file__).resolve().parents[1]
FROZEN = BASE / "artifacts/confirmatory/confirm_001/full_recovery_001"
PAPER = BASE / "paper"
FIGURES = PAPER / "figures"
SOURCE = PAPER / "source_data"
TABLES = PAPER / "tables"
QA = PAPER / "qa"

EXPECTED_FILES = {
    "actual_mask_scenarios.csv",
    "admission_gate.json",
    "availability_map.csv",
    "bootstrap_intervals.csv",
    "cohort_audit.csv",
    "factorial_contrasts.csv",
    "generator_vectors.parquet",
    "nuisance_fold_assignments.parquet",
    "predictions.parquet",
    "pxr_boundary_decision.json",
    "pxr_linkage_audit.csv",
    "pxr_overlap_summary.csv",
    "pxr_selection_strata.csv",
    "row_replay_seed41001.parquet",
    "run_manifest.json",
    "semisynthetic_scenarios.csv",
    "structure_audit.json",
    "summary.json",
    "verification.md",
}

COLORS = {
    "navy": "#23395B",
    "teal": "#2A7F8E",
    "amber": "#D39424",
    "red": "#B64B4B",
    "green": "#3D7D69",
    "ink": "#20262E",
    "grey": "#6B7280",
    "light": "#F4F6F8",
    "pale_teal": "#E5F2F3",
    "pale_amber": "#FAF1DD",
    "pale_red": "#F7E8E8",
    "white": "#FFFFFF",
}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(df, path):
    df.to_csv(path, index=False, float_format="%.17g")


def assert_close(actual, expected, atol=1e-12):
    assert np.isfinite(actual)
    assert abs(float(actual) - float(expected)) <= atol


def style():
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.7,
            "svg.fonttype": "none",
            "svg.hashsalt": "jpa-selective-reference-20260911",
            "pdf.fonttype": 42,
            "pdf.use14corefonts": False,
            "savefig.facecolor": "white",
        }
    )


def save_figure(fig, stem):
    paths = []
    for suffix, dpi in [("svg", 300), ("pdf", 300), ("png", 300), ("tiff", 600)]:
        path = FIGURES / f"{stem}.{suffix}"
        kwargs = {"dpi": dpi, "bbox_inches": "tight", "pad_inches": 0.04}
        if suffix == "tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        elif suffix == "svg":
            kwargs["metadata"] = {"Date": None, "Creator": "JPA reproducibility package"}
        elif suffix == "pdf":
            kwargs["metadata"] = {"Creator": "JPA reproducibility package", "CreationDate": None, "ModDate": None}
        fig.savefig(path, **kwargs)
        paths.append(path)
    plt.close(fig)
    return paths


def box(ax, x, y, w, h, title, lines, face, edge, title_color=None):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=1,
        facecolor=face,
        edgecolor=edge,
    )
    ax.add_patch(patch)
    ax.text(
        x + 0.03 * w,
        y + h - 0.15 * h,
        title,
        fontsize=8.8,
        fontweight="bold",
        color=title_color or edge,
        va="top",
    )
    ax.text(
        x + 0.03 * w,
        y + h - 0.35 * h,
        "\n".join(lines),
        fontsize=7.1,
        color=COLORS["ink"],
        va="top",
        linespacing=1.25,
    )


def arrow(ax, x1, y1, x2, y2, color=None):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=10,
            linewidth=1.1,
            color=color or COLORS["grey"],
        )
    )


def figure1(counts):
    fig, ax = plt.subplots(figsize=(7.2, 4.15))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.01, 0.98, "a", fontsize=10, fontweight="bold", va="top")
    ax.text(0.055, 0.975, "Non-exchangeable evidence roles", fontsize=9, fontweight="bold", va="top")
    box(
        ax,
        0.055,
        0.63,
        0.25,
        0.25,
        "Development",
        ["ExpansionRx", "Hypothesis generation only", "0 rows in confirmation"],
        COLORS["light"],
        COLORS["grey"],
    )
    box(
        ax,
        0.37,
        0.63,
        0.25,
        0.25,
        "Confirmation",
        ["Biogen public ADME", "3,521 structures", "2,463 train | 1,058 test"],
        COLORS["pale_teal"],
        COLORS["teal"],
    )
    box(
        ax,
        0.685,
        0.63,
        0.25,
        0.25,
        "Boundary case",
        ["OpenADMET PXR", "10,747 parent", "2,722 follow-up"],
        COLORS["pale_amber"],
        COLORS["amber"],
    )
    arrow(ax, 0.305, 0.755, 0.37, 0.755)
    ax.plot([0.652, 0.652], [0.64, 0.87], color=COLORS["grey"], lw=0.8, ls=(0, (2, 2)))
    ax.text(0.655, 0.9, "not a validation cohort", fontsize=6.5, color=COLORS["grey"], ha="center")

    ax.text(0.01, 0.56, "b", fontsize=10, fontweight="bold", va="top")
    ax.text(0.055, 0.555, "Frozen confirmation experiment", fontsize=9, fontweight="bold", va="top")
    box(
        ax,
        0.055,
        0.27,
        0.19,
        0.20,
        "Fixed losses",
        ["HLM: n=929", "RLM: n=922", "RF400 + XGB800"],
        COLORS["light"],
        COLORS["navy"],
    )
    box(
        ax,
        0.285,
        0.27,
        0.20,
        0.20,
        "Known truth",
        ["MCAR and benign MAR", "SOL + MDR1 anchors", "100 frozen masks"],
        COLORS["light"],
        COLORS["navy"],
    )
    box(
        ax,
        0.525,
        0.27,
        0.19,
        0.20,
        "2 × 2 nuisance",
        ["P: 11 or 2,266", "M: 11 or 2,266", "6,400 result rows"],
        COLORS["light"],
        COLORS["navy"],
    )
    box(
        ax,
        0.755,
        0.27,
        0.19,
        0.20,
        "Prespecified gates",
        ["H1/H2 confirmation", "16-cell admission", "direct + PXR rules"],
        COLORS["light"],
        COLORS["navy"],
    )
    arrow(ax, 0.245, 0.37, 0.285, 0.37)
    arrow(ax, 0.485, 0.37, 0.525, 0.37)
    arrow(ax, 0.715, 0.37, 0.755, 0.37)

    ax.text(0.01, 0.20, "c", fontsize=10, fontweight="bold", va="top")
    ax.text(0.055, 0.195, "Decision chain", fontsize=9, fontweight="bold", va="top")
    outcomes = [
        (0.055, 0.085, 0.18, "Fixed comparisons", "H1 4/4 | H2 4/4", COLORS["pale_teal"], COLORS["teal"]),
        (0.285, 0.085, 0.18, "Relative benefit", "MAR criterion: 0/8", COLORS["pale_red"], COLORS["red"]),
        (0.515, 0.085, 0.18, "Direct masks", "8/8; not adjudicated", COLORS["pale_amber"], COLORS["amber"]),
        (0.745, 0.085, 0.20, "PXR records", "Full risk not estimated", COLORS["pale_amber"], COLORS["amber"]),
    ]
    for x, y, w, title, text, face, edge in outcomes:
        p = FancyBboxPatch((x, y), w, 0.065, boxstyle="round,pad=0.008", facecolor=face, edgecolor=edge, lw=1)
        ax.add_patch(p)
        ax.text(x + 0.012, y + 0.045, title, fontsize=7.4, fontweight="bold", color=edge, va="center")
        ax.text(x + 0.012, y + 0.016, text, fontsize=6.4, color=COLORS["ink"], va="center")
    arrow(ax, 0.235, 0.118, 0.285, 0.118)
    arrow(ax, 0.465, 0.118, 0.515, 0.118)
    arrow(ax, 0.695, 0.118, 0.745, 0.118)
    ax.text(0.5, 0.012, "Post-result diagnostics distinguish relative benefit, absolute error and uncertainty.", ha="center", fontsize=8, fontweight="bold", color=COLORS["navy"])
    return fig


def clean_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#D8DDE3", lw=0.6, zorder=0)
    ax.tick_params(length=3, width=0.6)


def figure2(boot):
    stat_order = [
        ("oracle_to_fitted_aipw_error_ratio", "a", "Oracle / fitted-full error", "ratio (lower is better)", 0.5),
        ("low_to_full_fitted_aipw_error_ratio", "b", "Low-P / full-P error", "ratio (lower is better)", 0.5),
        ("full_minus_low_clip_fraction", "c", "Full-P minus low-P clipping", "absolute fraction (higher is better)", 0.2),
    ]
    labels = ["HLM | MCAR", "HLM | MAR", "RLM | MCAR", "RLM | MAR"]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 4.1), sharey=True)
    for ax, (stat, panel, title, xlabel, threshold) in zip(axes, stat_order):
        d = boot[boot["statistic"] == stat].copy()
        d["key"] = d["endpoint"].map({"LOG_HLM": "HLM", "LOG_RLM": "RLM"}) + " | " + d["mechanism"].str.upper()
        d = d.set_index("key").loc[labels].reset_index()
        y = np.arange(len(d))[::-1]
        x = d["point_estimate"].to_numpy()
        xerr = np.vstack([x - d["lower"].to_numpy(), d["upper"].to_numpy() - x])
        ax.errorbar(x, y, xerr=xerr, fmt="o", color=COLORS["navy"], ecolor=COLORS["teal"], capsize=2.5, ms=4.5, lw=1.2, zorder=3)
        ax.axvline(threshold, color=COLORS["red"], lw=1.1, ls=(0, (4, 2)))
        ax.set_title(title, loc="left", fontweight="bold", pad=8)
        ax.text(-0.12, 1.08, panel, transform=ax.transAxes, fontsize=10, fontweight="bold")
        ax.set_xlabel(xlabel)
        ax.set_yticks(y)
        if ax is axes[0]:
            ax.set_yticklabels(labels)
        else:
            ax.tick_params(labelleft=False)
        clean_axis(ax)
        xmax = max(float(d["upper"].max()) * 1.18, threshold * 1.18)
        ax.set_xlim(0, xmax)
        ax.text(threshold, -0.72, f"threshold {threshold:.2f}", ha="center", color=COLORS["red"], fontsize=6.5)
    fig.suptitle("Prespecified failure signatures reproduced in all endpoint–mechanism cells", y=1.015, fontsize=10.5, fontweight="bold", color=COLORS["ink"])
    fig.text(0.5, -0.02, "Points are frozen estimates; bars are 10,000 whole-seed bootstrap intervals. Thresholds, not significance tests, adjudicated H1/H2.", ha="center", fontsize=7, color=COLORS["grey"])
    fig.subplots_adjust(wspace=0.18, left=0.17, right=0.99, top=0.82, bottom=0.21)
    return fig


def figure3(admission, direct, strata, overlap):
    fig = plt.figure(figsize=(7.2, 5.18))
    outer = fig.add_gridspec(2, 2, height_ratios=[1.05, 1], width_ratios=[1.06, 0.94], hspace=0.78, wspace=0.34)
    axa = fig.add_subplot(outer[:, 0])
    axb = fig.add_subplot(outer[0, 1])
    cgrid = outer[1, 1].subgridspec(1, 2, width_ratios=[0.92, 1.08], wspace=0.42)
    axc1 = fig.add_subplot(cgrid[0, 0])
    axc2 = fig.add_subplot(cgrid[0, 1])

    adm = admission.copy()
    adm["row"] = adm["endpoint"].map({"LOG_HLM": "HLM", "LOG_RLM": "RLM"}) + "–" + adm["anchor"].map({"LOG_SOL": "SOL", "LOG_MDR1_ER": "MDR1"}) + "–" + adm["model"].map({"ecfp_rf": "RF", "ecfp_desc_xgb": "XGB"})
    order = [f"{e}–{a}–{m}" for e in ["HLM", "RLM"] for a in ["SOL", "MDR1"] for m in ["RF", "XGB"]]
    matrix = np.zeros((len(order), 2))
    notes = np.empty((len(order), 2), dtype=object)
    for i, row in enumerate(order):
        for j, mech in enumerate(["mcar", "mar"]):
            r = adm[(adm["row"] == row) & (adm["mechanism"] == mech)].iloc[0]
            matrix[i, j] = 1 if bool(r["pass"]) else 0
            notes[i, j] = "PASS" if bool(r["pass"]) else str(r["failure_reason"])
    cmap = mpl.colors.ListedColormap([COLORS["pale_red"], COLORS["pale_teal"]])
    axa.imshow(matrix, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    axa.set_xticks([0, 1], ["MCAR", "MAR"])
    axa.set_yticks(np.arange(len(order)), order)
    for i in range(len(order)):
        for j in range(2):
            axa.text(j, i, notes[i, j], ha="center", va="center", fontsize=6.4, color=COLORS["green"] if matrix[i, j] else COLORS["red"], fontweight="bold")
    axa.set_title("a  Known-truth admission: 8/16", loc="left", fontweight="bold", pad=25)
    axa.text(0.5, -0.09, "All eight MAR cells failed the error-ratio threshold; six also failed SMD.", transform=axa.transAxes, ha="center", fontsize=6.7, color=COLORS["red"])
    for spine in axa.spines.values():
        spine.set_visible(False)
    axa.set_xticks(np.arange(-0.5, 2, 1), minor=True)
    axa.set_yticks(np.arange(-0.5, len(order), 1), minor=True)
    axa.grid(which="minor", color="white", linewidth=2)
    axa.tick_params(which="minor", bottom=False, left=False)

    direct = direct.copy()
    direct["label"] = direct["endpoint"].map({"LOG_HLM": "HLM", "LOG_RLM": "RLM"}) + "–" + direct["anchor"].map({"LOG_SOL": "SOL", "LOG_MDR1_ER": "MDR1"}) + "–" + direct["model"].map({"ecfp_rf": "RF", "ecfp_desc_xgb": "XGB"})
    direct = direct.sort_values(["endpoint", "anchor", "model"], ascending=[True, True, True]).reset_index(drop=True)
    y = np.arange(len(direct))[::-1]
    for i, r in direct.iterrows():
        yy = y[i]
        axb.plot([r["cc_abs_error"], r["fitted_aipw_abs_error"]], [yy, yy], color="#B9C0C8", lw=1)
    axb.scatter(direct["cc_abs_error"], y, marker="o", s=22, color=COLORS["navy"], label="Complete case", zorder=3)
    axb.scatter(direct["fitted_aipw_abs_error"], y, marker="s", s=22, color=COLORS["amber"], label="Candidate AIPW", zorder=3)
    axb.set_xscale("log")
    axb.set_yticks(y, direct["label"])
    axb.set_xlabel("absolute finite-population risk error")
    axb.set_title("b  Direct masks: overlap 8/8", loc="left", fontweight="bold", pad=28)
    axb.legend(frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(0.5, 1.005), fontsize=6.4, handletextpad=0.4, columnspacing=0.8)
    axb.text(0.5, -0.38, "AIPW lower in 4/8; pooled reduction −34.1%\nNOT ADJUDICATED after admission failure", transform=axb.transAxes, ha="center", fontsize=6.5, color=COLORS["red"], fontweight="bold")
    clean_axis(axb)

    s = strata.sort_values("formal_hit", ascending=False).reset_index(drop=True)
    names = ["Formal hit", "Non-hit"]
    available = s["n_followup_available"].to_numpy()
    unavailable = s["n_followup_unavailable"].to_numpy()
    xx = np.arange(2)
    axc1.bar(xx, available, color=COLORS["teal"], width=0.65, label="Follow-up")
    axc1.bar(xx, unavailable, bottom=available, color="#DDE1E6", width=0.65, label="Unavailable")
    axc1.set_xticks(xx, names, rotation=22, ha="right")
    axc1.set_ylabel("parent compounds")
    axc1.set_title("c  PXR selection", loc="left", fontweight="bold", pad=7)
    axc1.legend(frameon=False, fontsize=6.2, loc="upper left")
    for i, r in s.iterrows():
        axc1.text(i, r["n_total"] * 1.02, f"{r['advancement_probability']:.1%}", ha="center", fontsize=6.5, color=COLORS["ink"])
    clean_axis(axc1)

    metrics = [("clip_fraction", 0.10, "Two-sided clip"), ("ess_fraction", 0.50, "ESS/obs"), ("max_abs_smd", 0.10, "max |SMD|")]
    specs = ["pxr_structure_10", "pxr_structure_screen_12"]
    x2 = np.arange(len(metrics))
    widths = [-0.16, 0.16]
    markers = ["o", "s"]
    colors = [COLORS["teal"], COLORS["red"]]
    for spec, dx, marker, color in zip(specs, widths, markers, colors):
        r = overlap[overlap["specification"] == spec].iloc[0]
        vals = [r[m] / t for m, t, _ in metrics]
        axc2.scatter(x2 + dx, vals, marker=marker, s=26, color=color, label="Structure" if spec.endswith("10") else "Structure + screen", zorder=3)
    axc2.axhline(1, color=COLORS["grey"], ls=(0, (3, 2)), lw=1)
    axc2.set_xticks(x2, [m[2] for m in metrics], rotation=22, ha="right")
    axc2.set_ylabel("diagnostic / threshold")
    axc2.set_yscale("log")
    axc2.set_ylim(0.25, 5.5)
    axc2.legend(frameon=False, fontsize=5.9, loc="upper left")
    clean_axis(axc2)
    axc2.text(0.5, -0.37, "PXR risk not estimated\nAssumptions and loss not established", transform=axc2.transAxes, ha="center", fontsize=6.7, color=COLORS["amber"], fontweight="bold")

    fig.subplots_adjust(left=0.13, right=0.99, top=0.92, bottom=0.14)
    return fig


def heatmap(ax, values, annotations, title, cmap, vmin=None, vmax=None):
    im = ax.imshow(values, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(4), ["low P\nlow M", "full P\nlow M", "low P\nfull M", "full P\nfull M"])
    ax.set_yticks(np.arange(4), ["HLM | MCAR", "HLM | MAR", "RLM | MCAR", "RLM | MAR"])
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            rgba = im.cmap(im.norm(values[i, j]))
            luminance = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
            color = COLORS["ink"] if luminance > 0.58 else "white"
            ax.text(j, i, annotations[i, j], ha="center", va="center", fontsize=6.1, color=color)
    ax.set_title(title, loc="left", fontweight="bold", pad=7)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, 4, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 4, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    return im


def figure_s1(factorial):
    rows = [("LOG_HLM", "mcar"), ("LOG_HLM", "mar"), ("LOG_RLM", "mcar"), ("LOG_RLM", "mar")]
    cols = [("lowdim_11", "lowdim_11"), ("full_2266", "lowdim_11"), ("lowdim_11", "full_2266"), ("full_2266", "full_2266")]

    def matrix(column):
        a = np.zeros((4, 4))
        for i, (endpoint, mechanism) in enumerate(rows):
            for j, (p, m) in enumerate(cols):
                a[i, j] = factorial[(factorial["endpoint"] == endpoint) & (factorial["mechanism"] == mechanism) & (factorial["propensity_representation"] == p) & (factorial["loss_representation"] == m)][column].iloc[0]
        return a

    err = matrix("median_abs_error")
    coverage = matrix("coverage")
    clip = matrix("median_clip")
    ess = matrix("median_ess_fraction")
    fig, axes = plt.subplots(2, 2, figsize=(5.9, 4.1))
    heatmap(axes[0, 0], np.log10(err), np.vectorize(lambda x: f"{x:.2g}")(err), "a  Median absolute estimation error", "magma")
    heatmap(axes[0, 1], coverage, np.vectorize(lambda x: f"{x:.3f}")(coverage), "b  Empirical 95% interval coverage", "YlGnBu", 0.90, 1.0)
    heatmap(axes[1, 0], clip, np.vectorize(lambda x: f"{x:.3f}")(clip), "c  Median propensity clipping", "YlOrRd", 0, 0.65)
    heatmap(axes[1, 1], ess, np.vectorize(lambda x: f"{x:.3f}")(ess), "d  Median ESS / observed", "YlGnBu", 0.25, 1.0)
    fig.suptitle("Complete propensity-by-loss nuisance factorial", y=1.015, fontsize=10.2, fontweight="bold")
    fig.text(0.5, -0.015, "Full-loss cells have errors of 165–485; their nominal coverage reflects extremely wide intervals, not superior recovery.", ha="center", fontsize=7, color=COLORS["red"], fontweight="bold")
    fig.subplots_adjust(left=0.17, right=0.99, top=0.88, bottom=0.16, hspace=0.46, wspace=0.32)
    return fig


def graphical_abstract():
    fig, ax = plt.subplots(figsize=(9.0, 3.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.5, 0.955, "Interpreting selective-reference evaluation of clearance models", ha="center", va="top", fontsize=14, fontweight="bold", color=COLORS["navy"])
    ax.text(0.5, 0.875, "Measured Biogen cohorts: HLM n=929 and RLM n=922 | Fixed RF and XGBoost models", ha="center", va="top", fontsize=8.2, color=COLORS["grey"])

    panels = [
        (0.03, "1  Measure selection bias", ["Original MAR: weak bias", "|Corr(probability, loss)| < 0.107", "Oracle also misses 25% benefit", "Benefit is distinct from accuracy"], COLORS["pale_teal"], COLORS["teal"]),
        (0.37, "2  Trace extrapolation", ["Single-seed Ipc diagnostic", "Limit only held-out Ipc range", "Low-P median error: 702.65 → 0.00595", "Configuration-specific finding"], COLORS["pale_red"], COLORS["red"]),
        (0.71, "3  Audit PXR records", ["Dose-table membership", "Hits: 869/873", "Non-hits: 1,853/9,874", "Full-library risk not estimated"], COLORS["pale_amber"], COLORS["amber"]),
    ]
    for x, title, lines, face, edge in panels:
        patch = FancyBboxPatch((x, 0.25), 0.26, 0.54, boxstyle="round,pad=0.012,rounding_size=0.025", linewidth=1.6, facecolor=face, edgecolor=edge)
        ax.add_patch(patch)
        ax.text(x + 0.02, 0.735, title, fontsize=10.2, fontweight="bold", color=edge, va="top")
        for i, line in enumerate(lines):
            weight = "bold" if i == len(lines) - 1 else "normal"
            color = edge if i == len(lines) - 1 else COLORS["ink"]
            ax.text(x + 0.015, 0.64 - i * 0.105, line, fontsize=7.3, fontweight=weight, color=color, va="top")
    arrow(ax, 0.30, 0.52, 0.36, 0.52, COLORS["grey"])
    arrow(ax, 0.64, 0.52, 0.70, 0.52, COLORS["grey"])
    ax.text(0.50, 0.205, "Original frozen results retained; bias, interval and extrapolation diagnostics are post-result", ha="center", fontsize=8.3, color=COLORS["navy"])
    banner = FancyBboxPatch((0.13, 0.045), 0.74, 0.095, boxstyle="round,pad=0.01,rounding_size=0.025", linewidth=0, facecolor=COLORS["navy"])
    ax.add_patch(banner)
    ax.text(0.50, 0.092, "Check bias strength, feature extrapolation and the uncertainty target separately", ha="center", va="center", fontsize=9.2, fontweight="bold", color=COLORS["white"])
    fig.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.01)
    return fig


def main():
    for path in [FIGURES, SOURCE, TABLES, QA]:
        path.mkdir(parents=True, exist_ok=True)
    actual = {p.name for p in FROZEN.iterdir() if p.is_file()}
    assert actual == EXPECTED_FILES
    assert not any(p.is_dir() for p in FROZEN.iterdir())
    fixed_hashes = {
        "run_manifest.json": "8da9b9c1d0425fff118d983a8491aceb57ae8b6085300eacf18cd8e0d7e6eac4",
        "summary.json": "784132d1df5ab6a99a934bf2c53d17c7d11f650beb4b733be222011f060986ae",
        "verification.md": "7ef5b2cdb2d1885b074250402efd313aa9758b337b4f486147ddc10e0042af99",
        "generator_vectors.parquet": "a6b817b9a3db207dd84f1a816a7374a5a21b186c853b73c2f2c422a1bfc830f7",
    }
    assert {name: sha256(FROZEN / name) for name in fixed_hashes} == fixed_hashes
    manifest = json.loads((FROZEN / "run_manifest.json").read_text(encoding="utf-8"))
    for rel, expected in manifest["output_hashes_before_manifest"].items():
        target = FROZEN / rel
        assert sha256(target) == expected

    summary = json.loads((FROZEN / "summary.json").read_text(encoding="utf-8"))
    admission_json = json.loads((FROZEN / "admission_gate.json").read_text(encoding="utf-8"))
    boundary = json.loads((FROZEN / "pxr_boundary_decision.json").read_text(encoding="utf-8"))
    structure = json.loads((FROZEN / "structure_audit.json").read_text(encoding="utf-8"))
    boot = pd.read_csv(FROZEN / "bootstrap_intervals.csv")
    cohort = pd.read_csv(FROZEN / "cohort_audit.csv")
    availability = pd.read_csv(FROZEN / "availability_map.csv")
    scenarios = pd.read_csv(FROZEN / "semisynthetic_scenarios.csv")
    direct_all = pd.read_csv(FROZEN / "actual_mask_scenarios.csv")
    contrasts = pd.read_csv(FROZEN / "factorial_contrasts.csv")
    strata = pd.read_csv(FROZEN / "pxr_selection_strata.csv")
    overlap = pd.read_csv(FROZEN / "pxr_overlap_summary.csv")
    predictions = pd.read_parquet(FROZEN / "predictions.parquet")

    assert summary["integrity"]["pass"] is True
    assert summary["decision"] == "GO_WRITE_FALSIFICATION_BOUNDARY"
    assert len(boot) == 12 and set(boot["resamples"]) == {10000}
    assert boot.groupby("statistic").size().to_dict() == {
        "full_minus_low_clip_fraction": 4,
        "low_to_full_fitted_aipw_error_ratio": 4,
        "oracle_to_fitted_aipw_error_ratio": 4,
    }
    assert boot["seed_index_matrix_sha256"].nunique() == 1

    admission = pd.DataFrame(admission_json["known_truth_admission"]["cells"])
    def failure_reason(r):
        failures = []
        if r["error_ratio"] > r["error_threshold"]:
            failures.append("ratio")
        if r["coverage"] < 0.90:
            failures.append("coverage")
        if r["ess_fraction_q05"] < 0.50:
            failures.append("ESS")
        if r["clip_fraction_median"] > 0.10 or r["clip_fraction_q95"] > 0.20:
            failures.append("clipping")
        if r["max_abs_smd_q95"] > 0.10:
            failures.append("SMD")
        if not r["fold_class_support"]:
            failures.append("fold")
        if not r["all_finite"]:
            failures.append("finite")
        return "+".join(failures) if failures else "PASS"
    admission["failure_reason"] = admission.apply(failure_reason, axis=1)
    assert len(admission) == 16
    assert admission.groupby("mechanism")["pass"].sum().to_dict() == {"mar": 0, "mcar": 8}
    assert admission["failure_reason"].str.contains("ratio").sum() == 8
    assert admission["failure_reason"].str.contains("SMD").sum() == 6

    direct = direct_all[(direct_all["propensity_representation"] == "lowdim_11") & (direct_all["loss_representation"] == "lowdim_11") & (direct_all["anchor"].isin(["LOG_SOL", "LOG_MDR1_ER"]))].copy()
    direct["aipw_better"] = direct["fitted_aipw_abs_error"] < direct["cc_abs_error"]
    assert len(direct_all) == 64 and len(direct) == 8
    assert direct["overlap_eligible"].all() and int(direct["aipw_better"].sum()) == 4
    assert_close(direct["cc_abs_error"].median(), 0.00459157014222647)
    assert_close(direct["fitted_aipw_abs_error"].median(), 0.006157635158421465)
    assert admission_json["direct_mask_translation"]["status"] == "NOT_ADJUDICATED_ADMISSION_FAILED"

    factorial = scenarios.groupby(["endpoint", "mechanism", "propensity_representation", "loss_representation"], sort=False).agg(
        rows=("seed", "size"),
        median_abs_error=("fitted_aipw_abs_error", "median"),
        coverage=("fitted_covered", "mean"),
        median_clip=("fitted_clip_fraction", "median"),
        median_ess_fraction=("fitted_ess_fraction", "median"),
    ).reset_index()
    assert len(scenarios) == 6400 and len(factorial) == 16 and set(factorial["rows"]) == {400}
    assert len(contrasts) == 56 and contrasts["status"].str.startswith("exploratory_only").all()
    assert len(strata) == 2 and int(strata["n_total"].sum()) == 10747 and int(strata["n_followup_available"].sum()) == 2722
    assert len(overlap) == 2
    assert overlap.loc[overlap["specification"] == "pxr_structure_screen_12", ["trigger_clip", "trigger_ess", "trigger_smd"]].all(axis=None)
    assert boundary["decision"] == "ABSTAIN_NONIDENTIFIABLE_FULL_LIBRARY_DOSE_RESPONSE"
    assert boundary["full_library_point_estimate_allowed"] is False
    assert len(predictions) == 3702
    assert np.max(np.abs(predictions["loss"] - np.abs(predictions["reference"] - predictions["prediction"]))) <= 1e-12

    counts = pd.DataFrame(
        [
            ["development", "ExpansionRx contribution to confirmation", 0, np.nan, 0, "rows", "summary.json", "forensic_role", "hypothesis generation only"],
            ["confirmation", "Biogen structures", 3521, 3521, 1, "fraction", "structure_audit.json", "biogen", "unique canonical structures"],
            ["confirmation", "Biogen derived train split", 2463, 3521, 2463 / 3521, "fraction", "cohort_audit.csv", "train", "not temporal or prospective"],
            ["confirmation", "Biogen derived test split", 1058, 3521, 1058 / 3521, "fraction", "cohort_audit.csv", "test", "not an official Biogen split"],
            ["confirmation", "HLM target cohort", 929, 1058, 929 / 1058, "fraction", "cohort_audit.csv", "LOG_HLM test", "finite released reference"],
            ["confirmation", "RLM target cohort", 922, 1058, 922 / 1058, "fraction", "cohort_audit.csv", "LOG_RLM test", "finite released reference"],
            ["experiment", "Semisynthetic result rows", 6400, 6400, 1, "fraction", "run_manifest.json", "grid_rows", "complete frozen grid"],
            ["decision", "H1 cells passed", 4, 4, 1, "fraction", "summary.json", "H1", "mechanism confirmation"],
            ["decision", "H2 cells passed", 4, 4, 1, "fraction", "summary.json", "H2", "mechanism confirmation"],
            ["decision", "Admission cells passed", 8, 16, 0.5, "fraction", "admission_gate.json", "known_truth_admission", "all MAR cells failed"],
            ["decision", "Direct regular overlap cells passed", 8, 8, 1, "fraction", "admission_gate.json", "direct_overlap", "translation not adjudicated"],
            ["boundary", "PXR follow-up available", 2722, 10747, 2722 / 10747, "fraction", "pxr_selection_strata.csv", "sum", "full-library risk not estimated; assumptions and loss unestablished"],
        ],
        columns=["stage", "item", "numerator", "denominator", "value", "unit", "source_file", "source_field", "interpretation"],
    )

    risk = predictions.groupby(["endpoint", "model"], as_index=False).agg(n=("loss", "size"), mean_absolute_risk=("loss", "mean"), median_absolute_loss=("loss", "median"), rmse=("loss", lambda x: float(np.sqrt(np.mean(np.square(x))))))
    table1_rows = []
    table1_rows.extend(
        [
            ["Biogen cohort", "All structures", "n", 3521, "3,521 unique canonical structures"],
            ["Biogen cohort", "Derived train split", "n", 2463, "structure-grouped; not official or prospective"],
            ["Biogen cohort", "Derived test split", "n", 1058, "zero canonical overlap with train"],
        ]
    )
    for _, r in risk.iterrows():
        endpoint = r["endpoint"].replace("LOG_", "")
        model = {"ecfp_rf": "RF400", "ecfp_desc_xgb": "XGB800"}[r["model"]]
        table1_rows.append(["Fixed predictor risk", f"{endpoint} | {model}", "mean absolute risk", r["mean_absolute_risk"], f"target n={int(r['n'])}"])
    for _, r in availability.iterrows():
        table1_rows.append(["Anchor availability", f"{r['endpoint'].replace('LOG_', '')} | {r['anchor'].replace('LOG_', '').replace('_ER', '')}", "available/total", r["coverage"], f"{int(r['n_available'])}/{int(r['n_total'])}"])
    table1_rows.extend(
        [
            ["PXR boundary", "Parent library", "n", 10747, "public-file-specific parent"],
            ["PXR boundary", "Dose-response follow-up", "fraction", 2722 / 10747, "2,722/10,747"],
        ]
    )
    table1 = pd.DataFrame(table1_rows, columns=["section", "item", "measure", "value", "display_note"])

    table2 = pd.DataFrame(
        [
            ["H1", "oracle/full-fitted error ratio", "PASS", "4/4 cells; threshold <=0.50"],
            ["H2", "low/full propensity error and clipping", "PASS", "4/4 cells; thresholds <=0.50 and >=0.20"],
            ["Original admission", "lowdim_11/lowdim_11", "Not met", "MCAR 8/8; MAR 0/8; relative-benefit rule"],
            ["Direct overlap", "regular empirical masks", "PASS", "8/8 cells"],
            ["Direct translation", "error improvement", "Not adjudicated", "original rule; descriptive improvement 4/8"],
            ["PXR", "selection-record audit", "Not estimated", "identification assumptions and endpoint-loss contract unestablished"],
        ],
        columns=["analysis", "estimand_or_test", "status", "evidence"],
    )

    pxr_s5 = pd.concat(
        [
            strata.assign(record_type="selection_stratum", entity=strata["formal_hit"].map({1: "formal_hit", 0: "non_hit"})),
            overlap.assign(record_type="overlap_specification", entity=overlap["specification"]),
        ],
        ignore_index=True,
        sort=False,
    )
    integrity = pd.DataFrame(
        [
            ["immutable inventory", "19 files; 0 directories", "PASS"],
            ["frozen inputs", "51/51 hashes", "PASS"],
            ["pre-manifest outputs", "18/18 hashes", "PASS"],
            ["row replay", "59,232 rows at seed 41,001", "PASS"],
            ["formula replay", "maximum absolute difference 0; tolerance 1e-12", "PASS"],
            ["external verification", "verification replay 002", "PASS"],
            ["verification scope", "file identity and computational replay; not statistical validity", "Qualified"],
        ],
        columns=["check", "result", "status"],
    )

    source_outputs = {
        SOURCE / "figure1_design_counts.csv": counts,
        SOURCE / "figure2_confirmation.csv": boot,
        SOURCE / "figure3_admission.csv": admission,
        SOURCE / "figure3_direct_masks.csv": direct,
        SOURCE / "figure3_pxr_strata.csv": strata,
        SOURCE / "figure3_pxr_overlap.csv": overlap,
        SOURCE / "figureS1_factorial_summary.csv": factorial,
        TABLES / "Table_1_cohorts_models_and_availability.csv": table1,
        TABLES / "Table_2_gate_summary.csv": table2,
        TABLES / "Table_S1_factorial_summary.csv": factorial,
        TABLES / "Table_S2_admission_matrix.csv": admission,
        TABLES / "Table_S3_direct_masks_all_cells.csv": direct_all,
        TABLES / "Table_S4_factorial_contrasts.csv": contrasts,
        TABLES / "Table_S5_pxr_aggregate_boundary.csv": pxr_s5,
        TABLES / "Table_S6_bootstrap_intervals.csv": boot,
        TABLES / "Table_S7_integrity_summary.csv": integrity,
    }
    for path, df in source_outputs.items():
        write_csv(df, path)

    style()
    figure_paths = []
    figure_paths.extend(save_figure(figure1(counts), "Figure_1_study_design"))
    figure_paths.extend(save_figure(figure2(boot), "Figure_2_confirmatory_signatures"))
    figure_paths.extend(save_figure(figure3(admission, direct, strata, overlap), "Figure_3_admission_translation_pxr"))
    figure_paths.extend(save_figure(figure_s1(factorial), "Figure_S1_nuisance_factorial"))
    figure_paths.extend(save_figure(graphical_abstract(), "Graphical_Abstract"))

    generated = list(source_outputs) + figure_paths
    report = {
        "status": "PASS",
        "script": str(Path(__file__).relative_to(BASE)).replace("\\", "/"),
        "script_sha256": sha256(Path(__file__)),
        "frozen_directory": str(FROZEN.relative_to(BASE)).replace("\\", "/"),
        "frozen_inventory": {"files": 19, "directories": 0},
        "accepted_hashes": fixed_hashes,
        "assertions": {
            "manifest_outputs_match": len(manifest["output_hashes_before_manifest"]),
            "semisynthetic_rows": len(scenarios),
            "admission_cells": len(admission),
            "direct_cells": len(direct_all),
            "factorial_summary_cells": len(factorial),
            "factorial_contrast_rows": len(contrasts),
            "bootstrap_rows": len(boot),
            "prediction_rows": len(predictions),
            "pxr_parent": int(strata["n_total"].sum()),
            "pxr_followup": int(strata["n_followup_available"].sum()),
        },
        "scientific_boundary": {
            "H1": "4/4 PASS",
            "H2": "4/4 PASS",
            "admission": "8/16; all MAR cells failed",
            "direct_translation": "NOT_ADJUDICATED_ADMISSION_FAILED",
            "PXR": "NOT_ESTIMATED_ASSUMPTIONS_AND_LOSS_UNESTABLISHED",
            "PXR_original_operational_state": "ABSTAIN_NONIDENTIFIABLE_FULL_LIBRARY_DOSE_RESPONSE",
        },
        "outputs": {str(p.relative_to(BASE)).replace("\\", "/"): sha256(p) for p in generated},
        "software": {"python": manifest["environment"]["python"], "pandas": pd.__version__, "numpy": np.__version__, "matplotlib": mpl.__version__},
    }
    report_path = QA / "paper_assets_build_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    md = [
        "# Paper assets build report",
        "",
        "Status: **PASS**",
        "",
        "The immutable accepted result directory was read only. Its exact 19-file inventory, four anchor hashes, and all pre-manifest output hashes matched before any paper asset was generated.",
        "",
        "## Scientific boundary retained",
        "",
        "- H1: 4/4 PASS.",
        "- H2: 4/4 PASS.",
        "- Admission: 8/16; every MAR cell failed.",
        "- Direct translation: NOT_ADJUDICATED_ADMISSION_FAILED.",
        "- PXR: not estimated; identification assumptions and the endpoint-loss contract remain unestablished. The original operational abstention label is retained only as historical provenance.",
        "",
        f"Generated {len(generated)} hashed figure, source-data and table files. Machine-readable hashes are in `{report_path.name}`.",
        "",
    ]
    (QA / "paper_assets_build_report.md").write_text("\n".join(md), encoding="utf-8")
    print(json.dumps({"status": "PASS", "generated": len(generated), "report": str(report_path)}, indent=2))


if __name__ == "__main__":
    main()
