from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "paper/source_data/Figure_5_multisplit_robustness.csv"
OUT = ROOT / "paper/figures/Figure_5_multisplit_robustness"


def main():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 9.5, "axes.labelsize": 10.5, "axes.edgecolor": "#6B7280",
        "axes.labelcolor": "#374151", "xtick.color": "#4B5563", "ytick.color": "#4B5563",
        "legend.fontsize": 8.5, "savefig.dpi": 300,
    })
    frame = pd.read_csv(DATA)
    x = np.arange(len(frame))
    colors = {"raw": "#B54A4A", "drop": "#2A7F8E", "low": "#6E5A9E"}
    markers = {"raw": "o", "drop": "s", "low": "^"}
    labels = {"raw": "Raw full loss", "drop": "Full loss without raw Ipc", "low": "Low-dimensional loss"}
    offsets = {"raw": -0.16, "drop": 0.0, "low": 0.16}

    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    for i, present in enumerate(frame.case_2945_in_target):
        if present:
            ax.axvspan(i - 0.43, i + 0.43, color="#F3D8CE", alpha=0.58, lw=0, zorder=0)
    for variant in ["raw", "drop", "low"]:
        ax.plot(
            x + offsets[variant], frame[variant], marker=markers[variant], markersize=5.8,
            color=colors[variant], lw=1.15, linestyle="-", label=labels[variant],
            markeredgecolor="white", markeredgewidth=0.55, zorder=3,
        )

    ax.set_yscale("log")
    ax.set_ylim(7e-4, 2.4e3)
    ax.set_xlim(-0.55, len(frame) - 0.45)
    ax.set_xticks(x)
    ax.set_xticklabels(frame.split_seed.astype(str), rotation=0)
    ax.set_xlabel("Prespecified structure-split seed")
    ax.set_ylabel("Median absolute AIPW risk error")
    ax.grid(axis="y", which="major", color="#D8DCE3", linewidth=0.6, alpha=0.75)
    ax.grid(axis="y", which="minor", visible=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3.5, width=0.8)
    patch = mpatches.Patch(facecolor="#F3D8CE", edgecolor="none", alpha=0.8, label="Mol2946 in target cohort")
    handles, legend_labels = ax.get_legend_handles_labels()
    ax.legend(
        handles + [patch], legend_labels + [patch.get_label()], loc="upper center", ncol=2,
        bbox_to_anchor=(0.5, 1.20),
        frameon=True, facecolor="white", edgecolor="#C8CDD5", framealpha=0.96,
        columnspacing=1.0, handletextpad=0.5, borderpad=0.45,
    )
    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.19, top=0.79)
    fig.savefig(OUT.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(OUT.with_suffix(".pdf"), facecolor="white")
    fig.savefig(OUT.with_suffix(".svg"), facecolor="white")
    fig.savefig(OUT.with_suffix(".tiff"), dpi=300, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(OUT.with_suffix(".png"))


if __name__ == "__main__":
    main()
