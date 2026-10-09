import argparse
from pathlib import Path

import pandas as pd

import build_paper_assets as assets
import build_revision_assets as revision


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output.resolve() if args.output else root / "replayed_figures"
    output.mkdir(parents=True, exist_ok=True)
    assets.FIGURES = output
    source = root / "paper/source_data"
    counts = pd.read_csv(source / "figure1_design_counts.csv")
    boot = pd.read_csv(source / "figure2_confirmation.csv")
    admission = pd.read_csv(source / "figure3_admission.csv")
    direct = pd.read_csv(source / "figure3_direct_masks.csv")
    strata = pd.read_csv(source / "figure3_pxr_strata.csv")
    overlap = pd.read_csv(source / "figure3_pxr_overlap.csv")
    factorial = pd.read_csv(source / "figureS1_factorial_summary.csv")
    calibration = pd.read_csv(source / "figure4_fixed_nuisance_calibration.csv", float_precision="round_trip")
    forensic = pd.read_csv(source / "figure4_ipc_sensitivity.csv", float_precision="round_trip")
    assert len(counts) == 12
    assert len(boot) == 12 and set(boot["resamples"]) == {10000}
    assert len(admission) == 16 and int(admission["pass"].sum()) == 8
    assert admission.groupby("mechanism")["pass"].sum().to_dict() == {"mar": 0, "mcar": 8}
    assert len(direct) == 8 and direct["overlap_eligible"].all()
    assert int(direct["aipw_better"].sum()) == 4
    assert len(strata) == 2 and int(strata["n_total"].sum()) == 10747
    assert int(strata["n_followup_available"].sum()) == 2722
    assert len(overlap) == 2 and not bool(overlap.loc[overlap["specification"] == "pxr_structure_screen_12", "overlap_pass"].iloc[0])
    assert len(factorial) == 16 and set(factorial["rows"]) == {400}
    assert len(calibration) == 32 and calibration["method"].eq("true_p_fixed_replay_m_aipw").all()
    assert len(forensic) == 3 and forensic["n_distinct"].eq(16).all()
    assets.style()
    generated = []
    generated += assets.save_figure(assets.figure1(counts), "Figure_1_study_design")
    generated += assets.save_figure(assets.figure2(boot), "Figure_2_confirmatory_signatures")
    generated += assets.save_figure(assets.figure3(admission, direct, strata, overlap), "Figure_3_admission_translation_pxr")
    generated += assets.save_figure(assets.figure_s1(factorial), "Figure_S1_nuisance_factorial")
    generated += assets.save_figure(assets.graphical_abstract(), "Graphical_Abstract")
    generated += assets.save_figure(revision.figure4(calibration, forensic), "Figure_4_diagnostic_mechanisms")
    assert len(generated) == 24 and all(p.stat().st_size > 0 for p in generated)
    print(f"PASS regenerated {len(generated)} files in {output}")


if __name__ == "__main__":
    main()
