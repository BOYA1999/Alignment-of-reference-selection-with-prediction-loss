import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts/revision_004_multisplit/full"
TABLES = ROOT / "paper/tables"
DATA = ROOT / "paper/source_data"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    manifest = json.loads((SOURCE / "manifest.json").read_text(encoding="utf-8"))
    verification = json.loads((SOURCE / "independent_verification.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "PASS_FULL_VERIFIED_AFTER_SUMMARY_RECOVERY"
    assert verification["status"] == "PASS_INDEPENDENT_AGGREGATE_VERIFICATION"
    estimates = pd.read_csv(SOURCE / "estimates.csv", float_precision="round_trip")
    predictions = pd.read_csv(SOURCE / "predictions.csv", float_precision="round_trip")

    aipw = estimates[(estimates.method == "aipw") & (estimates.mechanism == "mar")]
    figure = aipw.groupby(
        ["split_seed", "case_2945_in_target", "loss_variant"], as_index=False
    ).abs_error.median().pivot(
        index=["split_seed", "case_2945_in_target"], columns="loss_variant", values="abs_error"
    ).reset_index()
    figure.columns.name = None
    figure = figure[["split_seed", "case_2945_in_target", "raw", "drop", "low"]]
    figure.to_csv(DATA / "Figure_5_multisplit_robustness.csv", index=False, float_format="%.17g")

    risks = predictions.groupby(["split_seed", "endpoint", "model"], as_index=False).agg(
        target_n=("source_row_index", "size"), mae=("loss", "mean")
    )
    risk = risks.pivot(index=["split_seed", "endpoint", "target_n"], columns="model", values="mae").reset_index()
    risk.columns.name = None
    errors = aipw.groupby(["split_seed", "endpoint", "loss_variant"], as_index=False).abs_error.median().pivot(
        index=["split_seed", "endpoint"], columns="loss_variant", values="abs_error"
    ).reset_index()
    errors.columns.name = None
    display = risk.merge(errors, on=["split_seed", "endpoint"], validate="one_to_one")
    case_map = estimates.groupby("split_seed").case_2945_in_target.first()
    display["case"] = display.split_seed.map(case_map).map({True: "Target", False: "Training"})
    display["raw_drop"] = display["raw"] / display["drop"]
    display = display[[
        "split_seed", "endpoint", "case", "target_n", "ecfp_rf", "ecfp_desc_xgb",
        "raw", "drop", "low", "raw_drop",
    ]]
    display.columns = [
        "Split", "Endpoint", "Mol2946", "Target N", "RF MAE", "XGB MAE",
        "Raw error", "Drop-Ipc error", "Low-dim error", "Raw/Drop",
    ]
    display.to_csv(TABLES / "Table_S17_multisplit_display.csv", index=False, float_format="%.6g")

    detailed = pd.read_csv(SOURCE / "aipw_loss_comparisons.csv", float_precision="round_trip")
    detailed.to_csv(TABLES / "Table_S17_multisplit_cell_comparisons.csv", index=False, float_format="%.17g")
    report = {
        "status": "PASS_REVISION4_ASSETS",
        "source_manifest_sha256": sha(SOURCE / "manifest.json"),
        "source_verification_sha256": sha(SOURCE / "independent_verification.json"),
        "outputs": {
            str(path.relative_to(ROOT)).replace("\\", "/"): sha(path)
            for path in [
                DATA / "Figure_5_multisplit_robustness.csv",
                TABLES / "Table_S17_multisplit_display.csv",
                TABLES / "Table_S17_multisplit_cell_comparisons.csv",
            ]
        },
    }
    (SOURCE / "paper_assets_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
