import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


expected = {}
for line in (ROOT / "manifest-sha256.txt").read_text(encoding="utf-8").splitlines():
    value, name = line.split("  ", 1)
    expected[name] = value
for name, value in expected.items():
    assert digest(ROOT / name) == value

estimates = pd.read_csv(ROOT / "results/estimates.csv", float_precision="round_trip")
parts = pd.concat(
    [pd.read_csv(ROOT / f"results/by_split/estimates_{seed}.csv", float_precision="round_trip") for seed in range(61453, 61463)],
    ignore_index=True,
)
keys = ["split_seed", "endpoint", "anchor", "mechanism", "mask_seed", "model", "method", "loss_variant"]
assert len(estimates) == len(parts) == 80000
assert not estimates.duplicated(keys).any()
assert estimates.groupby(["split_seed", "endpoint", "anchor", "mechanism", "mask_seed"]).size().eq(10).all()
left = estimates.sort_values(keys).reset_index(drop=True)
right = parts.sort_values(keys).reset_index(drop=True)
assert left[keys].equals(right[keys])
assert np.allclose(left[["truth", "estimate", "error", "abs_error", "top1_abs_share"]], right[["truth", "estimate", "error", "abs_error", "top1_abs_share"]], rtol=0, atol=0, equal_nan=True)

aipw = estimates[(estimates.method == "aipw") & (estimates.mechanism == "mar")]
cell = aipw.groupby(["split_seed", "case_2945_in_target", "endpoint", "anchor", "model", "loss_variant"]).abs_error.median().unstack("loss_variant")
cell["raw_drop"] = cell["raw"] / cell["drop"]
cell["raw_low"] = cell["raw"] / cell["low"]
present = cell[cell.index.get_level_values("case_2945_in_target")]
absent = cell[~cell.index.get_level_values("case_2945_in_target")]
assert len(present) == 32 and ((present.raw_drop >= 10) & (present.raw_low >= 10)).all()
assert len(absent) == 48 and not (absent.raw_drop >= 10).any()
assert aipw[(aipw.loss_variant == "raw") & aipw.case_2945_in_target].case_abs_rank.eq(1).all()

decision = json.loads((ROOT / "results/decision.json").read_text(encoding="utf-8"))
ordering = pd.read_csv(ROOT / "results/predictor_ordering.csv")
assert decision["conditional_mechanism_recurrence"]["cells_meeting_both_10x_rules"] == 32
assert decision["broader_raw_ipc_instability_without_case"]["cells_raw_to_drop_at_least_10x"] == 0
assert ordering.lower_mae_model.eq("ecfp_desc_xgb").sum() == 20
assert len(pd.read_csv(ROOT / "paper_assets/Table_S17_multisplit_display.csv")) == 20
assert len(pd.read_csv(ROOT / "paper_assets/Table_S17_multisplit_cell_comparisons.csv")) == 160
assert len(pd.read_csv(ROOT / "paper_assets/Figure_5_multisplit_robustness.csv")) == 10
print("PASS_CILS_MULTISPLIT_ADDENDUM")
