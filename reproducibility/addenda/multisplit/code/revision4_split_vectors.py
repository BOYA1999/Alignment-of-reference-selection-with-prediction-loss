import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts/revision_004_multisplit/screen"
SOURCE = ROOT / "evidence/upstream/expansion_ml_comparison/biogen_adme_3521.csv"
SEEDS = range(61453, 61463)
EXPECTED_SOURCE_SHA256 = "cd310dd3911a9f9763d0581e6202bc20aaf8a810523d01e9ab5e403681a386b4"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def split_by_cluster(clusters, fraction, seed):
    sizes = pd.Series(clusters).value_counts()
    order = np.random.default_rng(seed).permutation(sizes.index.to_numpy())
    held_out, running = set(), 0
    for cluster in order:
        if running >= fraction * len(clusters):
            break
        held_out.add(cluster)
        running += int(sizes[cluster])
    return np.where(pd.Series(clusters).isin(held_out), "test", "train")


def main():
    assert importlib.metadata.version("numpy") == "2.5.2"
    assert importlib.metadata.version("pandas") == "3.0.5"
    assert sha(SOURCE) == EXPECTED_SOURCE_SHA256
    assert not OUT.exists()
    frame = pd.read_csv(SOURCE)
    rows, summaries = [], []
    for seed in SEEDS:
        labels = split_by_cluster(frame["cluster"].to_numpy(), 0.30, seed)
        if seed == 61453:
            assert np.array_equal(labels, frame["ds"].to_numpy())
        train_clusters = set(frame.loc[labels == "train", "cluster"])
        test_clusters = set(frame.loc[labels == "test", "cluster"])
        assert not train_clusters & test_clusters
        rows.extend(
            {"split_seed": seed, "source_row_index": i, "ds": label}
            for i, label in enumerate(labels)
        )
        summary = {
            "split_seed": seed,
            "train_rows": int(np.sum(labels == "train")),
            "test_rows": int(np.sum(labels == "test")),
            "held_out_clusters": len(test_clusters),
            "case_2945_in_test": bool(labels[2945] == "test"),
        }
        for endpoint in ["LOG_HLM", "LOG_RLM", "LOG_SOL", "LOG_MDR1_ER"]:
            finite = np.isfinite(pd.to_numeric(frame[endpoint], errors="coerce").to_numpy(float))
            summary[f"train_finite_{endpoint}"] = int(np.sum((labels == "train") & finite))
            summary[f"test_finite_{endpoint}"] = int(np.sum((labels == "test") & finite))
        summaries.append(summary)
    OUT.mkdir(parents=True)
    assignments = pd.DataFrame(rows)
    screen = pd.DataFrame(summaries)
    assignments.to_csv(OUT / "split_assignments.csv", index=False)
    screen.to_csv(OUT / "split_screen_counts.csv", index=False)
    manifest = {
        "status": "PASS_OUTCOME_BLIND_SPLIT_RECONSTRUCTION",
        "seeds": list(SEEDS),
        "rows_per_seed": len(frame),
        "assignment_rows": len(assignments),
        "original_seed_exact_match": True,
        "source_sha256": sha(SOURCE),
        "protocol_sha256": sha(ROOT / "artifacts/revision_004_multisplit/protocol.md"),
        "script_sha256": sha(Path(__file__)),
        "output_sha256": {
            "split_assignments.csv": sha(OUT / "split_assignments.csv"),
            "split_screen_counts.csv": sha(OUT / "split_screen_counts.csv"),
        },
    }
    (OUT / "split_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(screen.to_string(index=False))


if __name__ == "__main__":
    main()
