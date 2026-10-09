import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import time

sys.dont_write_bytecode = True

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts/revision_001/reconstruction/replay"
UPSTREAM = ROOT / "evidence/upstream/expansion_ml_comparison"
FROZEN = ROOT / "artifacts/confirmatory/confirm_001/full_recovery_001"
PINS = {"numpy": "2.5.2", "pandas": "3.0.5", "scipy": "1.17.1", "scikit-learn": "1.9.0", "rdkit": "2026.3.5", "bblean": "0.10.3"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def vector_sha(values):
    return hashlib.sha256("\n".join(str(value) for value in values).encode("utf-8")).hexdigest()


def main():
    started = time.perf_counter()
    assert sys.version_info[:2] == (3, 12), "Use the separate Python 3.12 split-replay environment."
    versions = {name: importlib.metadata.version(name) for name in PINS}
    assert versions == PINS, f"Historical split-replay versions are required: {versions}"
    assert not OUT.exists(), f"Output already exists; preserve it: {OUT}"
    cfg_path = ROOT / "artifacts/confirmatory/confirm_001/config.json"
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))["source_freeze"]["biogen"]
    provenance = cfg["split_provenance"]
    inputs = {UPSTREAM / name: digest for name, digest in provenance["construction_files_sha256"].items()}
    official = ROOT / cfg["official_data"]["local_path"]
    derivative = ROOT / cfg["analysis_derivative"]["path"]
    inputs[official] = cfg["official_data"]["sha256"]
    inputs[derivative] = cfg["analysis_derivative"]["sha256"]
    for path, expected in inputs.items():
        assert sha(path) == expected, f"Pinned input mismatch: {path}"
    preserved = {path: sha(path) for path in FROZEN.glob("*") if path.is_file()}
    if FROZEN.exists():
        assert len(preserved) == 19 and not any(path.is_dir() for path in FROZEN.iterdir())
    os.environ["ADME_DATASET"] = "biogen"
    sys.path.insert(0, str(UPSTREAM))
    spec = importlib.util.spec_from_file_location("pinned_biogen_preparation", UPSTREAM / "00b_prepare_biogen.py")
    preparation = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(preparation)
    rebuilt = preparation.clean(preparation.load_source(official, False))
    rebuilt["cluster"] = preparation.assign_clusters(rebuilt["SMILES"])
    rebuilt["ds"] = preparation.split_by_cluster(rebuilt["cluster"].to_numpy(), provenance["test_fraction"], provenance["seed"])
    accepted = pd.read_csv(derivative)
    assert len(rebuilt) == len(accepted) == 3521
    mismatches = {}
    for name in accepted:
        a, b = rebuilt[name].to_numpy(), accepted[name].to_numpy()
        equal = (a == b) | (pd.isna(a) & pd.isna(b))
        mismatches[name] = int((~equal).sum())
    amendment_path = ROOT / "data/external/source_freeze_20260911/source_manifest_amendment_001.json"
    amendment = json.loads(amendment_path.read_text(encoding="utf-8"))
    corrected = {item["vector"]: item["corrected_sha256"] for item in amendment["corrections"]}
    cluster_hash = vector_sha(rebuilt["cluster"])
    split_hash = vector_sha(rebuilt["ds"])
    hash_match = cluster_hash == corrected["ordered_cluster"] and split_hash == corrected["ordered_ds"]
    preserved_unchanged = all(sha(path) == digest for path, digest in preserved.items())
    inputs_unchanged = all(sha(path) == digest for path, digest in inputs.items())
    sizes = rebuilt["cluster"].value_counts()
    train_clusters = set(rebuilt.loc[rebuilt.ds.eq("train"), "cluster"])
    test_clusters = set(rebuilt.loc[rebuilt.ds.eq("test"), "cluster"])
    passed = not any(mismatches.values()) and hash_match and preserved_unchanged and inputs_unchanged
    OUT.mkdir(parents=True)
    local_csv = OUT / "biogen_adme_reconstructed_local_only.csv"
    rebuilt.to_csv(local_csv, index=False, lineterminator="\r\n")
    report = {
        "analysis_id": "revision_001_biogen_split_reconstruction",
        "status": "PASS_EXACT_SCIENTIFIC_COLUMNS_AND_VECTORS" if passed else "FAIL_RETAINED_MISMATCH",
        "scope": "Independent 3521-row upstream preparation and split reconstruction only; no fixed predictor or nuisance fit.",
        "python": platform.python_version(), "platform": platform.platform(), "package_versions": versions,
        "installed_distributions": {dist.metadata["Name"]: dist.version for dist in importlib.metadata.distributions()},
        "construction_revision": provenance["construction_revision"],
        "upstream_functions": ["clean", "assign_clusters", "split_by_cluster"],
        "clustering_contract": provenance["clustering"], "seed": provenance["seed"], "test_fraction": provenance["test_fraction"],
        "rows": len(rebuilt), "clusters": len(sizes), "singletons": int(sizes.eq(1).sum()), "largest_cluster": int(sizes.max()),
        "train_rows": int(rebuilt.ds.eq("train").sum()), "test_rows": int(rebuilt.ds.eq("test").sum()),
        "held_out_clusters": len(test_clusters), "clusters_straddling_split": len(train_clusters & test_clusters),
        "column_mismatch_rows": mismatches,
        "ordered_cluster_sha256_no_terminal_lf": cluster_hash, "ordered_ds_sha256_no_terminal_lf": split_hash,
        "accepted_vector_hashes_match": hash_match,
        "input_sha256": {path.relative_to(ROOT).as_posix(): digest for path, digest in inputs.items()},
        "config_sha256": sha(cfg_path), "amendment_sha256": sha(amendment_path), "script_sha256": sha(Path(__file__)),
        "source_inputs_unchanged": inputs_unchanged, "frozen_files_checked": len(preserved), "frozen_files_unchanged": preserved_unchanged,
        "reconstructed_csv_sha256": sha(local_csv), "accepted_csv_sha256": sha(derivative),
        "csv_bytes_equal": sha(local_csv) == sha(derivative),
        "redistribution": "Reconstructed CSV contains structures and assay values; keep local and exclude from the release package.",
        "elapsed_seconds": time.perf_counter() - started,
    }
    (OUT / "reconstruction_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in ("status", "rows", "clusters", "train_rows", "test_rows", "column_mismatch_rows", "accepted_vector_hashes_match", "csv_bytes_equal", "frozen_files_unchanged", "elapsed_seconds")}, indent=2))
    assert passed, "Reconstruction differs; report and local candidate were preserved."


if __name__ == "__main__":
    main()
