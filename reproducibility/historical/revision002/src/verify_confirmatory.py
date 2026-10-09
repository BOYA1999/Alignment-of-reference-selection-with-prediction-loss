import argparse
import hashlib
import json
import platform
import re
import sys
import warnings
from datetime import datetime
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import rdkit
import scipy
import sklearn
import xgboost
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.special import expit
from sklearn.ensemble import RandomForestRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor


EXPECTED_PROTOCOL_SHA256 = "1819fa04713ee1a37df66c847825e5de3adcba736fba47c5c40b1c8db7796d7a"
EXPECTED_CONFIG_SHA256 = "04fc27fc407a562853dd362d31a4c4cbc01bdb754199d5a18014cdfcbc7a5c23"
EXPECTED_RUNNER_SHA256 = "650a8b8f6b9bf65db69efc7ba60340aafa99c3c4b8b5bba957a8e6c6fa800c5c"
EXPECTED_AMENDMENT_HASHES = {
    "artifacts/confirmatory/confirm_001/protocol_amendment_001.md": "989eac8ebd6e26bdcb5cee454ec51f4c00e76ad93bf4b245b92e0a59e9ff4a49",
    "artifacts/confirmatory/confirm_001/config_amendment_001.json": "9008878d6a8c23aae059cf3b4aefcdc813bdd7387cbf1aa84fd6d38db620ec3a",
    "data/external/source_freeze_20260911/source_manifest_amendment_001.json": "5b6b3337a9b7b1a4d359062ab07807527accc37b5cb63798db69e5f23837dfe5",
}
EXPECTED_AMENDMENT_AUDIT = {
    "artifacts/confirmatory/confirm_001/amendment_001_independent_audit.md": "cc08015a2a8ea126b90cfd0a53af1a90d548a09f890ce86371490c0bdbc158cf"
}
EXPECTED_RECOVERY_HASHES = {
    "artifacts/confirmatory/confirm_001/protocol_amendment_002.md": "d34c1053da408a3e349c2ecdc20d1e49f8aa5cda9dd491746146d621c5988d9f",
    "artifacts/confirmatory/confirm_001/config_amendment_002.json": "144af9c376929f24699452051abed2ddb5243ab8cdf47c8f4024e219ac278984",
    "artifacts/confirmatory/confirm_001/failure_001.md": "a75a4b853bf9ce5469873870d0375e077f2dd2be7d5e6981febf5abd41667267",
    "artifacts/confirmatory/confirm_001/amendment_002_independent_audit.md": "7a266a7614729ccf5bb2ff13890eabb6d2e6e32cc0ebc66ea21a34c8f1a821c3",
}
DECISION_004_PATH = "artifacts/decisions/004_confirm_001_clear_to_smoke.md"
DECISION_005_PATH = "artifacts/decisions/005_confirm_001_clear_to_full.md"
STATIC_AUDIT_PATH = "artifacts/confirmatory/confirm_001/static_code_audit_001.md"
SMOKE_OUTPUT_PATH = "artifacts/confirmatory/confirm_001/smoke_001"
SMOKE_VERIFICATION_PATH = "artifacts/confirmatory/confirm_001/smoke_001_independent_verification.json"
DECISION_006_PATH = "artifacts/decisions/006_confirm_001_clear_to_recovery_smoke.md"
DECISION_007_PATH = "artifacts/decisions/007_confirm_001_clear_to_recovery_full.md"
RECOVERY_STATIC_AUDIT_PATH = "artifacts/confirmatory/confirm_001/recovery_001_static_code_audit.md"
RECOVERY_SMOKE_OUTPUT_PATH = "artifacts/confirmatory/confirm_001/recovery_smoke_001"
RECOVERY_SMOKE_VERIFICATION_PATH = "artifacts/confirmatory/confirm_001/recovery_smoke_001_independent_verification.json"
RECOVERY_FULL_OUTPUT_PATH = "artifacts/confirmatory/confirm_001/full_recovery_001"
RECOVERY_SMOKE_CHECK_KEYS = {
    "folds_complete",
    "fold_class_support",
    "raw_propensity_finite",
    "clipped_propensity_finite",
    "loss_nuisance_finite",
    "estimates_finite",
    "nuisance_dtype_float64",
    "csv_roundtrip_component_bit_exact",
    "csv_roundtrip_ratio_within_absolute_tolerance",
    "csv_roundtrip_large_ratio",
}
RECOVERY_SMOKE_SCOPE = {
    "independent_source_feature_prediction_generator_and_seed41001_formula_reconstruction": False,
    "independent_seed41001_crossfit_propensity_and_ridge_refit": False,
    "independent_remaining_99_seed_and_direct_or_pxr_nuisance_refit": False,
    "remaining_formal_nuisance_vectors_checked_by_hash_and_grid_invariants": False,
    "prelaunch_directory_absence_checked_from_runner_record_only": True,
}

ENDPOINTS = ("LOG_HLM", "LOG_RLM")
REGULAR_ANCHORS = ("LOG_SOL", "LOG_MDR1_ER")
STRESS_ANCHORS = ("LOG_HPPB", "LOG_RPPB")
MODELS = ("ecfp_rf", "ecfp_desc_xgb")
P_REPRESENTATIONS = ("lowdim_11", "full_2266")
M_REPRESENTATIONS = ("lowdim_11", "full_2266")
MECHANISMS = ("mcar", "mar")
SEEDS = tuple(range(41001, 41101))
LOWDIM_COLUMNS = (
    "MolWt",
    "MolLogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "RingCount",
    "FractionCSP3",
    "HeavyAtomCount",
    "MolMR",
    "max_train_similarity",
)
PXR_STRUCTURE_COLUMNS = LOWDIM_COLUMNS[:-1]
EXPECTED_ROWS = {
    "availability_map.csv": 8,
    "predictions.parquet": 3702,
    "nuisance_fold_assignments.parquet": 1851,
    "generator_vectors.parquet": 7404,
    "semisynthetic_scenarios.csv": 6400,
    "actual_mask_scenarios.csv": 64,
    "row_replay_seed41001.parquet": 59232,
    "factorial_contrasts.csv": 56,
    "bootstrap_intervals.csv": 12,
    "pxr_linkage_audit.csv": 10747,
    "pxr_selection_strata.csv": 2,
    "pxr_overlap_summary.csv": 2,
}
EXPECTED_OUTPUTS = {
    "cohort_audit.csv",
    "structure_audit.json",
    *EXPECTED_ROWS,
    "admission_gate.json",
    "pxr_boundary_decision.json",
    "summary.json",
    "verification.md",
}
FOLD_SUPPORT_COLUMNS = tuple(
    f"fold_{fold}_train_{state}"
    for fold in range(3)
    for state in ("observed", "unavailable")
)
SCHEMAS = {
    "cohort_audit.csv": (
        "dataset", "split", "endpoint", "n_rows", "n_finite", "n_unavailable", "coverage",
    ),
    "availability_map.csv": (
        "endpoint", "anchor", "n_total", "n_available", "n_unavailable", "coverage", "availability_sha256",
    ),
    "predictions.parquet": (
        "endpoint", "model", "source_row_index", "cohort_row_index", "molecule_id", "reference", "prediction", "loss",
    ),
    "nuisance_fold_assignments.parquet": (
        "endpoint", "source_row_index", "cohort_row_index", "molecule_id", "scaffold_group",
        "scaffold_sha256", "fold", "source_row_index_sha256", "scaffold_group_vector_sha256", "fold_sha256",
    ),
    "generator_vectors.parquet": (
        "endpoint", "anchor", "mechanism", "source_row_index", "cohort_row_index", "direct_coverage",
        "generator_score", "intercept_shift", "p_oracle", "p_oracle_sha256", "coverage_error", "all_finite",
    ),
    "semisynthetic_scenarios.csv": (
        "endpoint", "anchor", "mechanism", "seed", "model", "propensity_representation",
        "loss_representation", "n_total", "n_observed", "observed_fraction", "direct_anchor_coverage",
        "generator_intercept_shift", "source_row_index_sha256", "scaffold_group_vector_sha256", "fold_sha256",
        "R_sha256", "p_oracle_sha256", "p_raw_sha256", "p_fitted_sha256", "mhat_sha256", "full_truth", "cc",
        "or_mean", "fitted_ipw", "fitted_aipw", "fitted_se", "fitted_lcl", "fitted_ucl", "fitted_covered",
        "fitted_ess", "fitted_ess_fraction", "fitted_clip_fraction", "fitted_max_abs_smd",
        "fitted_max_normalized_weight", "fitted_smd_json", "oracle_ipw", "oracle_aipw", "oracle_se",
        "oracle_lcl", "oracle_ucl", "oracle_covered", "oracle_ess", "oracle_ess_fraction",
        "oracle_max_abs_smd", "oracle_max_normalized_weight", "oracle_smd_json", *FOLD_SUPPORT_COLUMNS,
        "fold_class_support", "all_finite", "cc_abs_error", "fitted_aipw_abs_error", "oracle_ipw_abs_error",
        "oracle_aipw_abs_error",
    ),
    "actual_mask_scenarios.csv": (
        "endpoint", "anchor", "mask_type", "model", "propensity_representation", "loss_representation",
        "n_total", "n_observed", "observed_fraction", "source_row_index_sha256", "scaffold_group_vector_sha256",
        "fold_sha256", "R_sha256", "p_raw_sha256", "p_fitted_sha256", "mhat_sha256", "full_truth", "cc",
        "or_mean", "fitted_ipw", "fitted_aipw", "fitted_se", "fitted_lcl", "fitted_ucl", "fitted_covered",
        "fitted_ess", "fitted_ess_fraction", "fitted_clip_fraction", "fitted_max_abs_smd",
        "fitted_max_normalized_weight", "fitted_smd_json", *FOLD_SUPPORT_COLUMNS, "fold_class_support", "all_finite",
        "overlap_eligible", "overlap_label", "cc_abs_error", "fitted_aipw_abs_error",
    ),
    "row_replay_seed41001.parquet": (
        "endpoint", "anchor", "mechanism", "seed", "model", "propensity_representation", "loss_representation",
        "source_row_index", "cohort_row_index", "molecule_id", "scaffold_group", "scaffold_sha256", "fold", "R",
        "loss", "mhat", "p_raw", "p_fitted", "p_oracle", *tuple(f"audit_{name}" for name in PXR_STRUCTURE_COLUMNS),
        "audit_max_train_similarity", "augmentation_fitted", "pseudo_outcome_fitted", "influence_fitted",
        "augmentation_oracle", "pseudo_outcome_oracle", "influence_oracle",
    ),
    "factorial_contrasts.csv": (
        "row_type", "endpoint", "anchor", "mechanism", "model", "loss_representation", "contrast", "estimate",
        "numerator", "denominator", "ratio", "rows", "status",
    ),
    "bootstrap_intervals.csv": (
        "hypothesis", "statistic", "endpoint", "mechanism", "point_estimate", "lower", "upper", "threshold",
        "direction", "point_pass", "denominator_valid", "resamples", "seed_index_matrix_sha256", "status",
    ),
    "pxr_linkage_audit.csv": (
        "source_row_index", "cohort_row_index", "OCNT_ID", "R_followup", "formal_hit", "smiles_parse_ok",
        "descriptor_10_all_finite", "canonical_smiles_sha256", "murcko_group_sha256", "fold", *PXR_STRUCTURE_COLUMNS,
        "log2_fc_estimate", "fdr_bh", "fdr_score", "p_raw_structure_10", "p_clipped_structure_10",
        "clipped_structure_10", "p_raw_structure_screen_12", "p_clipped_structure_screen_12",
        "clipped_structure_screen_12", "exact_smiles_match_if_linked", "batch_match_if_linked",
        "cross_file_exact_smiles_nunique", "cross_file_batch_nunique",
    ),
    "pxr_selection_strata.csv": (
        "formal_hit", "n_total", "n_followup_available", "n_followup_unavailable", "advancement_probability",
    ),
    "pxr_overlap_summary.csv": (
        "specification", "n_total", "n_observed", "n_unavailable", "coverage", "clip_fraction", "ess",
        "ess_fraction", "max_abs_smd", "max_normalized_weight", "smd_json", "zero_variance_features_json",
        *FOLD_SUPPORT_COLUMNS, "fold_class_support", "all_finite", "source_row_index_sha256",
        "murcko_group_vector_sha256", "fold_sha256", "p_raw_sha256", "p_clipped_sha256",
        "availability_sha256", "formal_hit_sha256", "trigger_clip", "trigger_ess", "trigger_smd", "overlap_pass",
    ),
}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def checked_fit(model, x, y, label):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning, module=r"sklearn\.linear_model\._logistic")
        warnings.simplefilter("error", ConvergenceWarning)
        fitted = model.fit(x, y)
    if isinstance(fitted, LogisticRegression):
        require(np.asarray(fitted.n_iter_).size > 0, f"missing iteration count: {label}")
        require(np.all(np.asarray(fitted.n_iter_) < fitted.max_iter), f"fit reached max_iter: {label}")
    return fitted


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def vector_hash(values, kind):
    if kind == "float":
        payload = np.asarray(values, dtype="<f8").tobytes(order="C")
    elif kind == "int":
        payload = np.asarray(values, dtype="<i8").tobytes(order="C")
    elif kind == "binary":
        payload = np.asarray(values, dtype="u1").tobytes(order="C")
    elif kind == "text":
        payload = "\n".join(str(value) for value in values).encode("utf-8")
    else:
        raise ValueError(f"unknown vector kind: {kind}")
    return hashlib.sha256(payload).hexdigest()


def row_text_hash(value):
    return hashlib.sha256(str(value).encode("utf-8")).hexdigest()


def close(left, right, tolerance, label):
    left = float(left)
    right = float(right)
    require(np.isfinite(left) and np.isfinite(right), f"nonfinite {label}")
    require(abs(left - right) <= tolerance, f"{label}: {left} != {right}")


def same_vector(left, right, tolerance, label):
    left = np.asarray(left)
    right = np.asarray(right)
    require(left.shape == right.shape, f"{label} shape mismatch")
    left_numeric = np.issubdtype(left.dtype, np.number) or np.issubdtype(left.dtype, np.bool_)
    right_numeric = np.issubdtype(right.dtype, np.number) or np.issubdtype(right.dtype, np.bool_)
    if left_numeric and right_numeric:
        require(
            np.allclose(left.astype(float), right.astype(float), rtol=0, atol=tolerance, equal_nan=True),
            f"{label} value mismatch",
        )
    else:
        require(np.array_equal(left.astype(str), right.astype(str)), f"{label} mismatch")


def compare_tree(expected, observed, tolerance, path="root"):
    if isinstance(expected, dict):
        require(isinstance(observed, dict), f"{path} is not an object")
        require(set(expected) == set(observed), f"{path} keys differ")
        for key in expected:
            compare_tree(expected[key], observed[key], tolerance, f"{path}.{key}")
        return
    if isinstance(expected, list):
        require(isinstance(observed, list), f"{path} is not a list")
        require(len(expected) == len(observed), f"{path} length differs")
        for index, (left, right) in enumerate(zip(expected, observed)):
            compare_tree(left, right, tolerance, f"{path}[{index}]")
        return
    if expected is None or isinstance(expected, (str, bool)):
        require(expected == observed, f"{path} differs")
        return
    if isinstance(expected, (int, float, np.integer, np.floating)):
        close(expected, observed, tolerance, path)
        return
    require(expected == observed, f"{path} differs")


def require_columns(frame, required, label):
    missing = set(required) - set(frame.columns)
    require(not missing, f"{label} missing columns: {sorted(missing)}")


def require_exact_grid(frame, keys, expected, label):
    require(not frame.duplicated(keys).any(), f"duplicate {label} key")
    observed = set(frame[list(keys)].itertuples(index=False, name=None))
    require(observed == set(expected), f"{label} Cartesian grid mismatch")


def finite_columns(frame, columns, label):
    require(
        np.isfinite(frame[list(columns)].to_numpy(dtype=float)).all(),
        f"nonfinite {label}",
    )


def canonical_and_group(smiles, source_row_index):
    molecule = Chem.MolFromSmiles(str(smiles))
    require(molecule is not None, f"SMILES parse failure at source row {source_row_index}")
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
    if scaffold.GetNumAtoms():
        group = "SCAFFOLD::" + Chem.MolToSmiles(
            scaffold, canonical=True, isomericSmiles=True
        )
    else:
        group = f"ACYCLIC_SOURCE_ROW_{int(source_row_index)}"
    return molecule, canonical, group


def fold_vector(groups):
    folds = np.full(len(groups), -1, dtype=np.int64)
    splitter = GroupKFold(n_splits=3, shuffle=False, random_state=None)
    for fold, (_, held_out) in enumerate(
        splitter.split(np.zeros(len(groups)), groups=np.asarray(groups, dtype=object))
    ):
        folds[held_out] = fold
    require(np.all(folds >= 0), "incomplete GroupKFold assignment")
    return folds


def fingerprints(molecules):
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048, includeChirality=True
    )
    return [generator.GetFingerprint(molecule) for molecule in molecules]


def fingerprint_matrix(bitvectors):
    matrix = np.empty((len(bitvectors), 2048), dtype=np.float32)
    for index, fingerprint in enumerate(bitvectors):
        row = np.zeros(2048, dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fingerprint, row)
        matrix[index] = row
    return matrix


def max_similarities(target_fingerprints, training_fingerprints):
    require(training_fingerprints, "empty endpoint prediction-training fingerprint set")
    return np.asarray(
        [max(DataStructs.BulkTanimotoSimilarity(fp, training_fingerprints)) for fp in target_fingerprints],
        dtype=float,
    )


def endpoint_count(cfg, endpoint):
    return int(cfg["source_freeze"]["biogen"]["exact_counts"]["test_finite"][endpoint])


def check_environment(cfg):
    expected = cfg["analysis_environment"]
    observed = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "rdkit": rdkit.__version__,
        "xgboost": xgboost.__version__,
    }
    for key, value in observed.items():
        require(value == expected[key], f"environment mismatch: {key}={value}, expected {expected[key]}")
    try:
        import pyarrow
    except ImportError as error:
        raise AssertionError("pyarrow unavailable") from error
    require(pyarrow.__version__ == expected["pyarrow"], "environment mismatch: pyarrow")
    require(platform.platform() == expected["platform"], "environment mismatch: platform")
    return observed | {"pyarrow": pyarrow.__version__, "platform": platform.platform()}


def amended_hashes(root):
    require(
        all(re.fullmatch(r"[0-9a-f]{64}", value) for value in EXPECTED_AMENDMENT_HASHES.values()),
        "independently cleared amendment hashes have not been installed in verifier",
    )
    for relative, wanted in EXPECTED_AMENDMENT_HASHES.items():
        path = root / relative
        require(path.is_file(), f"missing amendment: {relative}")
        require(sha256(path) == wanted, f"amendment hash mismatch: {relative}")
    for relative, wanted in EXPECTED_AMENDMENT_AUDIT.items():
        path = root / relative
        require(path.is_file(), f"missing independent amendment audit: {relative}")
        require(sha256(path) == wanted, f"amendment-audit hash mismatch: {relative}")
    config_overlay = json.loads(
        (root / "artifacts/confirmatory/confirm_001/config_amendment_001.json").read_text(
            encoding="utf-8"
        )
    )
    require(config_overlay["base_contract"]["protocol_sha256"] == EXPECTED_PROTOCOL_SHA256, "amendment base protocol mismatch")
    require(config_overlay["base_contract"]["config_sha256"] == EXPECTED_CONFIG_SHA256, "amendment base config mismatch")
    require(config_overlay["source_manifest_amendment"]["sha256"] == EXPECTED_AMENDMENT_HASHES["data/external/source_freeze_20260911/source_manifest_amendment_001.json"], "source amendment linkage mismatch")
    require(config_overlay["protocol_amendment"]["sha256"] == EXPECTED_AMENDMENT_HASHES["artifacts/confirmatory/confirm_001/protocol_amendment_001.md"], "protocol amendment linkage mismatch")
    require(
        config_overlay["formal_input_precedence"]
        == [
            "artifacts/confirmatory/confirm_001/protocol.md",
            "artifacts/confirmatory/confirm_001/config.json",
            "data/external/source_freeze_20260911/source_manifest.json",
            "data/external/source_freeze_20260911/source_manifest_amendment_001.json",
            "artifacts/confirmatory/confirm_001/protocol_amendment_001.md",
            "artifacts/confirmatory/confirm_001/config_amendment_001.json",
        ],
        "formal input precedence mismatch",
    )
    source_overlay = json.loads(
        (root / "data/external/source_freeze_20260911/source_manifest_amendment_001.json").read_text(
            encoding="utf-8"
        )
    )
    corrections = {row["vector"]: row["corrected_sha256"] for row in source_overlay["corrections"]}
    require(len(corrections) == 6, "source hash amendment must contain exactly six unique corrections")
    require(source_overlay["base_manifest"]["sha256"] == "be55555aed0b2a5bed095ae8d052500becd527c97c809d6db5b7a14e7548fd78", "source amendment base-manifest mismatch")
    return corrections, config_overlay


def load_recovery_amendment(root):
    for relative, wanted in EXPECTED_RECOVERY_HASHES.items():
        path = root / relative
        require(path.is_file() and sha256(path) == wanted, f"recovery evidence mismatch: {relative}")
    recovery = json.loads((root / "artifacts/confirmatory/confirm_001/config_amendment_002.json").read_text(encoding="utf-8"))
    require(recovery["study_id"] == "confirm_001" and recovery["amendment_id"] == "config_amendment_002", "recovery amendment identity mismatch")
    require(recovery["base_contract"]["protocol_sha256"] == EXPECTED_PROTOCOL_SHA256, "recovery base protocol mismatch")
    require(recovery["base_contract"]["config_sha256"] == EXPECTED_CONFIG_SHA256, "recovery base config mismatch")
    require(recovery["base_contract"]["config_amendment_001_sha256"] == EXPECTED_AMENDMENT_HASHES["artifacts/confirmatory/confirm_001/config_amendment_001.json"], "recovery predecessor mismatch")
    require(recovery["base_contract"]["protocol_amendment_002_sha256"] == EXPECTED_RECOVERY_HASHES["artifacts/confirmatory/confirm_001/protocol_amendment_002.md"], "recovery protocol mismatch")
    require(recovery["failed_execution"]["failure_record_sha256"] == EXPECTED_RECOVERY_HASHES["artifacts/confirmatory/confirm_001/failure_001.md"], "failure-record linkage mismatch")
    require(recovery["failed_execution"]["decision_004_sha256"] == "0e9ab1a55bbe4fa99409c17d1a2c8d9087b9f84f4eedc81ba819edc3b8af39d8", "historical Decision 004 mismatch")
    require(recovery["failed_execution"]["decision_005_sha256"] == "753068d66e5ca6dfd78dbcd1e2a56a81be87ea7908da214eb1a7ddd68675c025", "consumed Decision 005 mismatch")
    require(recovery["recovery"] == {
        "smoke_output_path": RECOVERY_SMOKE_OUTPUT_PATH,
        "smoke_verification_path": RECOVERY_SMOKE_VERIFICATION_PATH,
        "full_output_path": RECOVERY_FULL_OUTPUT_PATH,
        "runner_mode_smoke": "recovery-smoke",
        "runner_mode_full": "recovery-full",
        "paths_must_be_absent_before_launch": True,
        "write_once": True,
        "failed_output_is_evidence_only": True,
    }, "recovery path contract mismatch")
    contract = recovery["csv_roundtrip_contract"]
    require(
        contract["result_csv_write_float_format"] == "%.17g"
        and contract["result_csv_read_float_precision"] == "round_trip"
        and contract["apply_only_to_generated_result_csvs"] is True
        and contract["frozen_source_csv_readers_unchanged"] is True
        and contract["formula_replay_rtol"] == 0.0
        and contract["formula_replay_atol"] == 1e-12
        and contract["tolerance_relaxed"] is False
        and contract["recovery_smoke_requires_large_ratio_fixture"] is True
        and contract["recovery_smoke_minimum_ratio_magnitude"] == 100000.0,
        "recovery CSV contract mismatch",
    )
    failed = recovery["failed_execution"]
    failed_dir = root / failed["output_path"]
    require(failed_dir.is_dir(), "failed formal directory is missing")
    actual = {
        str(path.relative_to(failed_dir)).replace("\\", "/")
        for path in failed_dir.rglob("*")
        if path.is_file()
    }
    require(actual == set(failed["partial_inventory"]), "failed formal inventory changed")
    for name, record in failed["partial_inventory"].items():
        path = failed_dir / name
        require(path.stat().st_size == record["bytes"] and sha256(path) == record["sha256"], f"failed artifact changed: {name}")
    require(all(not (failed_dir / name).exists() for name in failed["required_absent_files"]), "failed run gained a completion artifact")
    for relative, wanted in failed["failed_code_snapshots"].items():
        path = root / relative
        require(path.is_file() and sha256(path) == wanted, f"failed code snapshot mismatch: {relative}")
    for key in ("decision_004", "decision_005", "failure_record"):
        path = root / failed[f"{key}_path"]
        require(path.is_file() and sha256(path) == failed[f"{key}_sha256"], f"failed-state evidence mismatch: {key}")
    audit = (root / "artifacts/confirmatory/confirm_001/amendment_002_independent_audit.md").read_text(encoding="utf-8")
    require("Verdict: CLEAR_TO_IMPLEMENT_RECOVERY" in audit, "recovery amendment was not cleared")
    return recovery


def validate_recovery_smoke_evidence(
    root, expected_inputs, expected_environment, runner_sha, verifier_sha, decision_006_sha, recovery,
):
    smoke_dir = root / RECOVERY_SMOKE_OUTPUT_PATH
    verification_path = root / RECOVERY_SMOKE_VERIFICATION_PATH
    require(smoke_dir.is_dir() and verification_path.is_file(), "recovery smoke evidence is missing")
    expected_files = {"roundtrip_fixture.csv", "smoke.json", "run_manifest.json"}
    entries = list(smoke_dir.iterdir())
    require({path.name for path in entries} == expected_files and all(path.is_file() for path in entries), "recovery smoke directory inventory mismatch")
    smoke_json_path = smoke_dir / "smoke.json"
    smoke_manifest_path = smoke_dir / "run_manifest.json"
    fixture_path = smoke_dir / "roundtrip_fixture.csv"
    hashes = {
        "Recovery smoke JSON SHA-256": sha256(smoke_json_path),
        "Recovery smoke run manifest SHA-256": sha256(smoke_manifest_path),
        "Recovery roundtrip fixture SHA-256": sha256(fixture_path),
        "Independent recovery smoke verification SHA-256": sha256(verification_path),
    }

    fixture = read_table(fixture_path)
    require(tuple(fixture.columns) == ("fixture_id", "numerator", "denominator", "ratio") and len(fixture) == 1, "recovery roundtrip fixture schema mismatch")
    numerator = 408.40724836686263
    denominator = 0.004053686965995372
    expected_fixture = pd.DataFrame({
        "fixture_id": ["synthetic_large_ratio"],
        "numerator": [numerator],
        "denominator": [denominator],
        "ratio": [numerator / denominator],
    })
    require(np.array_equal(fixture["fixture_id"].astype(str), expected_fixture["fixture_id"].astype(str)), "recovery roundtrip fixture identity mismatch")
    numeric_columns = ["numerator", "denominator", "ratio"]
    require(
        np.array_equal(
            fixture[numeric_columns].to_numpy(dtype="<f8").view("<u8"),
            expected_fixture[numeric_columns].to_numpy(dtype="<f8").view("<u8"),
        ),
        "recovery roundtrip fixture is not bit-exact",
    )
    stored_ratio = float(fixture.loc[0, "ratio"])
    ratio_difference = abs(stored_ratio - float(fixture.loc[0, "numerator"] / fixture.loc[0, "denominator"]))
    expected_contract = {
        "random_state": 20260911,
        "rows": 72,
        "groups": 24,
        "rows_per_group": 3,
        "availability_pattern": [0, 1],
        "propensity_representations": {"lowdim_11": [72, 11], "full_2266": [72, 2266]},
        "loss_columns": 2,
        "csv_roundtrip_fixture": {
            "synthetic_only": True,
            "float_format": "%.17g",
            "float_precision": "round_trip",
            "rows": 1,
            "ratio_magnitude": abs(stored_ratio),
            "derived_ratio_absolute_difference": ratio_difference,
        },
    }
    require(abs(stored_ratio) >= recovery["csv_roundtrip_contract"]["recovery_smoke_minimum_ratio_magnitude"] and ratio_difference <= 1e-12, "recovery roundtrip fixture threshold mismatch")

    smoke_json = json.loads(smoke_json_path.read_text(encoding="utf-8"))
    require(
        set(smoke_json) == {
            "study_id", "mode", "status", "non_adjudicative", "synthetic_contract", "checks",
            "fold_sha256", "input_hashes", "runner_sha256", "verifier_sha256",
        },
        "recovery smoke JSON top-level keys mismatch",
    )
    require(
        smoke_json["study_id"] == "confirm_001"
        and smoke_json["mode"] == "recovery-smoke"
        and smoke_json["status"] == "PASS"
        and smoke_json["non_adjudicative"] is True,
        "recovery smoke JSON identity mismatch",
    )
    require(smoke_json["runner_sha256"] == runner_sha and smoke_json["verifier_sha256"] == verifier_sha, "recovery smoke JSON code identity mismatch")
    require(smoke_json["input_hashes"] == expected_inputs, "recovery smoke JSON input hashes mismatch")
    require(smoke_json["synthetic_contract"] == expected_contract, "recovery smoke synthetic contract mismatch")
    require(set(smoke_json["checks"]) == RECOVERY_SMOKE_CHECK_KEYS and all(value is True for value in smoke_json["checks"].values()), "recovery smoke checks mismatch")
    groups = np.repeat([f"SMOKE_GROUP_{index}" for index in range(24)], 3)
    require(smoke_json["fold_sha256"] == vector_hash(fold_vector(groups), "int"), "recovery smoke fold hash mismatch")

    smoke_manifest = json.loads(smoke_manifest_path.read_text(encoding="utf-8"))
    require(
        set(smoke_manifest) == {
            "study_id", "mode", "status", "started_utc", "finished_utc", "duration_seconds", "command",
            "input_hashes", "runner_sha256", "verifier_sha256", "clearance_hashes", "environment",
            "synthetic_only", "output_hashes_before_manifest",
        },
        "recovery smoke manifest top-level keys mismatch",
    )
    require(smoke_manifest["study_id"] == "confirm_001" and smoke_manifest["mode"] == "recovery-smoke" and smoke_manifest["status"] == "PASS_NON_ADJUDICATIVE", "recovery smoke manifest identity mismatch")
    require(smoke_manifest["runner_sha256"] == runner_sha and smoke_manifest["verifier_sha256"] == verifier_sha, "recovery smoke manifest code identity mismatch")
    require(smoke_manifest["input_hashes"] == expected_inputs == smoke_json["input_hashes"], "recovery smoke manifest input hashes mismatch")
    require(smoke_manifest["environment"] == expected_environment, "recovery smoke manifest environment mismatch")
    require(smoke_manifest["clearance_hashes"] == {DECISION_006_PATH: decision_006_sha}, "recovery smoke clearance mismatch")
    require(smoke_manifest["synthetic_only"] is True, "recovery smoke manifest synthetic-only mismatch")
    require(
        smoke_manifest["output_hashes_before_manifest"]
        == {"roundtrip_fixture.csv": hashes["Recovery roundtrip fixture SHA-256"], "smoke.json": hashes["Recovery smoke JSON SHA-256"]},
        "recovery smoke manifest hash map mismatch",
    )
    command = smoke_manifest["command"]
    require(isinstance(command, list) and "--mode" in command and "recovery-smoke" in command and "--decision-006-sha" in command and decision_006_sha in command, "recovery smoke manifest command mismatch")
    require(float(smoke_manifest["duration_seconds"]) >= 0 and np.isfinite(float(smoke_manifest["duration_seconds"])), "recovery smoke manifest duration mismatch")
    started = datetime.fromisoformat(str(smoke_manifest["started_utc"]).replace("Z", "+00:00"))
    finished = datetime.fromisoformat(str(smoke_manifest["finished_utc"]).replace("Z", "+00:00"))
    require(finished >= started, "recovery smoke manifest timestamps reversed")

    verification = json.loads(verification_path.read_text(encoding="utf-8"))
    expected_verification = {
        "status": "PASS",
        "mode": "recovery-smoke",
        "non_adjudicative": True,
        "runner_sha256": runner_sha,
        "verifier_sha256": verifier_sha,
        "decision_006_sha256": decision_006_sha,
        "smoke_output_path": RECOVERY_SMOKE_OUTPUT_PATH,
        "smoke_json_sha256": hashes["Recovery smoke JSON SHA-256"],
        "smoke_run_manifest_sha256": hashes["Recovery smoke run manifest SHA-256"],
        "roundtrip_fixture_sha256": hashes["Recovery roundtrip fixture SHA-256"],
        "scope": RECOVERY_SMOKE_SCOPE,
    }
    require(set(verification) == set(expected_verification) and verification == expected_verification, "independent recovery smoke verification mismatch")
    return hashes


def recovery_authorization_hashes(
    root, mode, decision_006_sha, decision_007_sha, recovery, base_inputs, environment,
):
    def value(text, label):
        pattern = re.compile(rf"^\s*(?:[-*]\s*)?(?:\*\*)?{re.escape(label)}(?:\*\*)?\s*:\s*(.+?)\s*$")
        values = [match.group(1).strip().strip("`") for line in text.splitlines() if (match := pattern.match(line))]
        require(len(values) == 1 and values[0], f"missing or duplicate recovery clearance field: {label}")
        return values[0]

    require(mode in {"recovery-smoke", "recovery-full"}, "invalid recovery mode")
    require(re.fullmatch(r"[0-9a-f]{64}", decision_006_sha or ""), "trusted Decision 006 SHA-256 is required")
    decision_006 = root / DECISION_006_PATH
    require(decision_006.is_file() and sha256(decision_006) == decision_006_sha, "Decision 006 hash mismatch")
    text_006 = decision_006.read_text(encoding="utf-8")
    runner_sha = sha256(root / "src/run_confirmatory.py")
    verifier_sha = sha256(root / "src/verify_confirmatory.py")
    failed = recovery["failed_execution"]
    require(value(text_006, "Decision") == "CLEAR_TO_RECOVERY_SMOKE", "Decision 006 state mismatch")
    require(value(text_006, "Mode") == "recovery-smoke", "Decision 006 mode mismatch")
    require(value(text_006, "Runner SHA-256") == runner_sha, "Decision 006 runner mismatch")
    require(value(text_006, "Verifier SHA-256") == verifier_sha, "Decision 006 verifier mismatch")
    require(value(text_006, "Decision 004 SHA-256") == failed["decision_004_sha256"], "Decision 006 D004 mismatch")
    require(value(text_006, "Consumed Decision 005 SHA-256") == failed["decision_005_sha256"], "Decision 006 D005 mismatch")
    require(value(text_006, "Failed runner SHA-256") == failed["failed_code_snapshots"]["artifacts/confirmatory/confirm_001/failure_001_code_snapshot/run_confirmatory_failed.py"], "Decision 006 failed-runner mismatch")
    require(value(text_006, "Failed verifier SHA-256") == failed["failed_code_snapshots"]["artifacts/confirmatory/confirm_001/failure_001_code_snapshot/verify_confirmatory_failed.py"], "Decision 006 failed-verifier mismatch")
    fields = {
        "Failed output path": failed["output_path"],
        "Failure record path": failed["failure_record_path"],
        "Failure record SHA-256": failed["failure_record_sha256"],
        "Protocol Amendment 002 path": "artifacts/confirmatory/confirm_001/protocol_amendment_002.md",
        "Protocol Amendment 002 SHA-256": EXPECTED_RECOVERY_HASHES["artifacts/confirmatory/confirm_001/protocol_amendment_002.md"],
        "Config Amendment 002 path": "artifacts/confirmatory/confirm_001/config_amendment_002.json",
        "Config Amendment 002 SHA-256": EXPECTED_RECOVERY_HASHES["artifacts/confirmatory/confirm_001/config_amendment_002.json"],
        "Amendment 002 audit path": "artifacts/confirmatory/confirm_001/amendment_002_independent_audit.md",
        "Amendment 002 audit SHA-256": EXPECTED_RECOVERY_HASHES["artifacts/confirmatory/confirm_001/amendment_002_independent_audit.md"],
        "Amendment 002 audit verdict": "CLEAR_TO_IMPLEMENT_RECOVERY",
        "Recovery static audit path": RECOVERY_STATIC_AUDIT_PATH,
        "Recovery smoke output path": RECOVERY_SMOKE_OUTPUT_PATH,
        "Recovery full output path": RECOVERY_FULL_OUTPUT_PATH,
    }
    for label, expected in fields.items():
        require(value(text_006, label) == expected, f"Decision 006 field mismatch: {label}")
    static_path = root / RECOVERY_STATIC_AUDIT_PATH
    require(static_path.is_file(), "recovery static audit is missing")
    static_sha = sha256(static_path)
    require(value(text_006, "Recovery static audit SHA-256") == static_sha, "Decision 006 static-audit hash mismatch")
    require(value(text_006, "Recovery static audit verdict") == "CLEAR_TO_RECOVERY_SMOKE", "Decision 006 static-audit verdict mismatch")
    static_text = static_path.read_text(encoding="utf-8")
    require(value(static_text, "Verdict") == "CLEAR_TO_RECOVERY_SMOKE", "recovery static audit state mismatch")
    require(value(static_text, "Runner SHA-256") == runner_sha, "recovery static-audit runner mismatch")
    require(value(static_text, "Verifier SHA-256") == verifier_sha, "recovery static-audit verifier mismatch")
    records = {
        failed["decision_004_path"]: failed["decision_004_sha256"],
        failed["decision_005_path"]: failed["decision_005_sha256"],
        failed["failure_record_path"]: failed["failure_record_sha256"],
        **EXPECTED_RECOVERY_HASHES,
        RECOVERY_STATIC_AUDIT_PATH: static_sha,
        DECISION_006_PATH: decision_006_sha,
        "src/verify_confirmatory.py": verifier_sha,
        **failed["failed_code_snapshots"],
    }
    expected_smoke_inputs = {**base_inputs, **records, "src/run_confirmatory.py": runner_sha}
    if mode == "recovery-full":
        require(re.fullmatch(r"[0-9a-f]{64}", decision_007_sha or ""), "trusted Decision 007 SHA-256 is required")
        decision_007 = root / DECISION_007_PATH
        require(decision_007.is_file() and sha256(decision_007) == decision_007_sha, "Decision 007 hash mismatch")
        text_007 = decision_007.read_text(encoding="utf-8")
        require(value(text_007, "Decision") == "CLEAR_TO_RECOVERY_FULL", "Decision 007 state mismatch")
        require(value(text_007, "Mode") == "recovery-full", "Decision 007 mode mismatch")
        require(value(text_007, "Runner SHA-256") == runner_sha, "Decision 007 runner mismatch")
        require(value(text_007, "Verifier SHA-256") == verifier_sha, "Decision 007 verifier mismatch")
        require(value(text_007, "Decision 006 SHA-256") == decision_006_sha, "Decision 007 D006 mismatch")
        require(value(text_007, "Recovery smoke status") == "PASS", "Decision 007 smoke status mismatch")
        require(value(text_007, "Independent recovery smoke verifier status") == "PASS", "Decision 007 verifier status mismatch")
        require(value(text_007, "Recovery smoke output path") == RECOVERY_SMOKE_OUTPUT_PATH, "Decision 007 smoke path mismatch")
        require(value(text_007, "Recovery full output path") == RECOVERY_FULL_OUTPUT_PATH, "Decision 007 full path mismatch")
        require(value(text_007, "Failed output path") == failed["output_path"], "Decision 007 failed path mismatch")
        hashes = validate_recovery_smoke_evidence(
            root, expected_smoke_inputs, environment, runner_sha, verifier_sha, decision_006_sha, recovery,
        )
        for label, observed in hashes.items():
            require(value(text_007, label) == observed, f"Decision 007 evidence mismatch: {label}")
        require(value(text_007, "Independent recovery smoke verification path") == RECOVERY_SMOKE_VERIFICATION_PATH, "Decision 007 verification path mismatch")
        records.update({
            DECISION_007_PATH: decision_007_sha,
            f"{RECOVERY_SMOKE_OUTPUT_PATH}/smoke.json": hashes["Recovery smoke JSON SHA-256"],
            f"{RECOVERY_SMOKE_OUTPUT_PATH}/run_manifest.json": hashes["Recovery smoke run manifest SHA-256"],
            f"{RECOVERY_SMOKE_OUTPUT_PATH}/roundtrip_fixture.csv": hashes["Recovery roundtrip fixture SHA-256"],
            RECOVERY_SMOKE_VERIFICATION_PATH: hashes["Independent recovery smoke verification SHA-256"],
        })
    return records, runner_sha, verifier_sha


def authorization_hashes(root, mode, decision_004_sha, decision_005_sha):
    def clearance_value(text, label):
        pattern = re.compile(
            rf"^\s*(?:[-*]\s*)?(?:\*\*)?{re.escape(label)}(?:\*\*)?\s*:\s*(.+?)\s*$"
        )
        values = [
            match.group(1).strip().strip("`")
            for line in text.splitlines()
            if (match := pattern.match(line))
        ]
        require(len(values) == 1 and values[0], f"missing or duplicate clearance field: {label}")
        return values[0]

    require(re.fullmatch(r"[0-9a-f]{64}", decision_004_sha or ""), "a trusted Decision 004 SHA-256 must be supplied")
    decision_004 = root / DECISION_004_PATH
    require(decision_004.is_file(), "Decision 004 is missing")
    require(sha256(decision_004) == decision_004_sha, "Decision 004 hash mismatch")
    text_004 = decision_004.read_text(encoding="utf-8")
    runner_sha = sha256(root / "src/run_confirmatory.py")
    verifier_sha = sha256(root / "src/verify_confirmatory.py")
    require(clearance_value(text_004, "Decision") == "CLEAR_TO_SMOKE", "Decision 004 does not authorize smoke")
    require(clearance_value(text_004, "Mode") == "smoke", "Decision 004 mode mismatch")
    require(clearance_value(text_004, "Runner SHA-256") == runner_sha, "Decision 004 does not bind the current runner")
    require(clearance_value(text_004, "Verifier SHA-256") == verifier_sha, "Decision 004 does not bind the current verifier")
    require(clearance_value(text_004, "Smoke output path") == SMOKE_OUTPUT_PATH, "Decision 004 smoke path mismatch")
    require(clearance_value(text_004, "Static audit path") == STATIC_AUDIT_PATH, "Decision 004 static-audit path mismatch")
    static_audit = root / STATIC_AUDIT_PATH
    require(static_audit.is_file(), "static code audit is missing")
    static_audit_sha = sha256(static_audit)
    require(clearance_value(text_004, "Static audit SHA-256") == static_audit_sha, "Decision 004 static-audit hash mismatch")
    require(clearance_value(text_004, "Static audit verdict") == "CLEAR_TO_SMOKE", "Decision 004 static-audit verdict mismatch")
    static_text = static_audit.read_text(encoding="utf-8")
    require(clearance_value(static_text, "Verdict") == "CLEAR_TO_SMOKE", "static-audit state mismatch")
    require(clearance_value(static_text, "Runner SHA-256") == runner_sha, "static-audit runner hash mismatch")
    require(clearance_value(static_text, "Verifier SHA-256") == verifier_sha, "static-audit verifier hash mismatch")
    records = {
        DECISION_004_PATH: decision_004_sha,
        "src/verify_confirmatory.py": verifier_sha,
        STATIC_AUDIT_PATH: static_audit_sha,
    }
    if mode == "full":
        require(re.fullmatch(r"[0-9a-f]{64}", decision_005_sha or ""), "a trusted Decision 005 SHA-256 must be supplied")
        decision_005 = root / DECISION_005_PATH
        require(decision_005.is_file(), "Decision 005 is missing")
        require(sha256(decision_005) == decision_005_sha, "Decision 005 hash mismatch")
        text_005 = decision_005.read_text(encoding="utf-8")
        require(clearance_value(text_005, "Decision") == "CLEAR_TO_FULL", "Decision 005 does not authorize full run")
        require(clearance_value(text_005, "Mode") == "full", "Decision 005 mode mismatch")
        require(clearance_value(text_005, "Runner SHA-256") == runner_sha, "Decision 005 does not bind the current runner")
        require(clearance_value(text_005, "Verifier SHA-256") == verifier_sha, "Decision 005 does not bind the current verifier")
        require(clearance_value(text_005, "Decision 004 SHA-256") == decision_004_sha, "Decision 005 does not bind Decision 004")
        require(clearance_value(text_005, "Smoke status") == "PASS", "Decision 005 smoke status mismatch")
        require(clearance_value(text_005, "Independent verifier smoke status") == "PASS", "Decision 005 independent-verifier status mismatch")
        require(clearance_value(text_005, "Smoke output path") == SMOKE_OUTPUT_PATH, "Decision 005 smoke path mismatch")
        require(clearance_value(text_005, "Independent smoke verification path") == SMOKE_VERIFICATION_PATH, "Decision 005 smoke-verification path mismatch")
        smoke_output = (root / SMOKE_OUTPUT_PATH).resolve()
        smoke_json_path = smoke_output / "smoke.json"
        smoke_manifest_path = smoke_output / "run_manifest.json"
        smoke_verification_path = root / SMOKE_VERIFICATION_PATH
        require(smoke_output.is_dir() and smoke_json_path.is_file() and smoke_manifest_path.is_file(), "authorized smoke evidence is missing")
        require(smoke_verification_path.is_file(), "independent smoke verification is missing")
        smoke_json_sha = sha256(smoke_json_path)
        smoke_manifest_sha = sha256(smoke_manifest_path)
        smoke_verification_sha = sha256(smoke_verification_path)
        require(clearance_value(text_005, "Smoke JSON SHA-256") == smoke_json_sha, "Decision 005 smoke JSON hash mismatch")
        require(clearance_value(text_005, "Smoke run manifest SHA-256") == smoke_manifest_sha, "Decision 005 smoke manifest hash mismatch")
        require(clearance_value(text_005, "Independent smoke verification SHA-256") == smoke_verification_sha, "Decision 005 smoke-verification hash mismatch")
        smoke_json = json.loads(smoke_json_path.read_text(encoding="utf-8"))
        smoke_manifest = json.loads(smoke_manifest_path.read_text(encoding="utf-8"))
        smoke_verification = json.loads(smoke_verification_path.read_text(encoding="utf-8"))
        require(smoke_json["status"] == "PASS" and smoke_json["non_adjudicative"] is True, "smoke payload status mismatch")
        require(smoke_manifest["status"] == "PASS_NON_ADJUDICATIVE", "smoke manifest status mismatch")
        require(smoke_json["runner_sha256"] == runner_sha and smoke_manifest["runner_sha256"] == runner_sha, "smoke runner identity mismatch")
        require(smoke_json["verifier_sha256"] == verifier_sha and smoke_manifest["verifier_sha256"] == verifier_sha, "smoke verifier identity mismatch")
        require(smoke_manifest["clearance_hashes"] == {DECISION_004_PATH: decision_004_sha}, "smoke Decision 004 identity mismatch")
        require(smoke_manifest["output_hashes_before_manifest"].get("smoke.json") == smoke_json_sha, "smoke manifest payload hash mismatch")
        require(
            smoke_verification.get("status") == "PASS"
            and smoke_verification.get("mode") == "smoke"
            and smoke_verification.get("runner_sha256") == runner_sha
            and smoke_verification.get("verifier_sha256") == verifier_sha
            and smoke_verification.get("decision_004_sha256") == decision_004_sha
            and smoke_verification.get("smoke_json_sha256") == smoke_json_sha
            and smoke_verification.get("smoke_run_manifest_sha256") == smoke_manifest_sha,
            "independent smoke verification content mismatch",
        )
        records[DECISION_005_PATH] = decision_005_sha
        records[f"{SMOKE_OUTPUT_PATH}/smoke.json"] = smoke_json_sha
        records[f"{SMOKE_OUTPUT_PATH}/run_manifest.json"] = smoke_manifest_sha
        records[SMOKE_VERIFICATION_PATH] = smoke_verification_sha
    return records, runner_sha, verifier_sha


def expected_input_hashes(
    root, cfg, config_path, config_overlay, recovery, environment, mode,
    decision_004_sha, decision_005_sha, decision_006_sha, decision_007_sha,
):
    inputs = {
        "artifacts/confirmatory/confirm_001/config.json": EXPECTED_CONFIG_SHA256,
        "artifacts/confirmatory/confirm_001/protocol.md": EXPECTED_PROTOCOL_SHA256,
        "artifacts/confirmatory/confirm_001/protocol_audit.md": cfg["pre_result_protocol_audit"]["sha256"],
        cfg["decision_path"]: config_overlay["base_contract"]["decision_sha256"],
        cfg["collision_audit"]["path"]: cfg["collision_audit"]["sha256"],
        cfg["collision_audit"]["exact_claim_addendum"]["path"]: cfg["collision_audit"]["exact_claim_addendum"]["sha256"],
        cfg["source_freeze"]["audit_report"]["path"]: cfg["source_freeze"]["audit_report"]["sha256"],
        cfg["source_freeze"]["authority_manifest"]["path"]: cfg["source_freeze"]["authority_manifest"]["sha256"],
        cfg["source_freeze"]["biogen"]["official_data"]["local_path"]: cfg["source_freeze"]["biogen"]["official_data"]["sha256"],
        cfg["source_freeze"]["biogen"]["license"]["local_path"]: cfg["source_freeze"]["biogen"]["license"]["sha256"],
        cfg["source_freeze"]["biogen"]["official_readme"]["local_path"]: cfg["source_freeze"]["biogen"]["official_readme"]["sha256"],
        cfg["source_freeze"]["biogen"]["analysis_derivative"]["path"]: cfg["source_freeze"]["biogen"]["analysis_derivative"]["sha256"],
    }
    for record in cfg["source_freeze"]["pxr"]["files"].values():
        inputs[record["local_path"]] = record["sha256"]
    construction_root = "evidence/upstream/expansion_ml_comparison"
    for filename, wanted in cfg["source_freeze"]["biogen"]["split_provenance"]["construction_files_sha256"].items():
        inputs[f"{construction_root}/{filename}"] = wanted
    expansion = cfg["source_freeze"]["expansionrx_discovery"]
    inputs.update(
        {
            "data/raw/expansion_data_train.csv": expansion["train_sha256"],
            "data/raw/expansion_data_test.csv": expansion["test_sha256"],
            "data/raw/expansion_data_raw.csv": expansion["raw_sha256"],
        }
    )
    for relative, wanted in cfg["development_evidence"]["forensic_artifacts"].items():
        inputs[f"artifacts/diagnostics/forensic_001/{relative}"] = wanted
    inputs["artifacts/decisions/002_pilot_stop.md"] = cfg["development_evidence"]["pilot_decision_sha256"]
    inputs.update(EXPECTED_AMENDMENT_HASHES)
    inputs.update(EXPECTED_AMENDMENT_AUDIT)
    if mode in {"recovery-smoke", "recovery-full"}:
        inputs.update(EXPECTED_RECOVERY_HASHES)
        authorizations, runner_sha, verifier_sha = recovery_authorization_hashes(
            root, mode, decision_006_sha, decision_007_sha, recovery, dict(inputs), environment
        )
    else:
        authorizations, runner_sha, verifier_sha = authorization_hashes(
            root, mode, decision_004_sha, decision_005_sha
        )
    require(runner_sha == EXPECTED_RUNNER_SHA256, "runner differs from the frozen audited implementation")
    inputs.update(authorizations)
    inputs["src/run_confirmatory.py"] = runner_sha
    require(sha256(config_path) == EXPECTED_CONFIG_SHA256, "confirmatory config hash mismatch")
    for relative, wanted in inputs.items():
        path = root / relative
        require(path.is_file(), f"missing frozen input: {relative}")
        require(sha256(path) == wanted, f"frozen input hash mismatch: {relative}")
    return inputs, runner_sha, verifier_sha


def load_biogen_contexts(root, cfg, corrections):
    official_path = root / cfg["source_freeze"]["biogen"]["official_data"]["local_path"]
    derivative_path = root / cfg["source_freeze"]["biogen"]["analysis_derivative"]["path"]
    official = pd.read_csv(official_path)
    data = pd.read_csv(derivative_path)
    require(len(official) == len(data) == 3521, "Biogen row-count mismatch")
    expected_official_columns = cfg["source_freeze"]["biogen"]["official_data"]["column_names"]
    expected_derived_columns = cfg["source_freeze"]["biogen"]["analysis_derivative"]["column_names"]
    require(list(official.columns) == expected_official_columns, "Biogen official schema mismatch")
    require(list(data.columns) == expected_derived_columns, "Biogen derivative schema mismatch")
    require(
        vector_hash(
            official["Internal ID"].astype(str) + "\t" + official["SMILES"].astype(str),
            "text",
        )
        == corrections["ordered_internal_id_tab_smiles"],
        "amended Biogen ordered ID/SMILES hash mismatch",
    )
    mapping = cfg["source_freeze"]["biogen"]["source_to_derived_mapping"]
    for source, derived in mapping.items():
        if source in ("Internal ID", "SMILES"):
            require(
                np.array_equal(official[source].astype(str), data[derived].astype(str)),
                f"Biogen source-to-derivative mismatch: {source}",
            )
        else:
            left = pd.to_numeric(official[source], errors="coerce").to_numpy(float)
            right = pd.to_numeric(data[derived], errors="coerce").to_numpy(float)
            require(np.array_equal(left, right, equal_nan=True), f"Biogen mapping mismatch: {source}")
    require(vector_hash(data["cluster"].astype(str), "text") == corrections["ordered_cluster"], "amended cluster hash mismatch")
    require(vector_hash(data["ds"].astype(str), "text") == corrections["ordered_ds"], "amended split hash mismatch")
    require((data["ds"] == "train").sum() == 2463 and (data["ds"] == "test").sum() == 1058, "Biogen split counts mismatch")
    require(data.groupby("cluster")["ds"].nunique().max() == 1, "Biogen frozen cluster crosses split")

    molecules = []
    canonicals = []
    groups = []
    for source_row_index, smiles in enumerate(data["SMILES"]):
        molecule, canonical, group = canonical_and_group(smiles, source_row_index)
        molecules.append(molecule)
        canonicals.append(canonical)
        groups.append(group)
    require(len(set(canonicals)) == 3521, "Biogen canonical-SMILES uniqueness mismatch")
    train_canonical = {canonicals[index] for index in np.flatnonzero(data["ds"].eq("train"))}
    test_canonical = {canonicals[index] for index in np.flatnonzero(data["ds"].eq("test"))}
    require(not train_canonical.intersection(test_canonical), "canonical structure crosses Biogen split")
    all_fingerprints = fingerprints(molecules)
    all_fingerprint_matrix = fingerprint_matrix(all_fingerprints)
    descriptor_names = [name for name, _ in sorted(Descriptors._descList, key=lambda item: item[0])]
    require(len(descriptor_names) == 217, "RDKit descriptor count mismatch")
    require(
        vector_hash(descriptor_names, "text")
        == "84cea0bf91c3c74a704770fd5d5477dc689fc3c09e0997cd268d070acc61be73",
        "RDKit descriptor-name vector hash mismatch",
    )
    descriptor_functions = dict(Descriptors._descList)
    descriptor_matrix = np.asarray(
        [[descriptor_functions[name](molecule) for name in descriptor_names] for molecule in molecules],
        dtype=float,
    )
    require(np.isfinite(descriptor_matrix).all(), "nonfinite Biogen 217-descriptor block")
    descriptor_index = {name: index for index, name in enumerate(descriptor_names)}
    base_features = descriptor_matrix[:, [descriptor_index[name] for name in PXR_STRUCTURE_COLUMNS]]

    contexts = {}
    for endpoint in ENDPOINTS:
        y = pd.to_numeric(data[endpoint], errors="coerce").to_numpy(float)
        target_source = np.flatnonzero(data["ds"].eq("test").to_numpy() & np.isfinite(y))
        training_source = np.flatnonzero(data["ds"].eq("train").to_numpy() & np.isfinite(y))
        require(len(target_source) == endpoint_count(cfg, endpoint), f"{endpoint} target count mismatch")
        similarity = max_similarities(
            [all_fingerprints[index] for index in target_source],
            [all_fingerprints[index] for index in training_source],
        )
        lowdim = np.column_stack((base_features[target_source], similarity)).astype(np.float64)
        full = np.hstack(
            (
                all_fingerprint_matrix[target_source],
                descriptor_matrix[target_source],
                similarity[:, None],
            )
        ).astype(np.float64)
        endpoint_groups = np.asarray([groups[index] for index in target_source], dtype=object)
        folds = fold_vector(endpoint_groups)
        contexts[endpoint] = {
            "source_row_index": target_source.astype(np.int64),
            "cohort_row_index": np.arange(len(target_source), dtype=np.int64),
            "name": data.iloc[target_source]["Name"].astype(str).to_numpy(),
            "smiles": data.iloc[target_source]["SMILES"].astype(str).to_numpy(),
            "canonical": np.asarray([canonicals[index] for index in target_source], dtype=object),
            "groups": endpoint_groups,
            "fold": folds,
            "y": y[target_source],
            "lowdim": lowdim,
            "representations": {"lowdim_11": lowdim, "full_2266": full},
            "availability": {
                anchor: np.isfinite(pd.to_numeric(data[anchor], errors="coerce").to_numpy(float)[target_source]).astype(np.uint8)
                for anchor in REGULAR_ANCHORS + STRESS_ANCHORS
            },
            "source_hash": vector_hash(target_source, "int"),
            "canonical_hash": vector_hash([canonicals[index] for index in target_source], "text"),
            "group_hash": vector_hash(endpoint_groups, "text"),
            "fold_hash": vector_hash(folds, "int"),
            "lowdim_hash": vector_hash(lowdim.ravel(order="C"), "float"),
        }
    return data, contexts, {
        "ordered_canonical_smiles_sha256": vector_hash(canonicals, "text"),
        "unique_canonical_smiles": len(set(canonicals)),
        "cross_split_canonical_overlap": 0,
        "descriptor_names": descriptor_names,
        "descriptor_count": len(descriptor_names),
        "fingerprint_matrix": all_fingerprint_matrix,
        "descriptor_matrix": descriptor_matrix,
    }


def shifted_mar_probability(features, availability, cfg):
    scaled = StandardScaler().fit_transform(np.asarray(features, dtype=float))
    model_cfg = cfg["mask_contract"]["mar"]["generator_model"]
    model = checked_fit(LogisticRegression(
        penalty=model_cfg["penalty"],
        C=float(model_cfg["C"]),
        solver=model_cfg["solver"],
        class_weight=model_cfg["class_weight"],
        max_iter=int(model_cfg["max_iter"]),
        random_state=int(model_cfg["random_state"]),
    ), scaled, availability, "MAR generator")
    score = model.decision_function(scaled)
    solver_cfg = cfg["mask_contract"]["mar"]["intercept_shift_solver"]
    lo, hi = float(solver_cfg["initial_lower"]), float(solver_cfg["initial_upper"])
    clip_lo, clip_hi = map(float, cfg["mask_contract"]["mar"]["oracle_clip"])
    target = float(np.mean(availability))
    for _ in range(int(cfg["mask_contract"]["mar"]["intercept_shift_solver"]["updates"])):
        midpoint = (lo + hi) / 2
        if np.clip(expit(score + midpoint), clip_lo, clip_hi).mean() < target:
            lo = midpoint
        else:
            hi = midpoint
    delta = (lo + hi) / 2
    probability = np.clip(expit(score + delta), clip_lo, clip_hi)
    require(
        abs(float(probability.mean()) - target)
        <= float(cfg["mask_contract"]["mar"]["intercept_shift_absolute_tolerance"]),
        "MAR generator coverage-match failure",
    )
    return probability, float(delta), np.asarray(score, dtype=float)


def generator_contracts(contexts, cfg):
    contracts = {}
    for endpoint, context in contexts.items():
        for anchor in REGULAR_ANCHORS:
            available = context["availability"][anchor]
            coverage = float(available.mean())
            mcar = np.full(len(available), coverage, dtype=float)
            mar, delta, mar_score = shifted_mar_probability(context["lowdim"], available, cfg)
            for mechanism, probability, shift in (
                ("mcar", mcar, None),
                ("mar", mar, delta),
            ):
                contracts[(endpoint, anchor, mechanism)] = {
                    "probability": probability,
                    "intercept_shift": shift,
                    "generator_score": np.zeros(len(available), dtype=float) if mechanism == "mcar" else mar_score,
                    "coverage": coverage,
                    "probability_hash": vector_hash(probability, "float"),
                }
    return contracts


def regenerated_mask(context, contract, seed):
    probability = contract["probability"]
    uniform = np.random.default_rng(int(seed)).random(len(probability))
    mask = (uniform < probability).astype(np.uint8)
    return mask, vector_hash(mask, "binary")


def load_pxr_context(root, cfg, corrections):
    files = cfg["source_freeze"]["pxr"]["files"]
    screen = pd.read_csv(root / files["single_concentration_train"]["local_path"])
    dose = pd.read_csv(root / files["dose_response_train"]["local_path"])
    require(len(screen) == 21003 and len(dose) == 4139, "PXR source row-count mismatch")
    require(
        list(screen.columns) == files["single_concentration_train"]["column_names"],
        "PXR screen schema mismatch",
    )
    require(
        list(dose.columns) == files["dose_response_train"]["column_names"],
        "PXR dose-response schema mismatch",
    )
    screen["source_row_index"] = np.arange(len(screen), dtype=np.int64)
    parent = screen[
        screen["compound_class"].eq("Library")
        & screen["concentration_M"].eq(8.251e-06)
    ].copy()
    parent = parent.sort_values("source_row_index").reset_index(drop=True)
    parent["cohort_row_index"] = np.arange(len(parent), dtype=np.int64)
    require(len(parent) == 10747 and parent["OCNT_ID"].nunique() == 10747, "PXR parent mismatch")
    dose_ids = set(dose["OCNT_ID"].astype(str))
    available = parent["OCNT_ID"].astype(str).isin(dose_ids).astype(np.uint8).to_numpy()
    log2fc = pd.to_numeric(parent["log2_fc_estimate"], errors="coerce").to_numpy(float)
    fdr = pd.to_numeric(parent["fdr_bh"], errors="coerce").to_numpy(float)
    require(np.isfinite(log2fc).all() and np.isfinite(fdr).all(), "nonfinite PXR screen fields")
    formal_hit = ((log2fc > 1.0) & (fdr < 0.05)).astype(np.uint8)
    require(int(available.sum()) == 2722, "PXR availability count mismatch")
    require(int(formal_hit.sum()) == 873, "PXR formal-hit count mismatch")
    require(
        vector_hash(parent["OCNT_ID"].astype(str), "text")
        == corrections["ordered_parent_OCNT_ID"],
        "amended PXR parent ID hash mismatch",
    )
    require(
        vector_hash(available, "binary") == corrections["ordered_availability"],
        "amended PXR availability hash mismatch",
    )
    require(
        vector_hash(formal_hit, "binary") == corrections["ordered_formal_hit"],
        "amended PXR formal-hit hash mismatch",
    )

    molecules = []
    canonicals = []
    groups = []
    descriptor_rows = []
    descriptor_registry = dict(Descriptors._descList)
    for source_row_index, smiles in zip(parent["source_row_index"], parent["SMILES"]):
        molecule, canonical, group = canonical_and_group(smiles, source_row_index)
        molecules.append(molecule)
        canonicals.append(canonical)
        groups.append(group)
        descriptor_rows.append([descriptor_registry[name](molecule) for name in PXR_STRUCTURE_COLUMNS])
    descriptors = np.asarray(descriptor_rows, dtype=float)
    require(descriptors.shape == (10747, 10), "PXR descriptor shape mismatch")
    require(np.isfinite(descriptors).all(), "nonfinite PXR structure descriptor")
    require(len(set(canonicals)) == 10747, "PXR canonical-SMILES uniqueness mismatch")
    groups = np.asarray(groups, dtype=object)
    folds = fold_vector(groups)
    fdr_score = -np.log10(np.maximum(fdr, 1e-12))
    require(np.isfinite(fdr_score).all() and np.all((fdr_score >= 0) & (fdr_score <= 12)), "invalid PXR fdr_score")

    dose_by_id = dose.assign(__id=dose["OCNT_ID"].astype(str)).set_index("__id", drop=True)
    combined = pd.concat(
        [screen[["OCNT_ID", "SMILES", "OCNT Batch"]], dose[["OCNT_ID", "SMILES", "OCNT Batch"]]],
        ignore_index=True,
    )
    cardinality = combined.groupby("OCNT_ID").agg(
        cross_file_exact_smiles_nunique=("SMILES", "nunique"),
        cross_file_batch_nunique=("OCNT Batch", "nunique"),
    )
    exact_match = np.full(len(parent), None, dtype=object)
    batch_match = np.full(len(parent), None, dtype=object)
    smiles_nunique = parent["OCNT_ID"].map(cardinality["cross_file_exact_smiles_nunique"]).to_numpy(np.int64)
    batch_nunique = parent["OCNT_ID"].map(cardinality["cross_file_batch_nunique"]).to_numpy(np.int64)
    for index, row in parent.iterrows():
        key = str(row["OCNT_ID"])
        if key in dose_by_id.index:
            linked = dose_by_id.loc[key]
            exact_match[index] = bool(str(row["SMILES"]) == str(linked["SMILES"]))
            batch_match[index] = bool(str(row["OCNT Batch"]) == str(linked["OCNT Batch"]))
    require(sum(value is True for value in exact_match) == 2375, "PXR exact-SMILES link count mismatch")
    require(sum(value is True for value in batch_match) == 2700, "PXR batch link count mismatch")
    require(int((cardinality["cross_file_exact_smiles_nunique"] > 1).sum()) == 351, "PXR cross-file SMILES disagreement mismatch")
    require(int((cardinality["cross_file_batch_nunique"] > 1).sum()) == 23, "PXR cross-file batch disagreement mismatch")

    expected = pd.DataFrame(
        {
            "source_row_index": parent["source_row_index"].to_numpy(np.int64),
            "cohort_row_index": parent["cohort_row_index"].to_numpy(np.int64),
            "OCNT_ID": parent["OCNT_ID"].astype(str).to_numpy(),
            "R_followup": available,
            "formal_hit": formal_hit,
            "smiles_parse_ok": np.ones(len(parent), dtype=np.uint8),
            "descriptor_10_all_finite": np.ones(len(parent), dtype=np.uint8),
            "canonical_smiles_sha256": [row_text_hash(value) for value in canonicals],
            "murcko_group_sha256": [row_text_hash(value) for value in groups],
            "fold": folds,
            **{name: descriptors[:, index] for index, name in enumerate(PXR_STRUCTURE_COLUMNS)},
            "log2_fc_estimate": log2fc,
            "fdr_bh": fdr,
            "fdr_score": fdr_score,
            "exact_smiles_match_if_linked": exact_match,
            "batch_match_if_linked": batch_match,
            "cross_file_exact_smiles_nunique": smiles_nunique,
            "cross_file_batch_nunique": batch_nunique,
        }
    )
    return {
        "parent": parent,
        "expected_linkage": expected,
        "groups": groups,
        "fold": folds,
        "available": available,
        "formal_hit": formal_hit,
        "features": {
            "pxr_structure_10": descriptors,
            "pxr_structure_screen_12": np.column_stack((descriptors, log2fc, fdr_score)),
        },
        "hashes": {
            "ordered_source_row_index": vector_hash(parent["source_row_index"], "int"),
            "ordered_OCNT_ID": vector_hash(parent["OCNT_ID"].astype(str), "text"),
            "ordered_canonical_smiles": vector_hash(canonicals, "text"),
            "ordered_murcko_group": vector_hash(groups, "text"),
            "ordered_fold": vector_hash(folds, "int"),
            "ordered_availability": vector_hash(available, "binary"),
            "ordered_formal_hit": vector_hash(formal_hit, "binary"),
        },
    }


def cluster_se(influence, groups):
    codes = pd.factorize(np.asarray(groups, dtype=object), sort=False)[0]
    cluster_sums = np.bincount(codes, weights=np.asarray(influence, dtype=float))
    count = len(cluster_sums)
    require(count > 1, "fewer than two scaffold groups")
    return float(np.sqrt((count / (count - 1)) * np.square(cluster_sums).sum()) / len(influence))


def smd_diagnostics(mask, probability, features):
    mask = np.asarray(mask, dtype=float)
    probability = np.asarray(probability, dtype=float)
    features = np.asarray(features, dtype=float)
    weights = mask / probability
    require(weights.sum() > 0, "zero fitted observed weight")
    standard_deviation = features.std(axis=0, ddof=0)
    keep = standard_deviation > 0
    require(np.any(keep), "all audit features have zero variance")
    weighted = np.sum(weights[:, None] * features, axis=0) / weights.sum()
    smd = np.abs(weighted[keep] - features[:, keep].mean(axis=0)) / standard_deviation[keep]
    return float(smd.max()), int((~keep).sum())


def smd_map(mask, probability, features, names=LOWDIM_COLUMNS):
    mask = np.asarray(mask, dtype=float)
    probability = np.asarray(probability, dtype=float)
    features = np.asarray(features, dtype=float)
    weights = mask / probability
    standard_deviation = features.std(axis=0, ddof=0)
    weighted = np.sum(weights[:, None] * features, axis=0) / weights.sum()
    values = {}
    for index, name in enumerate(names):
        if standard_deviation[index] > 0:
            values[name] = float(abs(weighted[index] - features[:, index].mean()) / standard_deviation[index])
    return values


def estimate_from_rows(loss, mask, probability, mhat, groups):
    loss = np.asarray(loss, dtype=float)
    mask = np.asarray(mask, dtype=float)
    probability = np.asarray(probability, dtype=float)
    mhat = np.asarray(mhat, dtype=float)
    require(
        np.isfinite(np.concatenate((loss, mask, probability, mhat))).all(),
        "nonfinite estimator input",
    )
    require(np.isin(mask, [0.0, 1.0]).all(), "nonbinary estimator mask")
    require(np.all((probability > 0) & (probability <= 1)), "invalid propensity")
    observed = mask.astype(bool)
    require(observed.any(), "empty observed set")
    weights = mask / probability
    require(weights.sum() > 0 and np.square(weights).sum() > 0, "zero estimator weights")
    augmentation = mask * (loss - mhat) / probability
    pseudo = mhat + augmentation
    influence = pseudo - pseudo.mean()
    aipw = float(pseudo.mean())
    se = cluster_se(influence, groups)
    ess = float(weights.sum() ** 2 / np.square(weights).sum())
    return {
        "truth": float(loss.mean()),
        "cc": float(loss[observed].mean()),
        "or_mean": float(mhat.mean()),
        "ipw": float(np.sum(weights * loss) / weights.sum()),
        "aipw": aipw,
        "aipw_se": se,
        "aipw_lcl": aipw - 1.96 * se,
        "aipw_ucl": aipw + 1.96 * se,
        "coverage": bool(aipw - 1.96 * se <= loss.mean() <= aipw + 1.96 * se),
        "ess": ess,
        "ess_over_observed": ess / float(mask.sum()),
        "maximum_normalized_observed_weight": float(weights.max() / weights.sum()),
        "augmentation": augmentation,
        "pseudo_outcome": pseudo,
        "influence": influence,
    }


def safe_ratio(numerator, denominator):
    numerator = float(numerator)
    denominator = float(denominator)
    if not np.isfinite(numerator) or not np.isfinite(denominator) or denominator == 0:
        return None
    value = numerator / denominator
    return float(value) if np.isfinite(value) else None


def quantile(values, probability):
    values = np.asarray(values, dtype=float)
    require(np.isfinite(values).all(), "nonfinite quantile input")
    return float(np.quantile(values, probability, method="linear"))


def hypothesis_details(scenarios, cfg):
    selected_full = scenarios[
        scenarios["propensity_representation"].eq("full_2266")
        & scenarios["loss_representation"].eq("lowdim_11")
    ]
    selected_low = scenarios[
        scenarios["propensity_representation"].eq("lowdim_11")
        & scenarios["loss_representation"].eq("lowdim_11")
    ]
    h1_threshold = cfg["confirmatory_hypotheses"]["H1_PROPENSITY_INSTABILITY_CONFIRMED"]["maximum_ratio"]
    h2_cfg = cfg["confirmatory_hypotheses"]["H2_HIGH_DIMENSION_INSTABILITY_CONFIRMED"]
    h1_cells = []
    h2_cells = []
    for endpoint, mechanism in product(ENDPOINTS, MECHANISMS):
        full = selected_full[
            selected_full["endpoint"].eq(endpoint)
            & selected_full["mechanism"].eq(mechanism)
        ]
        low = selected_low[
            selected_low["endpoint"].eq(endpoint)
            & selected_low["mechanism"].eq(mechanism)
        ]
        require(len(full) == len(low) == 400, f"hypothesis aggregation multiplicity mismatch: {endpoint}/{mechanism}")
        full_keys = ["endpoint", "anchor", "mechanism", "seed", "model"]
        require(
            set(full[full_keys].itertuples(index=False, name=None))
            == set(low[full_keys].itertuples(index=False, name=None)),
            "H1/H2 selected row keys do not align",
        )
        b_oracle = float(np.median(np.abs(full["oracle_aipw"] - full["full_truth"])))
        b_full = float(np.median(np.abs(full["fitted_aipw"] - full["full_truth"])))
        b_low = float(np.median(np.abs(low["fitted_aipw"] - low["full_truth"])))
        h1_ratio = safe_ratio(b_oracle, b_full)
        h2_ratio = safe_ratio(b_low, b_full)
        support_valid = bool(full["fold_class_support"].all() and low["fold_class_support"].all())
        finite_valid = bool(full["all_finite"].all() and low["all_finite"].all())
        h1_adjudicative = bool(support_valid and finite_valid and h1_ratio is not None)
        h2_adjudicative = bool(support_valid and finite_valid and h2_ratio is not None)
        clip_reduction = float(np.median(full["fitted_clip_fraction"]) - np.median(low["fitted_clip_fraction"]))
        h1_cells.append(
            {
                "endpoint": endpoint,
                "mechanism": mechanism,
                "rows": 400,
                "B_oracle_aipw_fullP_lowM": b_oracle,
                "B_fitted_aipw_fullP_lowM": b_full,
                "ratio": h1_ratio,
                "threshold": float(h1_threshold),
                "denominator_valid": h1_ratio is not None,
                "fold_class_support": support_valid,
                "all_finite": finite_valid,
                "adjudicative": h1_adjudicative,
                "status": "OK" if h1_adjudicative else (
                    "NONADJUDICATIVE_FOLD_SUPPORT" if not support_valid else
                    "NONADJUDICATIVE_NONFINITE" if not finite_valid else
                    "NONADJUDICATIVE_DENOMINATOR"
                ),
                "pass": bool(h1_adjudicative and b_oracle <= h1_threshold * b_full),
            }
        )
        clip_full = float(np.median(full["fitted_clip_fraction"]))
        clip_low = float(np.median(low["fitted_clip_fraction"]))
        h2_cells.append(
            {
                "endpoint": endpoint,
                "mechanism": mechanism,
                "rows": 400,
                "B_fitted_aipw_lowP_lowM": b_low,
                "B_fitted_aipw_fullP_lowM": b_full,
                "error_ratio": h2_ratio,
                "error_threshold": float(h2_cfg["maximum_error_ratio"]),
                "denominator_valid": h2_ratio is not None,
                "fold_class_support": support_valid,
                "all_finite": finite_valid,
                "adjudicative": h2_adjudicative,
                "status": "OK" if h2_adjudicative else (
                    "NONADJUDICATIVE_FOLD_SUPPORT" if not support_valid else
                    "NONADJUDICATIVE_NONFINITE" if not finite_valid else
                    "NONADJUDICATIVE_DENOMINATOR"
                ),
                "clip_full_median": clip_full,
                "clip_low_median": clip_low,
                "clip_reduction": clip_reduction,
                "clip_threshold": float(h2_cfg["minimum_median_clip_fraction_reduction_full_minus_low"]),
                "error_pass": bool(h2_adjudicative and b_low <= h2_cfg["maximum_error_ratio"] * b_full),
                "clip_pass": bool(h2_adjudicative and clip_reduction >= h2_cfg["minimum_median_clip_fraction_reduction_full_minus_low"]),
                "pass": bool(
                    h2_adjudicative
                    and b_low <= h2_cfg["maximum_error_ratio"] * b_full
                    and clip_reduction >= h2_cfg["minimum_median_clip_fraction_reduction_full_minus_low"]
                ),
            }
        )
    return {
        "H1": {"pass": bool(all(row["pass"] for row in h1_cells)), "required_cells": 4, "cells": h1_cells},
        "H2": {"pass": bool(all(row["pass"] for row in h2_cells)), "required_cells": 4, "cells": h2_cells},
    }


def factorial_records(scenarios):
    records = []
    for endpoint, anchor, mechanism, model in product(
        ENDPOINTS, REGULAR_ANCHORS, MECHANISMS, MODELS
    ):
        part = scenarios[
            scenarios["endpoint"].eq(endpoint)
            & scenarios["anchor"].eq(anchor)
            & scenarios["mechanism"].eq(mechanism)
            & scenarios["model"].eq(model)
        ]
        require(len(part) == 400, "factorial stratum row count mismatch")
        cells = {}
        for p, m in product(P_REPRESENTATIONS, M_REPRESENTATIONS):
            cell = part[
                part["propensity_representation"].eq(p)
                & part["loss_representation"].eq(m)
            ]
            require(len(cell) == 100, "factorial cell seed count mismatch")
            cells[(p, m)] = float(np.median(np.abs(cell["fitted_aipw"] - cell["full_truth"])))
        low_low = cells[("lowdim_11", "lowdim_11")]
        low_full = cells[("lowdim_11", "full_2266")]
        full_low = cells[("full_2266", "lowdim_11")]
        full_full = cells[("full_2266", "full_2266")]
        values = {
            "P": 0.5 * ((full_low - low_low) + (full_full - low_full)),
            "M": 0.5 * ((low_full - low_low) + (full_full - full_low)),
            "P_x_M": (full_full - low_full) - (full_low - low_low),
        }
        records.extend(
            {
                "row_type": "factorial_contrast",
                "endpoint": endpoint,
                "anchor": anchor,
                "mechanism": mechanism,
                "model": model,
                "loss_representation": None,
                "contrast": contrast,
                "estimate": float(value),
                "numerator": None,
                "denominator": None,
                "ratio": None,
                "rows": 400,
                "status": "exploratory_only",
            }
            for contrast, value in values.items()
        )
    for endpoint, mechanism, m in product(ENDPOINTS, MECHANISMS, M_REPRESENTATIONS):
        part = scenarios[
            scenarios["endpoint"].eq(endpoint)
            & scenarios["mechanism"].eq(mechanism)
            & scenarios["propensity_representation"].eq("full_2266")
            & scenarios["loss_representation"].eq(m)
        ]
        require(len(part) == 400, "H3 aggregation row count mismatch")
        numerator = float(np.median(np.abs(part["oracle_aipw"] - part["full_truth"])))
        denominator = float(np.median(np.abs(part["oracle_ipw"] - part["full_truth"])))
        records.append(
            {
                "row_type": "augmentation_ratio",
                "endpoint": endpoint,
                "anchor": None,
                "mechanism": mechanism,
                "model": None,
                "loss_representation": m,
                "contrast": None,
                "estimate": None,
                "numerator": numerator,
                "denominator": denominator,
                "ratio": safe_ratio(numerator, denominator),
                "rows": 400,
                "status": "exploratory_only_no_gate",
            }
        )
    require(len(records) == 56, "factorial record count mismatch")
    return pd.DataFrame(records, columns=SCHEMAS["factorial_contrasts.csv"])


def bootstrap_seed_matrix(frame, column):
    ordered = frame.sort_values(["seed", "anchor", "model"])
    require(
        len(ordered) == 400
        and ordered.groupby("seed").size().eq(4).all()
        and ordered["seed"].drop_duplicates().tolist() == list(SEEDS),
        "bootstrap seed grid mismatch",
    )
    return ordered[column].to_numpy(float).reshape(100, 4)


def bootstrap_records(scenarios, cfg):
    bootstrap_cfg = cfg["confirmatory_hypotheses"]["bootstrap"]
    draws = np.random.default_rng(52001).integers(
        0, 100, size=(10000, 100), endpoint=False, dtype=np.int64
    )
    require(
        vector_hash(draws.ravel(order="C"), "int")
        == bootstrap_cfg["seed_index_matrix_sha256_numpy_1_26_4"],
        "bootstrap seed-index matrix hash mismatch",
    )
    records = []
    for endpoint, mechanism in product(ENDPOINTS, MECHANISMS):
        base = scenarios[
            scenarios["endpoint"].eq(endpoint)
            & scenarios["mechanism"].eq(mechanism)
            & scenarios["loss_representation"].eq("lowdim_11")
        ]
        full = base[base["propensity_representation"].eq("full_2266")]
        low = base[base["propensity_representation"].eq("lowdim_11")]
        require(len(full) == len(low) == 400, "bootstrap source multiplicity mismatch")
        support_valid = bool(full["fold_class_support"].all() and low["fold_class_support"].all())
        finite_valid = bool(full["all_finite"].all() and low["all_finite"].all())
        definitions = (
            (
                "H1_PROPENSITY_INSTABILITY_CONFIRMED",
                "oracle_to_fitted_aipw_error_ratio",
                bootstrap_seed_matrix(full, "oracle_aipw_abs_error"),
                bootstrap_seed_matrix(full, "fitted_aipw_abs_error"),
                float(cfg["confirmatory_hypotheses"]["H1_PROPENSITY_INSTABILITY_CONFIRMED"]["maximum_ratio"]),
                "less_or_equal",
                "ratio",
            ),
            (
                "H2_HIGH_DIMENSION_INSTABILITY_CONFIRMED",
                "low_to_full_fitted_aipw_error_ratio",
                bootstrap_seed_matrix(low, "fitted_aipw_abs_error"),
                bootstrap_seed_matrix(full, "fitted_aipw_abs_error"),
                float(cfg["confirmatory_hypotheses"]["H2_HIGH_DIMENSION_INSTABILITY_CONFIRMED"]["maximum_error_ratio"]),
                "less_or_equal",
                "ratio",
            ),
            (
                "H2_HIGH_DIMENSION_INSTABILITY_CONFIRMED",
                "full_minus_low_clip_fraction",
                bootstrap_seed_matrix(full, "fitted_clip_fraction"),
                bootstrap_seed_matrix(low, "fitted_clip_fraction"),
                float(cfg["confirmatory_hypotheses"]["H2_HIGH_DIMENSION_INSTABILITY_CONFIRMED"]["minimum_median_clip_fraction_reduction_full_minus_low"]),
                "greater_or_equal",
                "difference",
            ),
        )
        for hypothesis, statistic, numerator, denominator, threshold, direction, operation in definitions:
            numerator_point = float(np.median(numerator))
            denominator_point = float(np.median(denominator))
            numerator_boot = np.median(numerator[draws], axis=(1, 2))
            denominator_boot = np.median(denominator[draws], axis=(1, 2))
            require(np.isfinite(numerator_point) and np.isfinite(numerator_boot).all(), "nonfinite bootstrap numerator")
            if operation == "ratio":
                denominator_valid = bool(
                    np.isfinite(denominator_point)
                    and denominator_point > 0
                    and np.isfinite(denominator_boot).all()
                    and np.all(denominator_boot > 0)
                )
                if denominator_valid:
                    values = numerator_boot / denominator_boot
                    point = numerator_point / denominator_point
                    lower, upper = np.quantile(values, [0.025, 0.975], method="linear")
                    status = "OK"
                else:
                    point = lower = upper = 0.0
                    status = "NONADJUDICATIVE_DENOMINATOR"
            else:
                denominator_valid = bool(np.isfinite(denominator_point) and np.isfinite(denominator_boot).all())
                values = numerator_boot - denominator_boot
                require(denominator_valid and np.isfinite(values).all(), "nonfinite bootstrap difference")
                point = numerator_point - denominator_point
                lower, upper = np.quantile(values, [0.025, 0.975], method="linear")
                status = "OK"
            if not finite_valid:
                status = "NONADJUDICATIVE_NONFINITE"
            elif not support_valid:
                status = "NONADJUDICATIVE_FOLD_SUPPORT"
            records.append(
                {
                    "hypothesis": hypothesis,
                    "endpoint": endpoint,
                    "mechanism": mechanism,
                    "statistic": statistic,
                    "point_estimate": float(point),
                    "lower": float(lower),
                    "upper": float(upper),
                    "threshold": threshold,
                    "direction": direction,
                    "point_pass": bool(
                        denominator_valid and support_valid and finite_valid
                        and (point <= threshold if direction == "less_or_equal" else point >= threshold)
                    ),
                    "denominator_valid": denominator_valid,
                    "resamples": 10000,
                    "seed_index_matrix_sha256": bootstrap_cfg["seed_index_matrix_sha256_numpy_1_26_4"],
                    "status": status,
                }
            )
    require(len(records) == 12, "bootstrap record count mismatch")
    return pd.DataFrame(records, columns=SCHEMAS["bootstrap_intervals.csv"])


def admission_result(scenarios, cfg):
    candidate = scenarios[
        scenarios["propensity_representation"].eq("lowdim_11")
        & scenarios["loss_representation"].eq("lowdim_11")
    ]
    criteria = cfg["admission_gate"]["criteria"]
    cells = []
    for endpoint, anchor, model, mechanism in product(
        ENDPOINTS, REGULAR_ANCHORS, MODELS, MECHANISMS
    ):
        part = candidate[
            candidate["endpoint"].eq(endpoint)
            & candidate["anchor"].eq(anchor)
            & candidate["model"].eq(model)
            & candidate["mechanism"].eq(mechanism)
        ]
        require(len(part) == 100, "admission cell seed count mismatch")
        b_aipw = float(np.median(np.abs(part["fitted_aipw"] - part["full_truth"])))
        b_cc = float(np.median(np.abs(part["cc"] - part["full_truth"])))
        ratio = safe_ratio(b_aipw, b_cc)
        denominator_valid = bool(np.isfinite(b_cc) and b_cc > 0)
        coverage = float(part["fitted_covered"].astype(float).mean())
        ess_q05 = quantile(part["fitted_ess_fraction"], 0.05)
        clip_median = float(np.median(part["fitted_clip_fraction"]))
        clip_q95 = quantile(part["fitted_clip_fraction"], 0.95)
        smd_q95 = quantile(part["fitted_max_abs_smd"], 0.95)
        support = bool(part["fold_class_support"].astype(bool).all())
        all_finite = bool(part["all_finite"].astype(bool).all())
        ratio_limit = (
            criteria["maximum_median_aipw_to_complete_case_absolute_estimation_error_ratio_mcar"]
            if mechanism == "mcar"
            else criteria["maximum_median_aipw_to_complete_case_absolute_estimation_error_ratio_mar"]
        )
        passed = bool(
            ratio is not None
            and ratio <= ratio_limit
            and coverage >= criteria["minimum_empirical_95pct_coverage"]
            and ess_q05 >= criteria["minimum_fifth_percentile_ess_over_observed"]
            and clip_median <= criteria["maximum_median_clip_fraction"]
            and clip_q95 <= criteria["maximum_95th_percentile_clip_fraction"]
            and smd_q95 <= criteria["maximum_95th_percentile_max_abs_smd"]
            and support
            and all_finite
        )
        cells.append(
            {
                "endpoint": endpoint,
                "anchor": anchor,
                "model": model,
                "mechanism": mechanism,
                "rows": 100,
                "B_CC": b_cc,
                "B_AIPW": b_aipw,
                "error_ratio": ratio,
                "error_threshold": float(ratio_limit),
                "denominator_valid": denominator_valid,
                "coverage": coverage,
                "ess_fraction_q05": ess_q05,
                "clip_fraction_median": clip_median,
                "clip_fraction_q95": clip_q95,
                "max_abs_smd_q95": smd_q95,
                "fold_class_support": support,
                "all_finite": all_finite,
                "pass": passed,
            }
        )
    require(len(cells) == 16, "admission result cell count mismatch")
    return {
        "study_id": "confirm_001",
        "candidate_specification": {
            "propensity_representation": "lowdim_11",
            "loss_representation": "lowdim_11",
        },
        "known_truth_admission": {
            "pass": bool(all(row["pass"] for row in cells)),
            "required_cells": 16,
            "cells": cells,
        },
    }


def direct_result(actual, admission, cfg):
    candidate = actual[
        actual["propensity_representation"].eq("lowdim_11")
        & actual["loss_representation"].eq("lowdim_11")
        & actual["anchor"].isin(REGULAR_ANCHORS)
    ]
    require(len(candidate) == 8, "direct candidate grid mismatch")
    overlap_cfg = cfg["direct_mask_overlap_gate"]
    cells = []
    for row in candidate.sort_values(["endpoint", "anchor", "model"]).to_dict("records"):
        minimum_observed = min(int(row[f"fold_{fold}_train_observed"]) for fold in range(3))
        minimum_unavailable = min(int(row[f"fold_{fold}_train_unavailable"]) for fold in range(3))
        eligible = bool(
            minimum_observed >= overlap_cfg["minimum_observed_and_unavailable_per_fold"]
            and minimum_unavailable >= overlap_cfg["minimum_observed_and_unavailable_per_fold"]
            and float(row["fitted_clip_fraction"]) <= overlap_cfg["maximum_clip_fraction"]
            and float(row["fitted_ess_fraction"]) >= overlap_cfg["minimum_ess_over_observed"]
            and float(row["fitted_max_abs_smd"]) <= overlap_cfg["maximum_max_abs_smd"]
            and float(row["fitted_max_normalized_weight"])
            <= overlap_cfg["maximum_normalized_observed_weight"]
            and bool(row["fold_class_support"])
            and bool(row["all_finite"])
        )
        cells.append(
            {
                "endpoint": row["endpoint"],
                "anchor": row["anchor"],
                "model": row["model"],
                "fold_class_support": bool(row["fold_class_support"]),
                "clip_fraction": float(row["fitted_clip_fraction"]),
                "ess_fraction": float(row["fitted_ess_fraction"]),
                "max_abs_smd": float(row["fitted_max_abs_smd"]),
                "max_normalized_weight": float(row["fitted_max_normalized_weight"]),
                "pass": eligible,
                "label": "ELIGIBLE_OVERLAP" if eligible else overlap_cfg["failure_label"],
            }
        )
    all_overlap = bool(all(row["pass"] for row in cells))
    b_aipw = float(np.median(np.abs(candidate["fitted_aipw"] - candidate["full_truth"])))
    b_cc = float(np.median(np.abs(candidate["cc"] - candidate["full_truth"])))
    ratio = safe_ratio(b_aipw, b_cc)
    reduction = None if ratio is None else float(1 - ratio)
    improved = int(
        (np.abs(candidate["fitted_aipw"] - candidate["full_truth"]) < np.abs(candidate["cc"] - candidate["full_truth"])).sum()
    )
    translation_cfg = cfg["direct_mask_translation_gate"]
    if not admission["known_truth_admission"]["pass"]:
        status = translation_cfg["admission_failure_label"]
        passed = False
    elif not all_overlap:
        status = "NOT_ADJUDICATED_OVERLAP_FAILED"
        passed = False
    elif reduction is None:
        status = "NONADJUDICATIVE_FAIL"
        passed = False
    else:
        passed = bool(
            reduction >= translation_cfg["minimum_pooled_median_absolute_estimation_error_reduction_vs_complete_case"]
            and improved >= translation_cfg["minimum_cells_with_lower_aipw_absolute_estimation_error"]
        )
        status = "PASS" if passed else "FAIL"
    return {
        "direct_mask_overlap": {"pass": all_overlap, "required_cells": 8, "cells": cells},
        "direct_mask_translation": {
            "status": status,
            "pass": passed,
            "required_cells": 8,
            "cc_median_abs_error": b_cc,
            "aipw_median_abs_error": b_aipw,
            "denominator_valid": bool(ratio is not None),
            "error_reduction": reduction,
            "cells_aipw_better": improved,
            "required_aipw_better": int(
                translation_cfg["minimum_cells_with_lower_aipw_absolute_estimation_error"]
            ),
        },
    }


def verify_pxr_linkage(stored, context, cfg, tolerance):
    required = cfg["expected_outputs"]["pxr_linkage_audit.csv"]["required_fields"]
    require_columns(stored, required, "PXR linkage audit")
    require(len(stored) == 10747, "PXR linkage row count mismatch")
    require(not stored.duplicated(["source_row_index"]).any(), "duplicate PXR source row")
    stored = stored.sort_values("source_row_index").reset_index(drop=True)
    expected = context["expected_linkage"]
    for column in expected.columns:
        if column in ("exact_smiles_match_if_linked", "batch_match_if_linked"):
            left_missing = stored[column].isna().to_numpy()
            right_missing = expected[column].isna().to_numpy()
            require(np.array_equal(left_missing, right_missing), f"PXR {column} missingness mismatch")
            def semantic_bool(values):
                mapped = []
                for value in values:
                    if isinstance(value, (bool, np.bool_)):
                        mapped.append(bool(value))
                    else:
                        text_value = str(value).strip().lower()
                        require(text_value in ("true", "false"), f"PXR {column} invalid Boolean value")
                        mapped.append(text_value == "true")
                return np.asarray(mapped, dtype=bool)
            require(
                np.array_equal(
                    semantic_bool(stored.loc[~left_missing, column].to_numpy()),
                    semantic_bool(expected.loc[~right_missing, column].to_numpy()),
                ),
                f"PXR {column} mismatch",
            )
        elif pd.api.types.is_numeric_dtype(expected[column]):
            same_vector(stored[column].to_numpy(), expected[column].to_numpy(), tolerance, f"PXR {column}")
        else:
            same_vector(stored[column].astype(str), expected[column].astype(str), tolerance, f"PXR {column}")
    require(
        np.array_equal(stored["cohort_row_index"].to_numpy(np.int64), np.arange(10747)),
        "PXR cohort row order mismatch",
    )

    summary_rows = []
    thresholds = cfg["pxr_boundary"]["abstention_triggers_any"]
    for specification in ("pxr_structure_10", "pxr_structure_screen_12"):
        suffix = specification.removeprefix("pxr_")
        raw_column = f"p_raw_{suffix}"
        clipped_column = f"p_clipped_{suffix}"
        indicator_column = f"clipped_{suffix}"
        raw = stored[raw_column].to_numpy(float)
        clipped = stored[clipped_column].to_numpy(float)
        indicator = stored[indicator_column].to_numpy(np.uint8)
        require(np.isfinite(raw).all() and np.all((raw >= 0) & (raw <= 1)), f"invalid PXR raw propensity: {specification}")
        same_vector(clipped, np.clip(raw, 0.05, 0.95), tolerance, f"PXR clipped propensity: {specification}")
        require(
            np.array_equal(indicator, ((raw < 0.05) | (raw > 0.95)).astype(np.uint8)),
            f"PXR strict clipping indicator mismatch: {specification}",
        )
        mask = context["available"]
        weights = mask / clipped
        ess = float(weights.sum() ** 2 / np.square(weights).sum())
        ess_fraction = ess / float(mask.sum())
        maximum_weight = float(weights.max() / weights.sum())
        feature_names = (
            PXR_STRUCTURE_COLUMNS
            if specification == "pxr_structure_10"
            else (*PXR_STRUCTURE_COLUMNS, "log2_fc_estimate", "fdr_score")
        )
        smd = smd_map(mask, clipped, context["features"][specification], feature_names)
        zero_variance = [
            name
            for name, standard_deviation in zip(
                feature_names, context["features"][specification].std(axis=0, ddof=0)
            )
            if standard_deviation == 0
        ]
        max_smd = max(smd.values())
        clip_fraction = float(indicator.mean())
        clip_trigger = bool(clip_fraction > thresholds["maximum_clip_fraction"])
        ess_trigger = bool(ess_fraction < thresholds["minimum_ess_over_observed"])
        smd_trigger = bool(max_smd > thresholds["maximum_max_abs_smd"])
        summary_rows.append(
            {
                "specification": specification,
                "n_total": 10747,
                "n_observed": 2722,
                "n_unavailable": 8025,
                "coverage": float(mask.mean()),
                "clip_fraction": clip_fraction,
                "ess": ess,
                "ess_fraction": ess_fraction,
                "max_abs_smd": max_smd,
                "max_normalized_weight": maximum_weight,
                "smd_json": json.dumps(smd, sort_keys=True),
                "zero_variance_features_json": json.dumps(zero_variance),
                **fold_support(mask, context["fold"]),
                "fold_class_support": min(fold_support(mask, context["fold"]).values()) >= 20,
                "all_finite": True,
                "source_row_index_sha256": context["hashes"]["ordered_source_row_index"],
                "murcko_group_vector_sha256": context["hashes"]["ordered_murcko_group"],
                "fold_sha256": context["hashes"]["ordered_fold"],
                "p_raw_sha256": vector_hash(raw, "float"),
                "p_clipped_sha256": vector_hash(clipped, "float"),
                "availability_sha256": context["hashes"]["ordered_availability"],
                "formal_hit_sha256": context["hashes"]["ordered_formal_hit"],
                "trigger_clip": clip_trigger,
                "trigger_ess": ess_trigger,
                "trigger_smd": smd_trigger,
                "overlap_pass": bool(not (clip_trigger or ess_trigger or smd_trigger) and min(fold_support(mask, context["fold"]).values()) >= 20),
            }
        )
    return stored, pd.DataFrame(summary_rows, columns=SCHEMAS["pxr_overlap_summary.csv"])


def pxr_strata_and_boundary(context, overlap, cfg):
    rows = []
    minimum = int(cfg["pxr_boundary"]["abstention_triggers_any"]["selection_stratum_minimum_n"])
    lower, upper = map(
        float,
        cfg["pxr_boundary"]["abstention_triggers_any"]["selection_stratum_extreme_probability"],
    )
    for value in (1, 0):
        selected = context["formal_hit"] == value
        total = int(selected.sum())
        observed = int(context["available"][selected].sum())
        probability = observed / total
        trigger = bool(total >= minimum and (probability <= lower or probability >= upper))
        rows.append(
            {
                "formal_hit": value,
                "n_total": total,
                "n_followup_available": observed,
                "n_followup_unavailable": total - observed,
                "advancement_probability": probability,
            }
        )
    strata = pd.DataFrame(rows, columns=SCHEMAS["pxr_selection_strata.csv"])
    require(
        list(strata[["n_total", "n_followup_available", "n_followup_unavailable"]].itertuples(index=False, name=None))
        == [(873, 869, 4), (9874, 1853, 8021)],
        "PXR selection-stratum count mismatch",
    )
    selection_extreme = bool(
        (
            (strata["n_total"] >= minimum)
            & ((strata["advancement_probability"] <= lower) | (strata["advancement_probability"] >= upper))
        ).any()
    )
    boundary = {
        "study_id": "confirm_001",
        "decision": cfg["pxr_boundary"]["abstention_label"],
        "selection_strata_extreme": selection_extreme,
        "propensity_triggers": [
            {
                "specification": row["specification"],
                "clip": bool(row["trigger_clip"]),
                "ess": bool(row["trigger_ess"]),
                "smd": bool(row["trigger_smd"]),
                "any": bool(row["trigger_clip"] or row["trigger_ess"] or row["trigger_smd"]),
            }
            for row in overlap.to_dict("records")
        ],
        "documentary_triggers": {
            "missing_documented_promotion_branch_or_selection_variable": True,
            "public_files_do_not_reconstruct_official_parent_and_promoted_populations": True,
            "missing_nonadvanced_dose_response_truth_without_randomized_audit_or_identification_design": True,
        },
        "full_library_point_estimate_allowed": False,
        "allowed_output": cfg["pxr_boundary"]["allowed_output"],
    }
    return strata, boundary


def read_table(path):
    return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path, float_precision="round_trip")


def require_schema(frame, filename):
    require(filename in SCHEMAS, f"no frozen schema for {filename}")
    require(tuple(frame.columns) == SCHEMAS[filename], f"exact schema/order mismatch: {filename}")


def compare_frames(expected, observed, filename, keys, tolerance):
    require_schema(observed, filename)
    require(tuple(expected.columns) == SCHEMAS[filename], f"internal expected schema mismatch: {filename}")
    require(len(expected) == len(observed), f"row count mismatch: {filename}")
    left = expected.copy()
    right = observed.copy()
    left["__sort_key"] = left[list(keys)].fillna("<NULL>").astype(str).agg("\x1f".join, axis=1)
    right["__sort_key"] = right[list(keys)].fillna("<NULL>").astype(str).agg("\x1f".join, axis=1)
    require(not left["__sort_key"].duplicated().any(), f"duplicate expected key: {filename}")
    require(not right["__sort_key"].duplicated().any(), f"duplicate stored key: {filename}")
    left = left.sort_values("__sort_key").reset_index(drop=True)
    right = right.sort_values("__sort_key").reset_index(drop=True)
    require(np.array_equal(left["__sort_key"], right["__sort_key"]), f"key grid mismatch: {filename}")
    for column in expected.columns:
        if pd.api.types.is_numeric_dtype(left[column]) and pd.api.types.is_numeric_dtype(right[column]):
            require(
                np.allclose(
                    left[column].to_numpy(float), right[column].to_numpy(float),
                    rtol=0, atol=tolerance, equal_nan=True,
                ),
                f"numeric mismatch: {filename}/{column}",
            )
        else:
            require(
                np.array_equal(
                    left[column].fillna("<NULL>").astype(str),
                    right[column].fillna("<NULL>").astype(str),
                ),
                f"value mismatch: {filename}/{column}",
            )


def verify_availability_map(frame, contexts, tolerance):
    require_schema(frame, "availability_map.csv")
    require(len(frame) == 8, "availability-map row count mismatch")
    require_exact_grid(frame, ("endpoint", "anchor"), product(ENDPOINTS, REGULAR_ANCHORS + STRESS_ANCHORS), "availability")
    for row in frame.to_dict("records"):
        mask = contexts[row["endpoint"]]["availability"][row["anchor"]]
        require(int(row["n_total"]) == len(mask), "availability n_total mismatch")
        require(int(row["n_available"]) == int(mask.sum()), "availability n_available mismatch")
        require(int(row["n_unavailable"]) == len(mask) - int(mask.sum()), "availability n_unavailable mismatch")
        close(row["coverage"], mask.mean(), tolerance, "availability coverage")
        require(row["availability_sha256"] == vector_hash(mask, "binary"), "availability hash mismatch")


def verify_cohort_audit(frame, data, tolerance):
    require_schema(frame, "cohort_audit.csv")
    endpoints = ENDPOINTS + REGULAR_ANCHORS + STRESS_ANCHORS
    require(len(frame) == 18, "cohort-audit row count mismatch")
    require_exact_grid(
        frame,
        ("dataset", "split", "endpoint"),
        product(("Biogen_public_ADME",), ("all", "train", "test"), endpoints),
        "cohort audit",
    )
    for row in frame.to_dict("records"):
        split = np.ones(len(data), dtype=bool) if row["split"] == "all" else data["ds"].eq(row["split"]).to_numpy()
        finite = np.isfinite(pd.to_numeric(data[row["endpoint"]], errors="coerce").to_numpy(float)) & split
        n_rows = int(split.sum())
        n_finite = int(finite.sum())
        require(int(row["n_rows"]) == n_rows, "cohort n_rows mismatch")
        require(int(row["n_finite"]) == n_finite, "cohort n_finite mismatch")
        require(int(row["n_unavailable"]) == n_rows - n_finite, "cohort n_unavailable mismatch")
        close(row["coverage"], n_finite / n_rows, tolerance, "cohort coverage")


def expected_structure_audit(structure):
    return {
        "rows": 3521,
        "smiles_parse_failures": 0,
        "unique_canonical_smiles": 3521,
        "canonical_train_test_overlap": 0,
        "ordered_canonical_smiles_sha256": structure["ordered_canonical_smiles_sha256"],
        "descriptor_count": 217,
        "descriptor_names_sha256": vector_hash(structure["descriptor_names"], "text"),
        "descriptor_names": structure["descriptor_names"],
        "all_descriptors_finite": True,
    }


def verify_predictions(frame, data, contexts, model_features, cfg, tolerance):
    require_schema(frame, "predictions.parquet")
    require(len(frame) == 3702, "prediction row count mismatch")
    expected = {
        (endpoint, model, int(source))
        for endpoint in ENDPOINTS
        for model in MODELS
        for source in contexts[endpoint]["source_row_index"]
    }
    require_exact_grid(frame, ("endpoint", "model", "source_row_index"), expected, "prediction")
    losses = {}
    for endpoint, model in product(ENDPOINTS, MODELS):
        context = contexts[endpoint]
        part = frame[frame["endpoint"].eq(endpoint) & frame["model"].eq(model)].sort_values("source_row_index")
        same_vector(part["source_row_index"], context["source_row_index"], tolerance, "prediction source rows")
        same_vector(part["cohort_row_index"], context["cohort_row_index"], tolerance, "prediction cohort rows")
        same_vector(part["molecule_id"].astype(str), context["name"], tolerance, "prediction molecule IDs")
        same_vector(part["reference"], context["y"], tolerance, "prediction reference")
        prediction = part["prediction"].to_numpy(float)
        require(np.isfinite(prediction).all(), "nonfinite prediction")
        require(np.ptp(prediction) > 0, f"constant prediction: {endpoint}/{model}")
        all_y = pd.to_numeric(data[endpoint], errors="coerce").to_numpy(float)
        train = np.flatnonzero(data["ds"].eq("train").to_numpy() & np.isfinite(all_y))
        target = context["source_row_index"]
        if model == "ecfp_rf":
            model_cfg = cfg["prediction_models"][model]
            fitted_model = RandomForestRegressor(
                n_estimators=int(model_cfg["n_estimators"]),
                criterion=model_cfg["criterion"],
                min_samples_leaf=int(model_cfg["min_samples_leaf"]),
                max_features=model_cfg["max_features"],
                n_jobs=int(model_cfg["n_jobs"]),
                random_state=int(model_cfg["random_state"]),
            ).fit(model_features["fingerprint_matrix"][train], all_y[train])
            independently_fitted = fitted_model.predict(model_features["fingerprint_matrix"][target])
        else:
            model_cfg = {
                key: value
                for key, value in cfg["prediction_models"][model].items()
                if key != "features"
            }
            full = np.hstack(
                (model_features["fingerprint_matrix"], model_features["descriptor_matrix"])
            ).astype(np.float32)
            fitted_model = XGBRegressor(**model_cfg).fit(full[train], all_y[train])
            independently_fitted = fitted_model.predict(full[target])
        same_vector(prediction, independently_fitted, tolerance, f"independent prediction refit: {endpoint}/{model}")
        loss = np.abs(context["y"] - prediction)
        same_vector(part["loss"], loss, tolerance, "stored prediction loss")
        losses[(endpoint, model)] = loss
    return losses


def verify_folds(frame, contexts, tolerance):
    require_schema(frame, "nuisance_fold_assignments.parquet")
    require(len(frame) == 1851, "fold-assignment row count mismatch")
    expected = {
        (endpoint, int(source))
        for endpoint in ENDPOINTS
        for source in contexts[endpoint]["source_row_index"]
    }
    require_exact_grid(frame, ("endpoint", "source_row_index"), expected, "fold")
    for endpoint in ENDPOINTS:
        context = contexts[endpoint]
        part = frame[frame["endpoint"].eq(endpoint)].sort_values("source_row_index")
        same_vector(part["source_row_index"], context["source_row_index"], tolerance, "fold source rows")
        same_vector(part["cohort_row_index"], context["cohort_row_index"], tolerance, "fold cohort rows")
        same_vector(part["molecule_id"].astype(str), context["name"], tolerance, "fold molecule IDs")
        same_vector(part["scaffold_group"].astype(str), context["groups"], tolerance, "fold scaffold groups")
        same_vector(part["scaffold_sha256"].astype(str), [row_text_hash(value) for value in context["groups"]], tolerance, "fold row scaffold hashes")
        same_vector(part["fold"], context["fold"], tolerance, "fold IDs")
        for column, wanted in (
            ("source_row_index_sha256", context["source_hash"]),
            ("scaffold_group_vector_sha256", context["group_hash"]),
            ("fold_sha256", context["fold_hash"]),
        ):
            require(part[column].nunique() == 1 and part[column].iloc[0] == wanted, f"{endpoint} {column} mismatch")


def verify_generators(frame, contexts, contracts, tolerance):
    require_schema(frame, "generator_vectors.parquet")
    require(len(frame) == 7404, "generator row count mismatch")
    expected = {
        (endpoint, anchor, mechanism, int(source))
        for endpoint in ENDPOINTS
        for anchor in REGULAR_ANCHORS
        for mechanism in MECHANISMS
        for source in contexts[endpoint]["source_row_index"]
    }
    require_exact_grid(frame, ("endpoint", "anchor", "mechanism", "source_row_index"), expected, "generator")
    for endpoint, anchor, mechanism in product(ENDPOINTS, REGULAR_ANCHORS, MECHANISMS):
        context = contexts[endpoint]
        contract = contracts[(endpoint, anchor, mechanism)]
        part = frame[
            frame["endpoint"].eq(endpoint)
            & frame["anchor"].eq(anchor)
            & frame["mechanism"].eq(mechanism)
        ].sort_values("source_row_index")
        same_vector(part["source_row_index"], context["source_row_index"], tolerance, "generator source rows")
        same_vector(part["cohort_row_index"], context["cohort_row_index"], tolerance, "generator cohort rows")
        same_vector(part["p_oracle"], contract["probability"], tolerance, "generator oracle probability")
        close(part["direct_coverage"].iloc[0], contract["coverage"], tolerance, "generator direct coverage")
        require(part["direct_coverage"].nunique() == 1, "generator direct coverage not constant")
        require(part["p_oracle_sha256"].nunique() == 1 and part["p_oracle_sha256"].iloc[0] == contract["probability_hash"], "generator probability hash mismatch")
        close(part["coverage_error"].iloc[0], abs(np.mean(contract["probability"]) - contract["coverage"]), tolerance, "generator coverage error")
        require(part["coverage_error"].nunique() == 1, "generator coverage error not constant")
        require(part["all_finite"].astype(bool).all(), "generator finite flag failed")
        finite_columns(part, ("direct_coverage", "generator_score", "intercept_shift", "p_oracle", "coverage_error"), "generator values")
        if mechanism == "mar":
            same_vector(part["generator_score"], contract["generator_score"], tolerance, "MAR generator score")
            same_vector(
                part["intercept_shift"],
                np.full(len(part), contract["intercept_shift"], dtype=float),
                tolerance,
                "MAR intercept shift",
            )
        else:
            require(np.allclose(part["generator_score"], 0, rtol=0, atol=tolerance), "MCAR generator score must be zero")
            require(np.allclose(part["intercept_shift"], 0, rtol=0, atol=tolerance), "MCAR intercept shift must be zero")


def fold_support(mask, folds):
    result = {}
    mask = np.asarray(mask, dtype=np.uint8)
    folds = np.asarray(folds, dtype=np.int64)
    for fold in range(3):
        training = folds != fold
        result[f"fold_{fold}_train_observed"] = int(mask[training].sum())
        result[f"fold_{fold}_train_unavailable"] = int(training.sum() - mask[training].sum())
    return result


def verify_scenario_grid(frame, contexts, contracts, losses, cfg, tolerance):
    require_schema(frame, "semisynthetic_scenarios.csv")
    require(len(frame) == 6400, "semisynthetic row count mismatch")
    keys = ("endpoint", "anchor", "mechanism", "seed", "model", "propensity_representation", "loss_representation")
    require_exact_grid(
        frame,
        keys,
        product(ENDPOINTS, REGULAR_ANCHORS, MECHANISMS, SEEDS, MODELS, P_REPRESENTATIONS, M_REPRESENTATIONS),
        "semisynthetic scenario",
    )
    string_columns = {
        *keys[:-1], "loss_representation", "source_row_index_sha256", "scaffold_group_vector_sha256", "fold_sha256",
        "R_sha256", "p_oracle_sha256", "p_raw_sha256", "p_fitted_sha256", "mhat_sha256",
        "fitted_smd_json", "oracle_smd_json",
    }
    finite_columns(frame, [column for column in frame.columns if column not in string_columns], "semisynthetic scenario")
    require(frame["seed"].astype(int).isin(SEEDS).all(), "unexpected semisynthetic seed")
    for row in frame.to_dict("records"):
        endpoint, anchor, mechanism, seed, model, p_rep, m_rep = (row[key] for key in keys)
        context = contexts[endpoint]
        contract = contracts[(endpoint, anchor, mechanism)]
        mask, mask_hash = regenerated_mask(context, contract, int(seed))
        loss = losses[(endpoint, model)]
        support = fold_support(mask, context["fold"])
        require(int(row["n_total"]) == len(mask), "scenario n_total mismatch")
        require(int(row["n_observed"]) == int(mask.sum()), "scenario n_observed mismatch")
        close(row["observed_fraction"], mask.mean(), tolerance, "scenario observed fraction")
        close(row["direct_anchor_coverage"], contract["coverage"], tolerance, "scenario direct coverage")
        expected_shift = contract["intercept_shift"] if mechanism == "mar" else 0.0
        close(row["generator_intercept_shift"], expected_shift, tolerance, "scenario generator shift")
        require(row["source_row_index_sha256"] == context["source_hash"], "scenario source hash mismatch")
        require(row["scaffold_group_vector_sha256"] == context["group_hash"], "scenario scaffold hash mismatch")
        require(row["fold_sha256"] == context["fold_hash"], "scenario fold hash mismatch")
        require(row["R_sha256"] == mask_hash, "scenario mask hash mismatch")
        require(row["p_oracle_sha256"] == contract["probability_hash"], "scenario oracle-p hash mismatch")
        close(row["full_truth"], loss.mean(), tolerance, "scenario truth")
        close(row["cc"], loss[mask.astype(bool)].mean(), tolerance, "scenario complete case")
        close(row["cc_abs_error"], abs(row["cc"] - row["full_truth"]), tolerance, "scenario CC error")
        close(row["fitted_aipw_abs_error"], abs(row["fitted_aipw"] - row["full_truth"]), tolerance, "scenario fitted error")
        close(row["oracle_ipw_abs_error"], abs(row["oracle_ipw"] - row["full_truth"]), tolerance, "scenario oracle IPW error")
        close(row["oracle_aipw_abs_error"], abs(row["oracle_aipw"] - row["full_truth"]), tolerance, "scenario oracle AIPW error")
        for prefix in ("fitted", "oracle"):
            close(row[f"{prefix}_lcl"], row[f"{prefix}_aipw"] - 1.96 * row[f"{prefix}_se"], tolerance, f"{prefix} lower interval")
            close(row[f"{prefix}_ucl"], row[f"{prefix}_aipw"] + 1.96 * row[f"{prefix}_se"], tolerance, f"{prefix} upper interval")
            require(bool(row[f"{prefix}_covered"]) == bool(row[f"{prefix}_lcl"] <= row["full_truth"] <= row[f"{prefix}_ucl"]), f"{prefix} coverage mismatch")
        oracle_weights = mask / contract["probability"]
        oracle_ess = float(oracle_weights.sum() ** 2 / np.square(oracle_weights).sum())
        close(row["oracle_ess"], oracle_ess, tolerance, "scenario oracle ESS")
        close(row["oracle_ess_fraction"], oracle_ess / mask.sum(), tolerance, "scenario oracle ESS fraction")
        oracle_ipw = float(np.sum(oracle_weights * loss) / oracle_weights.sum())
        close(row["oracle_ipw"], oracle_ipw, tolerance, "scenario oracle IPW")
        oracle_smd = smd_map(mask, contract["probability"], context["lowdim"])
        stored_oracle_smd = json.loads(row["oracle_smd_json"])
        compare_tree(oracle_smd, stored_oracle_smd, tolerance, "scenario.oracle_smd")
        close(row["oracle_max_abs_smd"], max(oracle_smd.values()), tolerance, "scenario oracle max SMD")
        close(row["oracle_max_normalized_weight"], oracle_weights.max() / oracle_weights.sum(), tolerance, "scenario oracle max weight")
        stored_fitted_smd = json.loads(row["fitted_smd_json"])
        require(set(stored_fitted_smd) == set(LOWDIM_COLUMNS), "fitted SMD feature grid mismatch")
        require(np.isfinite(np.asarray(list(stored_fitted_smd.values()), dtype=float)).all(), "nonfinite fitted SMD")
        close(row["fitted_max_abs_smd"], max(stored_fitted_smd.values()), tolerance, "scenario fitted max SMD")
        for column, value in support.items():
            require(int(row[column]) == value, f"scenario {column} mismatch")
        expected_support = min(support.values()) >= 20
        require(bool(row["fold_class_support"]) == expected_support, "scenario fold-support flag mismatch")
        require(bool(row["all_finite"]), "scenario all_finite flag failed")
        require(-tolerance <= float(row["fitted_ess_fraction"]) <= 1 + tolerance, "invalid fitted ESS fraction")
        require(-tolerance <= float(row["fitted_clip_fraction"]) <= 1 + tolerance, "invalid fitted clip fraction")
        require(0 < float(row["fitted_max_normalized_weight"]) <= 1 + tolerance, "invalid fitted maximum weight")
    require(
        frame.groupby(["endpoint", "anchor", "mechanism", "seed", "propensity_representation"])[
            ["p_raw_sha256", "p_fitted_sha256", "fitted_ess", "fitted_ess_fraction", "fitted_clip_fraction"]
        ].nunique().max().max() == 1,
        "propensity outputs vary across model or loss representation",
    )
    require(
        frame.groupby(["endpoint", "anchor", "mechanism", "seed", "model", "loss_representation"])["mhat_sha256"].nunique().max() == 1,
        "loss nuisance varies across propensity representation",
    )
    return hypothesis_details(frame, cfg)


def verify_actual_grid(frame, contexts, losses, cfg, tolerance):
    require_schema(frame, "actual_mask_scenarios.csv")
    require(len(frame) == 64, "actual-mask row count mismatch")
    keys = ("endpoint", "anchor", "model", "propensity_representation", "loss_representation")
    require_exact_grid(
        frame,
        keys,
        product(ENDPOINTS, REGULAR_ANCHORS + STRESS_ANCHORS, MODELS, P_REPRESENTATIONS, M_REPRESENTATIONS),
        "actual mask scenario",
    )
    require(frame["mask_type"].eq("direct_empirical").all(), "actual mask_type mismatch")
    string_columns = {
        *keys, "mask_type", "source_row_index_sha256", "scaffold_group_vector_sha256", "fold_sha256", "R_sha256",
        "p_raw_sha256", "p_fitted_sha256", "mhat_sha256", "fitted_smd_json", "overlap_label",
    }
    finite_columns(frame, [column for column in frame.columns if column not in string_columns], "actual-mask scenario")
    overlap_cfg = cfg["direct_mask_overlap_gate"]
    for row in frame.to_dict("records"):
        endpoint, anchor, model, _, _ = (row[key] for key in keys)
        context = contexts[endpoint]
        mask = context["availability"][anchor]
        loss = losses[(endpoint, model)]
        support = fold_support(mask, context["fold"])
        require(int(row["n_total"]) == len(mask), "actual n_total mismatch")
        require(int(row["n_observed"]) == int(mask.sum()), "actual n_observed mismatch")
        close(row["observed_fraction"], mask.mean(), tolerance, "actual observed fraction")
        require(row["source_row_index_sha256"] == context["source_hash"], "actual source hash mismatch")
        require(row["scaffold_group_vector_sha256"] == context["group_hash"], "actual scaffold hash mismatch")
        require(row["fold_sha256"] == context["fold_hash"], "actual fold hash mismatch")
        require(row["R_sha256"] == vector_hash(mask, "binary"), "actual mask hash mismatch")
        close(row["full_truth"], loss.mean(), tolerance, "actual truth")
        close(row["cc"], loss[mask.astype(bool)].mean(), tolerance, "actual complete case")
        close(row["cc_abs_error"], abs(row["cc"] - row["full_truth"]), tolerance, "actual CC error")
        close(row["fitted_aipw_abs_error"], abs(row["fitted_aipw"] - row["full_truth"]), tolerance, "actual fitted error")
        close(row["fitted_lcl"], row["fitted_aipw"] - 1.96 * row["fitted_se"], tolerance, "actual lower interval")
        close(row["fitted_ucl"], row["fitted_aipw"] + 1.96 * row["fitted_se"], tolerance, "actual upper interval")
        require(bool(row["fitted_covered"]) == bool(row["fitted_lcl"] <= row["full_truth"] <= row["fitted_ucl"]), "actual coverage mismatch")
        stored_smd = json.loads(row["fitted_smd_json"])
        require(set(stored_smd) == set(LOWDIM_COLUMNS), "actual SMD feature grid mismatch")
        require(np.isfinite(np.asarray(list(stored_smd.values()), dtype=float)).all(), "nonfinite actual SMD")
        close(row["fitted_max_abs_smd"], max(stored_smd.values()), tolerance, "actual max SMD")
        for column, value in support.items():
            require(int(row[column]) == value, f"actual {column} mismatch")
        expected_support = min(support.values()) >= 20
        require(bool(row["fold_class_support"]) == expected_support, "actual fold-support flag mismatch")
        require(bool(row["all_finite"]), "actual all_finite flag failed")
        eligible = bool(
            expected_support
            and float(row["fitted_clip_fraction"]) <= overlap_cfg["maximum_clip_fraction"]
            and float(row["fitted_ess_fraction"]) >= overlap_cfg["minimum_ess_over_observed"]
            and float(row["fitted_max_abs_smd"]) <= overlap_cfg["maximum_max_abs_smd"]
            and float(row["fitted_max_normalized_weight"]) <= overlap_cfg["maximum_normalized_observed_weight"]
        )
        require(bool(row["overlap_eligible"]) == eligible, "actual overlap flag mismatch")
        require(row["overlap_label"] == ("ELIGIBLE_OVERLAP" if eligible else overlap_cfg["failure_label"]), "actual overlap label mismatch")
    require(
        frame.groupby(["endpoint", "anchor", "propensity_representation"])[
            ["p_raw_sha256", "p_fitted_sha256", "fitted_ess", "fitted_ess_fraction", "fitted_clip_fraction"]
        ].nunique().max().max() == 1,
        "actual propensity output varies across model or loss representation",
    )
    require(
        frame.groupby(["endpoint", "anchor", "model", "loss_representation"])["mhat_sha256"].nunique().max() == 1,
        "actual loss nuisance varies across propensity representation",
    )


def refit_seed_nuisance(context, mask, loss_matrix, cfg):
    propensity = {}
    outcome = {}
    for representation, features in context["representations"].items():
        require(features.dtype == np.float64, f"non-float64 nuisance representation: {representation}")
        p_raw = np.empty(len(mask), dtype=float)
        mhat = np.empty_like(loss_matrix, dtype=float)
        for fold in range(3):
            training = np.flatnonzero(context["fold"] != fold)
            held_out = np.flatnonzero(context["fold"] == fold)
            scaler = StandardScaler().fit(features[training])
            x_training = scaler.transform(features[training])
            x_held_out = scaler.transform(features[held_out])
            require(len(np.unique(mask[training])) == 2, "seed-41001 propensity fold has one class")
            p_cfg = cfg["nuisance_models"]["propensity"]
            fitted_propensity = checked_fit(
                LogisticRegression(
                    penalty=p_cfg["penalty"],
                    C=float(p_cfg["C"]),
                    solver=p_cfg["solver"],
                    class_weight=p_cfg["class_weight"],
                    max_iter=int(p_cfg["max_iter"]),
                    random_state=int(p_cfg["random_state"]),
                ),
                x_training,
                mask[training],
                f"seed-41001 propensity {representation} fold {fold}",
            )
            p_raw[held_out] = fitted_propensity.predict_proba(x_held_out)[:, 1]
            observed = mask[training].astype(bool)
            require(observed.any(), "seed-41001 loss-nuisance fold has no observed rows")
            m_cfg = cfg["nuisance_models"]["loss"]
            fitted_loss = checked_fit(
                Ridge(
                    alpha=float(m_cfg["alpha"]),
                    solver=m_cfg["solver"],
                    tol=float(m_cfg["tol"]),
                ),
                x_training[observed],
                loss_matrix[training][observed],
                f"seed-41001 loss nuisance {representation} fold {fold}",
            )
            mhat[held_out] = fitted_loss.predict(x_held_out)
        require(np.isfinite(p_raw).all() and np.isfinite(mhat).all(), "nonfinite independent seed-41001 nuisance refit")
        propensity[representation] = {
            "p_raw": p_raw,
            "p_fitted": np.clip(p_raw, 0.05, 0.95),
        }
        outcome[representation] = mhat
    return {"propensity": propensity, "outcome": outcome}


def verify_seed_replay(replay, scenarios, contexts, contracts, losses, cfg, tolerance):
    require_schema(replay, "row_replay_seed41001.parquet")
    require(len(replay) == 59232, "seed-41001 replay row count mismatch")
    keys = ("endpoint", "anchor", "mechanism", "seed", "model", "propensity_representation", "loss_representation")
    expected_cells = set(product(ENDPOINTS, REGULAR_ANCHORS, MECHANISMS, (41001,), MODELS, P_REPRESENTATIONS, M_REPRESENTATIONS))
    require(set(replay[list(keys)].drop_duplicates().itertuples(index=False, name=None)) == expected_cells, "seed-41001 replay cell grid mismatch")
    require(not replay.duplicated([*keys, "source_row_index"]).any(), "duplicate seed replay row")
    finite_columns(
        replay,
        [
            "source_row_index", "cohort_row_index", "fold", "R", "loss", "mhat", "p_raw", "p_fitted", "p_oracle",
            *tuple(f"audit_{name}" for name in PXR_STRUCTURE_COLUMNS), "audit_max_train_similarity",
            "augmentation_fitted", "pseudo_outcome_fitted", "influence_fitted",
            "augmentation_oracle", "pseudo_outcome_oracle", "influence_oracle",
        ],
        "seed replay",
    )
    scenario_index = scenarios[scenarios["seed"].eq(41001)].set_index(list(keys))
    maximum_difference = {}
    nuisance_difference = {"p_raw": 0.0, "p_fitted": 0.0, "mhat": 0.0}

    refitted = {}
    for endpoint, anchor, mechanism in product(ENDPOINTS, REGULAR_ANCHORS, MECHANISMS):
        context = contexts[endpoint]
        mask, _ = regenerated_mask(context, contracts[(endpoint, anchor, mechanism)], 41001)
        loss_matrix = np.column_stack([losses[(endpoint, model)] for model in MODELS])
        refitted[(endpoint, anchor, mechanism)] = refit_seed_nuisance(
            context, mask, loss_matrix, cfg
        )

    def track(label, left, right):
        difference = np.abs(np.asarray(left, dtype=float) - np.asarray(right, dtype=float))
        maximum = float(difference.max()) if difference.size else 0.0
        maximum_difference[label] = max(maximum_difference.get(label, 0.0), maximum)
        require(np.isfinite(difference).all() and maximum <= tolerance, f"replay formula mismatch: {label}")

    def track_nuisance(label, left, right):
        difference = np.abs(np.asarray(left, dtype=float) - np.asarray(right, dtype=float))
        maximum = float(difference.max()) if difference.size else 0.0
        nuisance_difference[label] = max(nuisance_difference[label], maximum)
        require(
            np.isfinite(difference).all() and maximum <= tolerance,
            f"independent seed-41001 nuisance mismatch: {label}",
        )

    for cell, part in replay.groupby(list(keys), sort=False):
        endpoint, anchor, mechanism, seed, model, p_rep, m_rep = cell
        context = contexts[endpoint]
        contract = contracts[(endpoint, anchor, mechanism)]
        part = part.sort_values("source_row_index")
        require(len(part) == len(context["source_row_index"]), f"replay cell length mismatch: {cell}")
        same_vector(part["source_row_index"], context["source_row_index"], tolerance, "replay source rows")
        same_vector(part["cohort_row_index"], context["cohort_row_index"], tolerance, "replay cohort rows")
        same_vector(part["molecule_id"].astype(str), context["name"], tolerance, "replay molecule IDs")
        same_vector(part["scaffold_group"].astype(str), context["groups"], tolerance, "replay scaffold groups")
        same_vector(part["scaffold_sha256"].astype(str), [row_text_hash(value) for value in context["groups"]], tolerance, "replay row scaffold hashes")
        same_vector(part["fold"], context["fold"], tolerance, "replay folds")
        mask, mask_hash = regenerated_mask(context, contract, seed)
        same_vector(part["R"], mask, tolerance, "replay mask")
        same_vector(part["loss"], losses[(endpoint, model)], tolerance, "replay loss")
        same_vector(part["p_oracle"], contract["probability"], tolerance, "replay oracle propensity")
        track("p_fitted_clip", part["p_fitted"], np.clip(part["p_raw"].to_numpy(float), 0.05, 0.95))
        independent = refitted[(endpoint, anchor, mechanism)]
        track_nuisance("p_raw", part["p_raw"], independent["propensity"][p_rep]["p_raw"])
        track_nuisance("p_fitted", part["p_fitted"], independent["propensity"][p_rep]["p_fitted"])
        track_nuisance(
            "mhat",
            part["mhat"],
            independent["outcome"][m_rep][:, MODELS.index(model)],
        )
        replay_features = np.column_stack(
            [part[f"audit_{name}"].to_numpy(float) for name in PXR_STRUCTURE_COLUMNS]
            + [part["audit_max_train_similarity"].to_numpy(float)]
        )
        same_vector(replay_features, context["lowdim"], tolerance, "replay audit features")

        stored = scenario_index.loc[cell]
        require(stored["R_sha256"] == mask_hash, "replay-to-scenario mask hash mismatch")
        require(stored["p_oracle_sha256"] == vector_hash(part["p_oracle"], "float"), "replay oracle-p hash mismatch")
        require(stored["p_raw_sha256"] == vector_hash(part["p_raw"], "float"), "replay raw-p hash mismatch")
        require(stored["p_fitted_sha256"] == vector_hash(part["p_fitted"], "float"), "replay fitted-p hash mismatch")
        require(stored["mhat_sha256"] == vector_hash(part["mhat"], "float"), "replay mhat hash mismatch")
        fitted = estimate_from_rows(part["loss"], part["R"], part["p_fitted"], part["mhat"], context["groups"])
        oracle = estimate_from_rows(part["loss"], part["R"], part["p_oracle"], part["mhat"], context["groups"])
        track("truth", stored["full_truth"], part["loss"].mean())
        track("complete_case", stored["cc"], part.loc[part["R"].eq(1), "loss"].mean())
        track("outcome_regression", stored["or_mean"], part["mhat"].mean())
        track(
            "clip_fraction",
            stored["fitted_clip_fraction"],
            np.mean((part["p_raw"] < 0.05) | (part["p_raw"] > 0.95)),
        )
        for prefix, values in (("fitted", fitted), ("oracle", oracle)):
            for stored_name, computed_name in (
                ("ipw", "ipw"), ("aipw", "aipw"), ("se", "aipw_se"), ("lcl", "aipw_lcl"),
                ("ucl", "aipw_ucl"), ("ess", "ess"), ("ess_fraction", "ess_over_observed"),
                ("max_normalized_weight", "maximum_normalized_observed_weight"),
            ):
                track(f"{prefix}_{stored_name}", stored[f"{prefix}_{stored_name}"], values[computed_name])
            for name in ("augmentation", "pseudo_outcome", "influence"):
                expected_name = "pseudo_outcome" if name == "pseudo_outcome" else name
                track(f"{prefix}_{expected_name}", part[f"{name}_{prefix}"], values[name])
        fitted_smd = smd_map(part["R"], part["p_fitted"], replay_features)
        oracle_smd = smd_map(part["R"], part["p_oracle"], replay_features)
        compare_tree(fitted_smd, json.loads(stored["fitted_smd_json"]), tolerance, "replay.fitted_smd")
        compare_tree(oracle_smd, json.loads(stored["oracle_smd_json"]), tolerance, "replay.oracle_smd")
        for prefix, values in (("fitted", fitted_smd), ("oracle", oracle_smd)):
            track(f"{prefix}_max_abs_smd", stored[f"{prefix}_max_abs_smd"], max(values.values()))
            for name in LOWDIM_COLUMNS:
                track(f"{prefix}_smd_{name}", json.loads(stored[f"{prefix}_smd_json"])[name], values[name])
    require(
        replay.groupby(["endpoint", "anchor", "mechanism", "seed", "propensity_representation", "source_row_index"])[
            ["p_raw", "p_fitted"]
        ].nunique().max().max() == 1,
        "replay propensity varies across model or loss representation",
    )
    require(
        replay.groupby(["endpoint", "anchor", "mechanism", "seed", "model", "loss_representation", "source_row_index"])["mhat"].nunique().max() == 1,
        "replay mhat varies across propensity representation",
    )
    require(
        replay.groupby(["endpoint", "anchor", "mechanism", "seed", "source_row_index"])[["R", "p_oracle"]].nunique().max().max() == 1,
        "replay mask/oracle probability not shared across factorial cells",
    )
    return maximum_difference, nuisance_difference


def verify_common_manifest(manifest, mode, status, inputs, runner_sha, verifier_sha, environment):
    require(manifest["study_id"] == "confirm_001", "manifest study identity mismatch")
    require(manifest["mode"] == mode and manifest["status"] == status, "manifest mode/status mismatch")
    require(manifest["input_hashes"] == inputs, "manifest input-hash map mismatch")
    require(manifest["runner_sha256"] == runner_sha, "manifest runner hash mismatch")
    require(manifest["verifier_sha256"] == verifier_sha, "manifest verifier hash mismatch")
    require(manifest["environment"] == environment, "manifest environment mismatch")
    require(isinstance(manifest["command"], list) and manifest["command"], "manifest command missing")
    require(f"--mode" in manifest["command"] and mode in manifest["command"], "manifest command mode mismatch")
    require(float(manifest["duration_seconds"]) >= 0 and np.isfinite(float(manifest["duration_seconds"])), "manifest duration invalid")
    started = datetime.fromisoformat(str(manifest["started_utc"]).replace("Z", "+00:00"))
    finished = datetime.fromisoformat(str(manifest["finished_utc"]).replace("Z", "+00:00"))
    require(finished >= started, "manifest timestamps reversed")


def require_exact_output_files(output, expected_without_manifest):
    actual = {
        str(path.relative_to(output)).replace("\\", "/")
        for path in output.rglob("*")
        if path.is_file()
    }
    require(actual == set(expected_without_manifest) | {"run_manifest.json"}, "missing, extra, or nested output file")


def verify_output_hashes(output, stored, expected_names):
    require(set(stored) == set(expected_names), "manifest output-hash grid mismatch")
    for relative, wanted in stored.items():
        require(re.fullmatch(r"[0-9a-f]{64}", wanted or "") is not None, f"invalid output hash: {relative}")
        require(sha256(output / relative) == wanted, f"output hash mismatch: {relative}")


def expected_integrity(maximum_difference, tolerance):
    return {
        "pass": True,
        "status": "PASS",
        "checks": {
            "input_hashes_and_environment": True,
            "cohort_and_structure": True,
            "exact_cartesian_grids": True,
            "shared_prediction_mask_and_fold_vectors": True,
            "row_replay_seed41001": True,
            "formula_replay_within_tolerance": True,
            "gates_recomputed_from_reloaded_outputs": True,
            "formal_directory_absent_before_launch": True,
            "formal_directory_written_once": True,
        },
        "replay": {
            "pass": True,
            "rows": 59232,
            "cells": 64,
            "tolerance": tolerance,
            "maximum_absolute_difference": maximum_difference,
        },
    }


def expected_decision(hypotheses, admission, boundary):
    if not hypotheses["H1"]["pass"] or not hypotheses["H2"]["pass"]:
        return "STOP_NO_EXTERNAL_CONFIRMATION"
    boundary_complete = bool(
        boundary["decision"] == "ABSTAIN_NONIDENTIFIABLE_FULL_LIBRARY_DOSE_RESPONSE"
        and boundary["full_library_point_estimate_allowed"] is False
    )
    require(boundary_complete, "PXR boundary incomplete")
    if (
        admission["known_truth_admission"]["pass"]
        and admission["direct_mask_overlap"]["pass"]
        and admission["direct_mask_translation"]["pass"]
    ):
        return "GO_WRITE_CORRECTION_BOUNDED"
    return "GO_WRITE_FALSIFICATION_BOUNDARY"


def expected_verification_text(integrity, hypotheses, decision):
    return (
        "# confirm_001 runner verification\n\n"
        f"- Internal durable-grid and seed-41001 formula replay: `{'PASS' if integrity['pass'] else 'FAIL'}`.\n"
        f"- H1 all four cells: `{'PASS' if hypotheses['H1']['pass'] else 'FAIL'}`.\n"
        f"- H2 all four cells: `{'PASS' if hypotheses['H2']['pass'] else 'FAIL'}`.\n"
        f"- Frozen decision: `{decision}`.\n\n"
        "Scope: this file records the runner's reload-and-replay check. It is not the independent full verifier, "
        "does not alter `STOP_PILOT_NEGATIVE`, and does not identify full-library PXR dose-response performance.\n"
    )


def verify_full_manifest(
    output,
    manifest,
    inputs,
    runner_sha,
    verifier_sha,
    environment,
    config_overlay,
    recovery,
    cfg,
    contexts,
    contracts,
    structure,
    pxr_context,
    integrity,
    decision,
    decision_006_sha,
    decision_007_sha,
):
    required = {
        "study_id", "mode", "status", "started_utc", "finished_utc", "duration_seconds", "command",
        "input_precedence", "input_hashes", "runner_sha256", "verifier_sha256", "clearance_hashes",
        "environment", "config_paths", "seeds", "feature_contract", "grid_rows", "vector_hashes",
        "pxr_structure_audit", "integrity", "decision", "output_hashes_before_manifest",
    }
    require(set(manifest) == required, "full run-manifest top-level keys mismatch")
    verify_common_manifest(
        manifest, "recovery-full", "COMPLETED_PENDING_INDEPENDENT_FULL_VERIFICATION",
        inputs, runner_sha, verifier_sha, environment,
    )
    require(
        manifest["clearance_hashes"]
        == {DECISION_006_PATH: decision_006_sha, DECISION_007_PATH: decision_007_sha},
        "full clearance hashes mismatch",
    )
    require(manifest["input_precedence"] == recovery["formal_input_precedence"], "manifest input precedence mismatch")
    require(
        manifest["config_paths"]
        == {
            "base_protocol": "artifacts/confirmatory/confirm_001/protocol.md",
            "base_config": "artifacts/confirmatory/confirm_001/config.json",
            "source_manifest": "data/external/source_freeze_20260911/source_manifest.json",
            "source_manifest_amendment": "data/external/source_freeze_20260911/source_manifest_amendment_001.json",
            "protocol_amendment_001": "artifacts/confirmatory/confirm_001/protocol_amendment_001.md",
            "config_amendment_001": "artifacts/confirmatory/confirm_001/config_amendment_001.json",
            "failure_record": "artifacts/confirmatory/confirm_001/failure_001.md",
            "protocol_amendment_002": "artifacts/confirmatory/confirm_001/protocol_amendment_002.md",
            "config_amendment_002": "artifacts/confirmatory/confirm_001/config_amendment_002.json",
        },
        "manifest config paths mismatch",
    )
    bootstrap_hash = cfg["confirmatory_hypotheses"]["bootstrap"]["seed_index_matrix_sha256_numpy_1_26_4"]
    require(
        manifest["seeds"]
        == {
            "mask_seeds": {"first": 41001, "last": 41100, "count": 100},
            "row_replay_seed": 41001,
            "bootstrap_seed": 52001,
            "bootstrap_resamples": 10000,
            "bootstrap_seed_index_matrix_sha256": bootstrap_hash,
        },
        "manifest seed contract mismatch",
    )
    expected_features = {
        "descriptor_registry_count": 217,
        "descriptor_names": structure["descriptor_names"],
        "descriptor_names_sha256": vector_hash(structure["descriptor_names"], "text"),
        "nuisance_representations": {
            endpoint: {
                "lowdim_11": {"shape": [len(contexts[endpoint]["source_row_index"]), 11], "dtype": "float64"},
                "full_2266": {"shape": [len(contexts[endpoint]["source_row_index"]), 2266], "dtype": "float64"},
            }
            for endpoint in ENDPOINTS
        },
        "lowdim_names": list(LOWDIM_COLUMNS),
        "mcar_generator_neutral_storage": {"generator_score": 0.0, "intercept_shift": 0.0},
    }
    compare_tree(expected_features, manifest["feature_contract"], 0, "manifest.feature_contract")
    expected_grid = {"cohort_audit.csv": 18, **EXPECTED_ROWS}
    require(manifest["grid_rows"] == expected_grid, "manifest row-grid map mismatch")
    biogen_hashes = {
        endpoint: {
            "source_row_index_sha256": contexts[endpoint]["source_hash"],
            "scaffold_group_vector_sha256": contexts[endpoint]["group_hash"],
            "fold_sha256": contexts[endpoint]["fold_hash"],
            "generator_p_oracle_sha256": {
                f"{anchor}_{mechanism}": contracts[(endpoint, anchor, mechanism)]["probability_hash"]
                for anchor in REGULAR_ANCHORS
                for mechanism in MECHANISMS
            },
        }
        for endpoint in ENDPOINTS
    }
    pxr_hashes = {
        "ordered_source_row_index": pxr_context["hashes"]["ordered_source_row_index"],
        "ordered_murcko_group": pxr_context["hashes"]["ordered_murcko_group"],
        "ordered_fold": pxr_context["hashes"]["ordered_fold"],
        "ordered_canonical_smiles": pxr_context["hashes"]["ordered_canonical_smiles"],
        "ordered_parent_OCNT_ID": pxr_context["hashes"]["ordered_OCNT_ID"],
        "ordered_availability": pxr_context["hashes"]["ordered_availability"],
        "ordered_formal_hit": pxr_context["hashes"]["ordered_formal_hit"],
    }
    require(
        manifest["vector_hashes"]
        == {"biogen": biogen_hashes, "pxr": pxr_hashes, "bootstrap_seed_index_matrix_sha256": bootstrap_hash},
        "manifest vector-hash contract mismatch",
    )
    expected_pxr_structure = {
        "rows": 10747,
        "smiles_parse_failures": 0,
        "unique_canonical_smiles": 10747,
        "descriptor_10_all_finite_rows": 10747,
        "ordered_canonical_smiles_sha256": pxr_hashes["ordered_canonical_smiles"],
        "ordered_murcko_group_sha256": pxr_hashes["ordered_murcko_group"],
    }
    require(manifest["pxr_structure_audit"] == expected_pxr_structure, "manifest PXR structure audit mismatch")
    compare_tree(integrity, manifest["integrity"], cfg["integrity"]["replay_tolerance"], "manifest.integrity")
    require(manifest["decision"] == decision, "manifest scientific decision mismatch")
    require_exact_output_files(output, EXPECTED_OUTPUTS)
    verify_output_hashes(output, manifest["output_hashes_before_manifest"], EXPECTED_OUTPUTS)


def smoke_crossfit(features, mask, losses, folds):
    prepared = {}
    for name, matrix in features.items():
        prepared[name] = []
        for fold in range(3):
            training = np.flatnonzero(folds != fold)
            held_out = np.flatnonzero(folds == fold)
            scaler = StandardScaler().fit(matrix[training])
            prepared[name].append((training, held_out, scaler.transform(matrix[training]), scaler.transform(matrix[held_out])))
    propensity = {}
    outcome = {}
    for name in P_REPRESENTATIONS:
        raw = np.empty(len(mask), dtype=float)
        mhat = np.empty_like(losses, dtype=float)
        for training, held_out, x_training, x_held_out in prepared[name]:
            require(len(np.unique(mask[training])) == 2, "smoke propensity fold has one class")
            fitted_propensity = checked_fit(LogisticRegression(
                penalty="l2", C=1.0, solver="liblinear", class_weight=None,
                max_iter=1000, random_state=20260910,
            ), x_training, mask[training], "smoke propensity")
            raw[held_out] = fitted_propensity.predict_proba(x_held_out)[:, 1]
            observed = mask[training].astype(bool)
            fitted_loss = checked_fit(
                Ridge(alpha=10.0, solver="lsqr", tol=1e-6),
                x_training[observed],
                losses[training][observed],
                "smoke loss nuisance",
            )
            mhat[held_out] = fitted_loss.predict(x_held_out)
        propensity[name] = {"raw": raw, "clipped": np.clip(raw, 0.05, 0.95)}
        outcome[name] = mhat
    return propensity, outcome


def verify_smoke(output, manifest, inputs, runner_sha, verifier_sha, environment, decision_006_sha, recovery):
    required_manifest = {
        "study_id", "mode", "status", "started_utc", "finished_utc", "duration_seconds", "command",
        "input_hashes", "runner_sha256", "verifier_sha256", "clearance_hashes", "environment",
        "synthetic_only", "output_hashes_before_manifest",
    }
    require(set(manifest) == required_manifest, "smoke manifest top-level keys mismatch")
    verify_common_manifest(
        manifest, "recovery-smoke", "PASS_NON_ADJUDICATIVE", inputs, runner_sha, verifier_sha, environment
    )
    require(manifest["clearance_hashes"] == {DECISION_006_PATH: decision_006_sha}, "smoke clearance hash mismatch")
    require(manifest["synthetic_only"] is True, "smoke is not marked synthetic-only")
    require_exact_output_files(output, {"roundtrip_fixture.csv", "smoke.json"})
    verify_output_hashes(output, manifest["output_hashes_before_manifest"], {"roundtrip_fixture.csv", "smoke.json"})

    smoke = json.loads((output / "smoke.json").read_text(encoding="utf-8"))
    require(
        set(smoke)
        == {
            "study_id", "mode", "status", "non_adjudicative", "synthetic_contract", "checks",
            "fold_sha256", "input_hashes", "runner_sha256", "verifier_sha256",
        },
        "smoke JSON top-level keys mismatch",
    )
    rng = np.random.default_rng(20260911)
    low = rng.normal(size=(72, 11)).astype(np.float64)
    full = np.column_stack((low, rng.normal(size=(72, 2255)))).astype(np.float64)
    groups = np.repeat([f"SMOKE_GROUP_{index}" for index in range(24)], 3)
    folds = fold_vector(groups)
    mask = np.tile([0, 1], 36).astype(np.uint8)
    losses = np.abs(rng.normal(size=(72, 2)))
    propensity, outcome = smoke_crossfit(
        {"lowdim_11": low, "full_2266": full}, mask, losses, folds
    )
    estimates = []
    for p_rep, m_rep, model_index in product(P_REPRESENTATIONS, M_REPRESENTATIONS, range(2)):
        result = estimate_from_rows(
            losses[:, model_index], mask, propensity[p_rep]["clipped"], outcome[m_rep][:, model_index], groups
        )
        estimates.append([result[name] for name in ("ipw", "aipw", "aipw_se", "ess", "ess_over_observed")])
    fixture = read_table(output / "roundtrip_fixture.csv")
    require(tuple(fixture.columns) == ("fixture_id", "numerator", "denominator", "ratio") and len(fixture) == 1, "roundtrip fixture schema mismatch")
    numerator = 408.40724836686263
    denominator = 0.004053686965995372
    expected_fixture = pd.DataFrame({
        "fixture_id": ["synthetic_large_ratio"],
        "numerator": [numerator],
        "denominator": [denominator],
        "ratio": [numerator / denominator],
    })
    require(np.array_equal(fixture["fixture_id"].astype(str), expected_fixture["fixture_id"].astype(str)), "roundtrip fixture identity mismatch")
    numeric_columns = ["numerator", "denominator", "ratio"]
    component_bit_exact = bool(np.array_equal(
        fixture[numeric_columns].to_numpy(dtype="<f8").view("<u8"),
        expected_fixture[numeric_columns].to_numpy(dtype="<f8").view("<u8"),
    ))
    stored_ratio = float(fixture.loc[0, "ratio"])
    recomputed_ratio = float(fixture.loc[0, "numerator"] / fixture.loc[0, "denominator"])
    ratio_difference = abs(stored_ratio - recomputed_ratio)
    support = fold_support(mask, folds)
    checks = {
        "folds_complete": bool(np.array_equal(np.unique(folds), np.arange(3))),
        "fold_class_support": bool(min(support.values()) >= 20),
        "raw_propensity_finite": bool(all(np.isfinite(value["raw"]).all() for value in propensity.values())),
        "clipped_propensity_finite": bool(all(np.isfinite(value["clipped"]).all() for value in propensity.values())),
        "loss_nuisance_finite": bool(all(np.isfinite(value).all() for value in outcome.values())),
        "estimates_finite": bool(np.isfinite(np.asarray(estimates, dtype=float)).all()),
        "nuisance_dtype_float64": bool(low.dtype == full.dtype == np.float64),
        "csv_roundtrip_component_bit_exact": component_bit_exact,
        "csv_roundtrip_ratio_within_absolute_tolerance": bool(ratio_difference <= 1e-12),
        "csv_roundtrip_large_ratio": bool(abs(stored_ratio) >= recovery["csv_roundtrip_contract"]["recovery_smoke_minimum_ratio_magnitude"]),
    }
    expected = {
        "study_id": "confirm_001",
        "mode": "recovery-smoke",
        "status": "PASS",
        "non_adjudicative": True,
        "synthetic_contract": {
            "random_state": 20260911,
            "rows": 72,
            "groups": 24,
            "rows_per_group": 3,
            "availability_pattern": [0, 1],
            "propensity_representations": {"lowdim_11": [72, 11], "full_2266": [72, 2266]},
            "loss_columns": 2,
            "csv_roundtrip_fixture": {
                "synthetic_only": True,
                "float_format": "%.17g",
                "float_precision": "round_trip",
                "rows": 1,
                "ratio_magnitude": abs(stored_ratio),
                "derived_ratio_absolute_difference": ratio_difference,
            },
        },
        "checks": checks,
        "fold_sha256": vector_hash(folds, "int"),
        "input_hashes": inputs,
        "runner_sha256": runner_sha,
        "verifier_sha256": verifier_sha,
    }
    require(all(checks.values()), "independent smoke reconstruction failed")
    compare_tree(expected, smoke, 0, "smoke")


def verify_full(
    root, output, cfg, corrections, config_overlay, recovery, inputs, runner_sha, verifier_sha,
    environment, decision_006_sha, decision_007_sha,
):
    require_exact_output_files(output, EXPECTED_OUTPUTS)
    manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    verify_output_hashes(output, manifest.get("output_hashes_before_manifest", {}), EXPECTED_OUTPUTS)
    tolerance = float(cfg["integrity"]["replay_tolerance"])
    data, contexts, structure = load_biogen_contexts(root, cfg, corrections)
    contracts = generator_contracts(contexts, cfg)
    pxr_context = load_pxr_context(root, cfg, corrections)
    cohort = read_table(output / "cohort_audit.csv")
    availability = read_table(output / "availability_map.csv")
    predictions = read_table(output / "predictions.parquet")
    folds = read_table(output / "nuisance_fold_assignments.parquet")
    generators = read_table(output / "generator_vectors.parquet")
    scenarios = read_table(output / "semisynthetic_scenarios.csv")
    actual = read_table(output / "actual_mask_scenarios.csv")
    replay = read_table(output / "row_replay_seed41001.parquet")
    contrasts = read_table(output / "factorial_contrasts.csv")
    bootstrap = read_table(output / "bootstrap_intervals.csv")
    pxr_linkage = read_table(output / "pxr_linkage_audit.csv")
    pxr_strata = read_table(output / "pxr_selection_strata.csv")
    pxr_overlap = read_table(output / "pxr_overlap_summary.csv")

    verify_cohort_audit(cohort, data, tolerance)
    stored_structure = json.loads((output / "structure_audit.json").read_text(encoding="utf-8"))
    compare_tree(expected_structure_audit(structure), stored_structure, tolerance, "structure_audit")
    verify_availability_map(availability, contexts, tolerance)
    losses = verify_predictions(predictions, data, contexts, structure, cfg, tolerance)
    verify_folds(folds, contexts, tolerance)
    verify_generators(generators, contexts, contracts, tolerance)
    hypotheses = verify_scenario_grid(scenarios, contexts, contracts, losses, cfg, tolerance)
    verify_actual_grid(actual, contexts, losses, cfg, tolerance)
    maximum_difference, nuisance_difference = verify_seed_replay(
        replay, scenarios, contexts, contracts, losses, cfg, tolerance
    )
    compare_frames(
        factorial_records(scenarios), contrasts, "factorial_contrasts.csv",
        ("row_type", "endpoint", "anchor", "mechanism", "model", "loss_representation", "contrast"), tolerance,
    )
    compare_frames(
        bootstrap_records(scenarios, cfg), bootstrap, "bootstrap_intervals.csv",
        ("hypothesis", "statistic", "endpoint", "mechanism"), tolerance,
    )
    admission = admission_result(scenarios, cfg)
    admission.update(direct_result(actual, admission, cfg))
    stored_admission = json.loads((output / "admission_gate.json").read_text(encoding="utf-8"))
    compare_tree(admission, stored_admission, tolerance, "admission_gate")
    _, expected_overlap = verify_pxr_linkage(pxr_linkage, pxr_context, cfg, tolerance)
    compare_frames(expected_overlap, pxr_overlap, "pxr_overlap_summary.csv", ("specification",), tolerance)
    expected_strata, boundary = pxr_strata_and_boundary(pxr_context, expected_overlap, cfg)
    compare_frames(expected_strata, pxr_strata, "pxr_selection_strata.csv", ("formal_hit",), tolerance)
    stored_boundary = json.loads((output / "pxr_boundary_decision.json").read_text(encoding="utf-8"))
    compare_tree(boundary, stored_boundary, tolerance, "pxr_boundary_decision")

    integrity = expected_integrity(maximum_difference, tolerance)
    decision = expected_decision(hypotheses, admission, boundary)
    summary = {
        "study_id": "confirm_001",
        "pilot_decision_unchanged": "STOP_PILOT_NEGATIVE",
        "forensic_role": "development_only_no_confirmatory_aggregation",
        "integrity": integrity,
        "H1": hypotheses["H1"],
        "H2": hypotheses["H2"],
        "conjunctive_confirmation": bool(hypotheses["H1"]["pass"] and hypotheses["H2"]["pass"]),
        "admission": admission["known_truth_admission"],
        "direct_overlap": admission["direct_mask_overlap"],
        "direct_translation": admission["direct_mask_translation"],
        "pxr_boundary": boundary,
        "decision": decision,
        "anti_claims": cfg["anti_claims"],
    }
    stored_summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    compare_tree(summary, stored_summary, tolerance, "summary")
    require(
        (output / "verification.md").read_text(encoding="utf-8")
        == expected_verification_text(integrity, hypotheses, decision),
        "runner verification narrative mismatch",
    )
    verify_full_manifest(
        output, manifest, inputs, runner_sha, verifier_sha, environment, config_overlay, recovery, cfg,
        contexts, contracts, structure, pxr_context, integrity, decision,
        decision_006_sha, decision_007_sha,
    )
    return decision, nuisance_difference


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("recovery-smoke", "recovery-full"), required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--decision-004-sha")
    parser.add_argument("--decision-005-sha")
    parser.add_argument("--decision-006-sha", required=True)
    parser.add_argument("--decision-007-sha")
    parser.add_argument("--verification-output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    config_path = root / "artifacts/confirmatory/confirm_001/config.json"
    protocol_path = root / "artifacts/confirmatory/confirm_001/protocol.md"
    require(sha256(protocol_path) == EXPECTED_PROTOCOL_SHA256, "confirmatory protocol hash mismatch")
    require(sha256(config_path) == EXPECTED_CONFIG_SHA256, "confirmatory config hash mismatch")
    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    require(
        cfg["study_id"] == "confirm_001"
        and cfg["status"] == "FROZEN_PRE_RESULT"
        and cfg["formal_run_authorized"] is True
        and cfg["protocol_sha256"] == EXPECTED_PROTOCOL_SHA256,
        "base confirmatory authorization contract mismatch",
    )
    corrections, config_overlay = amended_hashes(root)
    recovery = load_recovery_amendment(root)
    environment = check_environment(cfg)
    inputs, runner_sha, verifier_sha = expected_input_hashes(
        root, cfg, config_path, config_overlay, recovery, environment, args.mode,
        args.decision_004_sha, args.decision_005_sha,
        args.decision_006_sha, args.decision_007_sha,
    )
    recovery_full = (root / RECOVERY_FULL_OUTPUT_PATH).resolve()
    if args.mode == "recovery-full":
        require(args.verification_output is None, "--verification-output is reserved for recovery-smoke")
        require(args.decision_007_sha is not None, "recovery-full verification requires Decision 007")
        output = recovery_full if args.output is None else args.output.resolve()
        require(output == recovery_full, "recovery-full verification output path mismatch")
        require(output.is_dir(), "recovery-full output directory is missing")
        decision, nuisance_difference = verify_full(
            root, output, cfg, corrections, config_overlay, recovery, inputs, runner_sha, verifier_sha,
            environment, args.decision_006_sha, args.decision_007_sha,
        )
        result = {
            "status": "PASS",
            "mode": "recovery-full",
            "decision": decision,
            "seed41001_independent_nuisance_refit_maximum_absolute_difference": nuisance_difference,
        }
    else:
        require(args.decision_007_sha is None, "Decision 007 is not accepted for recovery-smoke verification")
        require(args.verification_output is not None, "recovery-smoke verification requires --verification-output")
        verification_output = args.verification_output.resolve()
        require(verification_output == (root / RECOVERY_SMOKE_VERIFICATION_PATH).resolve(), "recovery-smoke verification output path mismatch")
        require(not (root / DECISION_007_PATH).exists(), "Decision 007 must be absent before recovery smoke verification")
        require(not recovery_full.exists(), "recovery full output exists before recovery smoke verification")
        require(not verification_output.exists(), "recovery smoke verification record must be absent before verification")
        require(args.output is not None, "recovery-smoke verification requires --output")
        output = args.output.resolve()
        require(output == (root / RECOVERY_SMOKE_OUTPUT_PATH).resolve(), "recovery-smoke output path mismatch")
        require(output != recovery_full and recovery_full not in output.parents, "recovery-smoke output is within recovery-full directory")
        require(output.is_dir(), "recovery-smoke output directory is missing")
        manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
        verify_smoke(
            output, manifest, inputs, runner_sha, verifier_sha, environment,
            args.decision_006_sha, recovery,
        )
        result = {
            "status": "PASS",
            "mode": "recovery-smoke",
            "non_adjudicative": True,
            "runner_sha256": runner_sha,
            "verifier_sha256": verifier_sha,
            "decision_006_sha256": args.decision_006_sha,
            "smoke_output_path": RECOVERY_SMOKE_OUTPUT_PATH,
            "smoke_json_sha256": sha256(output / "smoke.json"),
            "smoke_run_manifest_sha256": sha256(output / "run_manifest.json"),
            "roundtrip_fixture_sha256": sha256(output / "roundtrip_fixture.csv"),
        }
    result["scope"] = RECOVERY_SMOKE_SCOPE if args.mode == "recovery-smoke" else {
        "independent_source_feature_prediction_generator_and_seed41001_formula_reconstruction": True,
        "independent_seed41001_crossfit_propensity_and_ridge_refit": True,
        "independent_remaining_99_seed_and_direct_or_pxr_nuisance_refit": False,
        "remaining_formal_nuisance_vectors_checked_by_hash_and_grid_invariants": True,
        "prelaunch_directory_absence_checked_from_runner_record_only": True,
    }
    if args.mode == "recovery-smoke":
        with verification_output.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(result, handle, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False)
            handle.write("\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(json.dumps({"status": "FAIL", "error": str(error)}, ensure_ascii=False), file=sys.stderr)
        raise
