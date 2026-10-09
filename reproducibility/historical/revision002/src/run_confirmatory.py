import argparse
import hashlib
import json
import platform
import re
import sys
import time
import warnings
from datetime import datetime, timezone
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow
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


BASE_PROTOCOL_SHA = "1819fa04713ee1a37df66c847825e5de3adcba736fba47c5c40b1c8db7796d7a"
BASE_CONFIG_SHA = "04fc27fc407a562853dd362d31a4c4cbc01bdb754199d5a18014cdfcbc7a5c23"
DECISION_003_SHA = "0b784a409faa5858fecdd70bf77e58fbb71ab85076bf343b89bc8c7a138df9ec"
PROTOCOL_AUDIT_SHA = "a4b4c630785c2cd6991872696d2bc2cf23a14235cc19c7331ee1f70731fb5b16"
SOURCE_MANIFEST_SHA = "be55555aed0b2a5bed095ae8d052500becd527c97c809d6db5b7a14e7548fd78"
SOURCE_AMENDMENT_SHA = "5b6b3337a9b7b1a4d359062ab07807527accc37b5cb63798db69e5f23837dfe5"
PROTOCOL_AMENDMENT_SHA = "989eac8ebd6e26bdcb5cee454ec51f4c00e76ad93bf4b245b92e0a59e9ff4a49"
CONFIG_AMENDMENT_SHA = "9008878d6a8c23aae059cf3b4aefcdc813bdd7387cbf1aa84fd6d38db620ec3a"
AMENDMENT_AUDIT_SHA = "cc08015a2a8ea126b90cfd0a53af1a90d548a09f890ce86371490c0bdbc158cf"
PROTOCOL_AMENDMENT_002_SHA = "d34c1053da408a3e349c2ecdc20d1e49f8aa5cda9dd491746146d621c5988d9f"
CONFIG_AMENDMENT_002_SHA = "144af9c376929f24699452051abed2ddb5243ab8cdf47c8f4024e219ac278984"
FAILURE_001_SHA = "a75a4b853bf9ce5469873870d0375e077f2dd2be7d5e6981febf5abd41667267"
AMENDMENT_002_AUDIT_SHA = "7a266a7614729ccf5bb2ff13890eabb6d2e6e32cc0ebc66ea21a34c8f1a821c3"
DESCRIPTOR_NAMES_SHA = "84cea0bf91c3c74a704770fd5d5477dc689fc3c09e0997cd268d070acc61be73"
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

ENDPOINTS = ["LOG_HLM", "LOG_RLM"]
REGULAR_ANCHORS = ["LOG_SOL", "LOG_MDR1_ER"]
STRESS_ANCHORS = ["LOG_HPPB", "LOG_RPPB"]
MODELS = ["ecfp_rf", "ecfp_desc_xgb"]
REPRESENTATIONS = ["lowdim_11", "full_2266"]
MECHANISMS = ["mcar", "mar"]
LOW_DESCRIPTOR_NAMES = [
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
]

COHORT_COLUMNS = ["dataset", "split", "endpoint", "n_rows", "n_finite", "n_unavailable", "coverage"]
AVAILABILITY_COLUMNS = [
    "endpoint", "anchor", "n_total", "n_available", "n_unavailable", "coverage", "availability_sha256"
]
PREDICTION_COLUMNS = [
    "endpoint", "model", "source_row_index", "cohort_row_index", "molecule_id", "reference", "prediction", "loss"
]
FOLD_COLUMNS = [
    "endpoint", "source_row_index", "cohort_row_index", "molecule_id", "scaffold_group",
    "scaffold_sha256", "fold", "source_row_index_sha256", "scaffold_group_vector_sha256", "fold_sha256",
]
GENERATOR_COLUMNS = [
    "endpoint", "anchor", "mechanism", "source_row_index", "cohort_row_index", "direct_coverage",
    "generator_score", "intercept_shift", "p_oracle", "p_oracle_sha256", "coverage_error", "all_finite",
]
FOLD_SUPPORT_COLUMNS = [
    f"fold_{fold}_train_{state}" for fold in range(3) for state in ["observed", "unavailable"]
]
SEMISYNTHETIC_COLUMNS = [
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
]
ACTUAL_COLUMNS = [
    "endpoint", "anchor", "mask_type", "model", "propensity_representation", "loss_representation",
    "n_total", "n_observed", "observed_fraction", "source_row_index_sha256", "scaffold_group_vector_sha256",
    "fold_sha256", "R_sha256", "p_raw_sha256", "p_fitted_sha256", "mhat_sha256", "full_truth", "cc",
    "or_mean", "fitted_ipw", "fitted_aipw", "fitted_se", "fitted_lcl", "fitted_ucl", "fitted_covered",
    "fitted_ess", "fitted_ess_fraction", "fitted_clip_fraction", "fitted_max_abs_smd",
    "fitted_max_normalized_weight", "fitted_smd_json", *FOLD_SUPPORT_COLUMNS, "fold_class_support", "all_finite",
    "overlap_eligible", "overlap_label", "cc_abs_error", "fitted_aipw_abs_error",
]
REPLAY_COLUMNS = [
    "endpoint", "anchor", "mechanism", "seed", "model", "propensity_representation", "loss_representation",
    "source_row_index", "cohort_row_index", "molecule_id", "scaffold_group", "scaffold_sha256", "fold", "R",
    "loss", "mhat", "p_raw", "p_fitted", "p_oracle", *[f"audit_{name}" for name in LOW_DESCRIPTOR_NAMES],
    "audit_max_train_similarity", "augmentation_fitted", "pseudo_outcome_fitted", "influence_fitted",
    "augmentation_oracle", "pseudo_outcome_oracle", "influence_oracle",
]
CONTRAST_COLUMNS = [
    "row_type", "endpoint", "anchor", "mechanism", "model", "loss_representation", "contrast", "estimate",
    "numerator", "denominator", "ratio", "rows", "status",
]
BOOTSTRAP_COLUMNS = [
    "hypothesis", "statistic", "endpoint", "mechanism", "point_estimate", "lower", "upper", "threshold",
    "direction", "point_pass", "denominator_valid", "resamples", "seed_index_matrix_sha256", "status",
]
PXR_LINKAGE_COLUMNS = [
    "source_row_index", "cohort_row_index", "OCNT_ID", "R_followup", "formal_hit", "smiles_parse_ok",
    "descriptor_10_all_finite", "canonical_smiles_sha256", "murcko_group_sha256", "fold", *LOW_DESCRIPTOR_NAMES,
    "log2_fc_estimate", "fdr_bh", "fdr_score", "p_raw_structure_10", "p_clipped_structure_10",
    "clipped_structure_10", "p_raw_structure_screen_12", "p_clipped_structure_screen_12",
    "clipped_structure_screen_12", "exact_smiles_match_if_linked", "batch_match_if_linked",
    "cross_file_exact_smiles_nunique", "cross_file_batch_nunique",
]
PXR_STRATA_COLUMNS = [
    "formal_hit", "n_total", "n_followup_available", "n_followup_unavailable", "advancement_probability"
]
PXR_OVERLAP_COLUMNS = [
    "specification", "n_total", "n_observed", "n_unavailable", "coverage", "clip_fraction", "ess",
    "ess_fraction", "max_abs_smd", "max_normalized_weight", "smd_json", "zero_variance_features_json",
    *FOLD_SUPPORT_COLUMNS, "fold_class_support", "all_finite", "source_row_index_sha256",
    "murcko_group_vector_sha256", "fold_sha256", "p_raw_sha256", "p_clipped_sha256",
    "availability_sha256", "formal_hit_sha256", "trigger_clip", "trigger_ess", "trigger_smd", "overlap_pass",
]


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def text_hash(values):
    return hashlib.sha256("\n".join(str(value) for value in values).encode("utf-8")).hexdigest()


def string_hash(value):
    return hashlib.sha256(str(value).encode("utf-8")).hexdigest()


def float_hash(values):
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes(order="C")).hexdigest()


def int_hash(values):
    return hashlib.sha256(np.asarray(values, dtype="<i8").tobytes(order="C")).hexdigest()


def binary_hash(values):
    return hashlib.sha256(np.asarray(values, dtype="u1").tobytes(order="C")).hexdigest()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def write_result_csv(frame, path):
    frame.to_csv(path, index=False, float_format="%.17g")


def read_result_csv(path):
    return pd.read_csv(path, float_precision="round_trip")


def verify_environment(cfg):
    observed = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "rdkit": rdkit.__version__,
        "xgboost": xgboost.__version__,
        "pyarrow": pyarrow.__version__,
        "platform": platform.platform(),
    }
    expected = {key: value for key, value in cfg["analysis_environment"].items() if key != "split_provenance_environment_is_separate"}
    require(observed == expected, f"analysis environment mismatch: {observed}")
    return observed


def clearance_value(text, label):
    values = []
    for line in text.splitlines():
        cleaned = re.sub(r"^\s*[-*]\s+", "", line.strip()).replace("**", "").replace("`", "")
        if cleaned.startswith(f"{label}:"):
            values.append(cleaned.split(":", 1)[1].strip())
    require(len(values) == 1 and values[0], f"missing or duplicate clearance field: {label}")
    return values[0]


def verify_clearance(root, mode, decision_004_sha, decision_005_sha):
    require(decision_004_sha is not None, "--decision-004-sha is required")
    decision_004_path = root / DECISION_004_PATH
    require(decision_004_path.is_file() and file_hash(decision_004_path) == decision_004_sha, "Decision 004 hash mismatch")
    text_004 = decision_004_path.read_text(encoding="utf-8")
    runner_sha = file_hash(Path(__file__))
    verifier_path = root / "src" / "verify_confirmatory.py"
    require(verifier_path.is_file(), "independent verifier is missing")
    verifier_sha = file_hash(verifier_path)
    require(clearance_value(text_004, "Decision") == "CLEAR_TO_SMOKE", "Decision 004 state mismatch")
    require(clearance_value(text_004, "Mode") == "smoke", "Decision 004 mode mismatch")
    require(clearance_value(text_004, "Runner SHA-256") == runner_sha, "Decision 004 runner hash mismatch")
    require(clearance_value(text_004, "Verifier SHA-256") == verifier_sha, "Decision 004 verifier hash mismatch")
    require(clearance_value(text_004, "Smoke output path") == SMOKE_OUTPUT_PATH, "Decision 004 smoke path mismatch")
    require(clearance_value(text_004, "Static audit path") == STATIC_AUDIT_PATH, "Decision 004 static-audit path mismatch")
    static_audit_path = root / STATIC_AUDIT_PATH
    require(static_audit_path.is_file(), "static code audit is missing")
    static_audit_sha = file_hash(static_audit_path)
    require(clearance_value(text_004, "Static audit SHA-256") == static_audit_sha, "Decision 004 static-audit hash mismatch")
    require(clearance_value(text_004, "Static audit verdict") == "CLEAR_TO_SMOKE", "Decision 004 static-audit verdict mismatch")
    static_text = static_audit_path.read_text(encoding="utf-8")
    require(clearance_value(static_text, "Verdict") == "CLEAR_TO_SMOKE", "static-audit state mismatch")
    require(clearance_value(static_text, "Runner SHA-256") == runner_sha, "static-audit runner hash mismatch")
    require(clearance_value(static_text, "Verifier SHA-256") == verifier_sha, "static-audit verifier hash mismatch")
    hashes = {
        DECISION_004_PATH: decision_004_sha,
        "src/verify_confirmatory.py": verifier_sha,
        STATIC_AUDIT_PATH: static_audit_sha,
    }
    if mode == "full":
        require(decision_005_sha is not None, "--decision-005-sha is required for full mode")
        decision_005_path = root / DECISION_005_PATH
        require(decision_005_path.is_file() and file_hash(decision_005_path) == decision_005_sha, "Decision 005 hash mismatch")
        text_005 = decision_005_path.read_text(encoding="utf-8")
        require(clearance_value(text_005, "Decision") == "CLEAR_TO_FULL", "Decision 005 state mismatch")
        require(clearance_value(text_005, "Mode") == "full", "Decision 005 mode mismatch")
        require(clearance_value(text_005, "Runner SHA-256") == runner_sha, "Decision 005 runner hash mismatch")
        require(clearance_value(text_005, "Verifier SHA-256") == verifier_sha, "Decision 005 verifier hash mismatch")
        require(clearance_value(text_005, "Decision 004 SHA-256") == decision_004_sha, "Decision 005 does not bind Decision 004")
        require(clearance_value(text_005, "Smoke status") == "PASS", "Decision 005 smoke status mismatch")
        require(clearance_value(text_005, "Independent verifier smoke status") == "PASS", "Decision 005 verifier status mismatch")
        require(clearance_value(text_005, "Smoke output path") == SMOKE_OUTPUT_PATH, "Decision 005 smoke path mismatch")
        require(clearance_value(text_005, "Independent smoke verification path") == SMOKE_VERIFICATION_PATH, "Decision 005 smoke-verification path mismatch")
        smoke_output = (root / SMOKE_OUTPUT_PATH).resolve()
        require(smoke_output.is_dir(), "authorized smoke output is missing")
        smoke_json_path = smoke_output / "smoke.json"
        smoke_manifest_path = smoke_output / "run_manifest.json"
        smoke_verification_path = root / SMOKE_VERIFICATION_PATH
        require(smoke_json_path.is_file() and smoke_manifest_path.is_file(), "authorized smoke files are missing")
        require(smoke_verification_path.is_file(), "independent smoke verification is missing")
        smoke_json_sha = file_hash(smoke_json_path)
        smoke_manifest_sha = file_hash(smoke_manifest_path)
        smoke_verification_sha = file_hash(smoke_verification_path)
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
        hashes[DECISION_005_PATH] = decision_005_sha
        hashes[f"{SMOKE_OUTPUT_PATH}/smoke.json"] = smoke_json_sha
        hashes[f"{SMOKE_OUTPUT_PATH}/run_manifest.json"] = smoke_manifest_sha
        hashes[SMOKE_VERIFICATION_PATH] = smoke_verification_sha
    return hashes


def load_recovery_amendment(root):
    paths = {
        "artifacts/confirmatory/confirm_001/protocol_amendment_002.md": PROTOCOL_AMENDMENT_002_SHA,
        "artifacts/confirmatory/confirm_001/config_amendment_002.json": CONFIG_AMENDMENT_002_SHA,
        "artifacts/confirmatory/confirm_001/failure_001.md": FAILURE_001_SHA,
        "artifacts/confirmatory/confirm_001/amendment_002_independent_audit.md": AMENDMENT_002_AUDIT_SHA,
    }
    for relative, expected in paths.items():
        path = root / relative
        require(path.is_file() and file_hash(path) == expected, f"recovery amendment hash mismatch: {relative}")
    recovery = json.loads((root / "artifacts/confirmatory/confirm_001/config_amendment_002.json").read_text(encoding="utf-8"))
    require(recovery["study_id"] == "confirm_001" and recovery["amendment_id"] == "config_amendment_002", "recovery amendment identity mismatch")
    require(recovery["base_contract"]["protocol_sha256"] == BASE_PROTOCOL_SHA, "recovery base protocol mismatch")
    require(recovery["base_contract"]["config_sha256"] == BASE_CONFIG_SHA, "recovery base config mismatch")
    require(recovery["base_contract"]["config_amendment_001_sha256"] == CONFIG_AMENDMENT_SHA, "recovery predecessor mismatch")
    require(recovery["base_contract"]["protocol_amendment_002_sha256"] == PROTOCOL_AMENDMENT_002_SHA, "recovery protocol amendment mismatch")
    require(recovery["failed_execution"]["failure_record_sha256"] == FAILURE_001_SHA, "failure-record linkage mismatch")
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
    csv_contract = recovery["csv_roundtrip_contract"]
    require(
        csv_contract["result_csv_write_float_format"] == "%.17g"
        and csv_contract["result_csv_read_float_precision"] == "round_trip"
        and csv_contract["apply_only_to_generated_result_csvs"] is True
        and csv_contract["frozen_source_csv_readers_unchanged"] is True
        and csv_contract["formula_replay_rtol"] == 0.0
        and csv_contract["formula_replay_atol"] == 1e-12
        and csv_contract["tolerance_relaxed"] is False
        and csv_contract["recovery_smoke_requires_large_ratio_fixture"] is True
        and csv_contract["recovery_smoke_minimum_ratio_magnitude"] == 100000.0,
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
        require(path.stat().st_size == record["bytes"] and file_hash(path) == record["sha256"], f"failed artifact changed: {name}")
    require(all(not (failed_dir / name).exists() for name in failed["required_absent_files"]), "failed run gained a completion artifact")
    for relative, expected in failed["failed_code_snapshots"].items():
        path = root / relative
        require(path.is_file() and file_hash(path) == expected, f"failed code snapshot mismatch: {relative}")
    for key in ["decision_004", "decision_005", "failure_record"]:
        path = root / failed[f"{key}_path"]
        require(path.is_file() and file_hash(path) == failed[f"{key}_sha256"], f"failed-state evidence mismatch: {key}")
    audit_text = (root / "artifacts/confirmatory/confirm_001/amendment_002_independent_audit.md").read_text(encoding="utf-8")
    require(clearance_value(audit_text, "Verdict") == "CLEAR_TO_IMPLEMENT_RECOVERY", "recovery amendment was not cleared")
    return recovery, paths


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
        "Recovery smoke JSON SHA-256": file_hash(smoke_json_path),
        "Recovery smoke run manifest SHA-256": file_hash(smoke_manifest_path),
        "Recovery roundtrip fixture SHA-256": file_hash(fixture_path),
        "Independent recovery smoke verification SHA-256": file_hash(verification_path),
    }

    fixture = read_result_csv(fixture_path)
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
    require(smoke_json["fold_sha256"] == int_hash(fold_contract(groups)[1]), "recovery smoke fold hash mismatch")

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


def verify_recovery_clearance(
    root, mode, decision_006_sha, decision_007_sha, recovery, base_inputs, environment,
):
    require(mode in {"recovery-smoke", "recovery-full"}, "invalid recovery mode")
    require(re.fullmatch(r"[0-9a-f]{64}", decision_006_sha or ""), "trusted Decision 006 SHA-256 is required")
    decision_006 = root / DECISION_006_PATH
    require(decision_006.is_file() and file_hash(decision_006) == decision_006_sha, "Decision 006 hash mismatch")
    text_006 = decision_006.read_text(encoding="utf-8")
    runner_sha = file_hash(Path(__file__))
    verifier_path = root / "src" / "verify_confirmatory.py"
    require(verifier_path.is_file(), "recovery verifier is missing")
    verifier_sha = file_hash(verifier_path)
    failed = recovery["failed_execution"]
    require(clearance_value(text_006, "Decision") == "CLEAR_TO_RECOVERY_SMOKE", "Decision 006 state mismatch")
    require(clearance_value(text_006, "Mode") == "recovery-smoke", "Decision 006 mode mismatch")
    require(clearance_value(text_006, "Runner SHA-256") == runner_sha, "Decision 006 runner hash mismatch")
    require(clearance_value(text_006, "Verifier SHA-256") == verifier_sha, "Decision 006 verifier hash mismatch")
    require(clearance_value(text_006, "Decision 004 SHA-256") == failed["decision_004_sha256"], "Decision 006 historical D004 mismatch")
    require(clearance_value(text_006, "Consumed Decision 005 SHA-256") == failed["decision_005_sha256"], "Decision 006 consumed D005 mismatch")
    require(clearance_value(text_006, "Failed runner SHA-256") == failed["failed_code_snapshots"]["artifacts/confirmatory/confirm_001/failure_001_code_snapshot/run_confirmatory_failed.py"], "Decision 006 failed-runner mismatch")
    require(clearance_value(text_006, "Failed verifier SHA-256") == failed["failed_code_snapshots"]["artifacts/confirmatory/confirm_001/failure_001_code_snapshot/verify_confirmatory_failed.py"], "Decision 006 failed-verifier mismatch")
    required_fields = {
        "Failed output path": failed["output_path"],
        "Failure record path": failed["failure_record_path"],
        "Failure record SHA-256": failed["failure_record_sha256"],
        "Protocol Amendment 002 path": "artifacts/confirmatory/confirm_001/protocol_amendment_002.md",
        "Protocol Amendment 002 SHA-256": PROTOCOL_AMENDMENT_002_SHA,
        "Config Amendment 002 path": "artifacts/confirmatory/confirm_001/config_amendment_002.json",
        "Config Amendment 002 SHA-256": CONFIG_AMENDMENT_002_SHA,
        "Amendment 002 audit path": "artifacts/confirmatory/confirm_001/amendment_002_independent_audit.md",
        "Amendment 002 audit SHA-256": AMENDMENT_002_AUDIT_SHA,
        "Amendment 002 audit verdict": "CLEAR_TO_IMPLEMENT_RECOVERY",
        "Recovery static audit path": RECOVERY_STATIC_AUDIT_PATH,
        "Recovery smoke output path": RECOVERY_SMOKE_OUTPUT_PATH,
        "Recovery full output path": RECOVERY_FULL_OUTPUT_PATH,
    }
    for label, expected in required_fields.items():
        require(clearance_value(text_006, label) == expected, f"Decision 006 field mismatch: {label}")
    static_path = root / RECOVERY_STATIC_AUDIT_PATH
    require(static_path.is_file(), "recovery static audit is missing")
    static_sha = file_hash(static_path)
    require(clearance_value(text_006, "Recovery static audit SHA-256") == static_sha, "Decision 006 static-audit hash mismatch")
    require(clearance_value(text_006, "Recovery static audit verdict") == "CLEAR_TO_RECOVERY_SMOKE", "Decision 006 static-audit verdict mismatch")
    static_text = static_path.read_text(encoding="utf-8")
    require(clearance_value(static_text, "Verdict") == "CLEAR_TO_RECOVERY_SMOKE", "recovery static audit state mismatch")
    require(clearance_value(static_text, "Runner SHA-256") == runner_sha, "recovery static-audit runner mismatch")
    require(clearance_value(static_text, "Verifier SHA-256") == verifier_sha, "recovery static-audit verifier mismatch")
    records = {
        failed["decision_004_path"]: failed["decision_004_sha256"],
        failed["decision_005_path"]: failed["decision_005_sha256"],
        failed["failure_record_path"]: failed["failure_record_sha256"],
        "artifacts/confirmatory/confirm_001/protocol_amendment_002.md": PROTOCOL_AMENDMENT_002_SHA,
        "artifacts/confirmatory/confirm_001/config_amendment_002.json": CONFIG_AMENDMENT_002_SHA,
        "artifacts/confirmatory/confirm_001/amendment_002_independent_audit.md": AMENDMENT_002_AUDIT_SHA,
        RECOVERY_STATIC_AUDIT_PATH: static_sha,
        DECISION_006_PATH: decision_006_sha,
        "src/verify_confirmatory.py": verifier_sha,
        **failed["failed_code_snapshots"],
    }
    expected_smoke_inputs = {**base_inputs, **records}
    if mode == "recovery-full":
        require(re.fullmatch(r"[0-9a-f]{64}", decision_007_sha or ""), "trusted Decision 007 SHA-256 is required")
        decision_007 = root / DECISION_007_PATH
        require(decision_007.is_file() and file_hash(decision_007) == decision_007_sha, "Decision 007 hash mismatch")
        text_007 = decision_007.read_text(encoding="utf-8")
        require(clearance_value(text_007, "Decision") == "CLEAR_TO_RECOVERY_FULL", "Decision 007 state mismatch")
        require(clearance_value(text_007, "Mode") == "recovery-full", "Decision 007 mode mismatch")
        require(clearance_value(text_007, "Runner SHA-256") == runner_sha, "Decision 007 runner hash mismatch")
        require(clearance_value(text_007, "Verifier SHA-256") == verifier_sha, "Decision 007 verifier hash mismatch")
        require(clearance_value(text_007, "Decision 006 SHA-256") == decision_006_sha, "Decision 007 does not bind Decision 006")
        require(clearance_value(text_007, "Recovery smoke status") == "PASS", "Decision 007 smoke status mismatch")
        require(clearance_value(text_007, "Independent recovery smoke verifier status") == "PASS", "Decision 007 independent smoke status mismatch")
        require(clearance_value(text_007, "Recovery smoke output path") == RECOVERY_SMOKE_OUTPUT_PATH, "Decision 007 smoke path mismatch")
        require(clearance_value(text_007, "Recovery full output path") == RECOVERY_FULL_OUTPUT_PATH, "Decision 007 full path mismatch")
        require(clearance_value(text_007, "Failed output path") == failed["output_path"], "Decision 007 failed path mismatch")
        hashes = validate_recovery_smoke_evidence(
            root, expected_smoke_inputs, environment, runner_sha, verifier_sha, decision_006_sha, recovery,
        )
        for label, observed in hashes.items():
            require(clearance_value(text_007, label) == observed, f"Decision 007 evidence mismatch: {label}")
        require(clearance_value(text_007, "Independent recovery smoke verification path") == RECOVERY_SMOKE_VERIFICATION_PATH, "Decision 007 verification path mismatch")
        records.update({
            DECISION_007_PATH: decision_007_sha,
            f"{RECOVERY_SMOKE_OUTPUT_PATH}/smoke.json": hashes["Recovery smoke JSON SHA-256"],
            f"{RECOVERY_SMOKE_OUTPUT_PATH}/run_manifest.json": hashes["Recovery smoke run manifest SHA-256"],
            f"{RECOVERY_SMOKE_OUTPUT_PATH}/roundtrip_fixture.csv": hashes["Recovery roundtrip fixture SHA-256"],
            RECOVERY_SMOKE_VERIFICATION_PATH: hashes["Independent recovery smoke verification SHA-256"],
        })
    return records


def verify_inputs(root, mode, decision_004_sha=None, decision_005_sha=None, decision_006_sha=None, decision_007_sha=None):
    paths = {
        "artifacts/confirmatory/confirm_001/protocol.md": BASE_PROTOCOL_SHA,
        "artifacts/confirmatory/confirm_001/config.json": BASE_CONFIG_SHA,
        "artifacts/decisions/003_go_confirm_001.md": DECISION_003_SHA,
        "artifacts/confirmatory/confirm_001/protocol_audit.md": PROTOCOL_AUDIT_SHA,
        "data/external/source_freeze_20260911/source_manifest.json": SOURCE_MANIFEST_SHA,
        "data/external/source_freeze_20260911/source_manifest_amendment_001.json": SOURCE_AMENDMENT_SHA,
        "artifacts/confirmatory/confirm_001/protocol_amendment_001.md": PROTOCOL_AMENDMENT_SHA,
        "artifacts/confirmatory/confirm_001/config_amendment_001.json": CONFIG_AMENDMENT_SHA,
        "artifacts/confirmatory/confirm_001/amendment_001_independent_audit.md": AMENDMENT_AUDIT_SHA,
    }
    cfg_path = root / "artifacts/confirmatory/confirm_001/config.json"
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    overlay = json.loads((root / "artifacts/confirmatory/confirm_001/config_amendment_001.json").read_text(encoding="utf-8"))
    source_manifest = json.loads((root / "data/external/source_freeze_20260911/source_manifest.json").read_text(encoding="utf-8"))
    source_amendment = json.loads((root / "data/external/source_freeze_20260911/source_manifest_amendment_001.json").read_text(encoding="utf-8"))
    paths[cfg["collision_audit"]["path"]] = cfg["collision_audit"]["sha256"]
    addendum = cfg["collision_audit"]["exact_claim_addendum"]
    paths[addendum["path"]] = addendum["sha256"]
    source_audit = cfg["source_freeze"]["audit_report"]
    paths[source_audit["path"]] = source_audit["sha256"]
    paths["artifacts/decisions/002_pilot_stop.md"] = cfg["development_evidence"]["pilot_decision_sha256"]
    for relative, expected in cfg["development_evidence"]["forensic_artifacts"].items():
        paths[f"artifacts/diagnostics/forensic_001/{relative}"] = expected
    derivative = cfg["source_freeze"]["biogen"]["analysis_derivative"]
    paths[derivative["path"]] = derivative["sha256"]
    construction_root = "evidence/upstream/expansion_ml_comparison"
    for filename, expected in cfg["source_freeze"]["biogen"]["split_provenance"]["construction_files_sha256"].items():
        paths[f"{construction_root}/{filename}"] = expected
    expansion = cfg["source_freeze"]["expansionrx_discovery"]
    paths.update({
        "data/raw/expansion_data_train.csv": expansion["train_sha256"],
        "data/raw/expansion_data_test.csv": expansion["test_sha256"],
        "data/raw/expansion_data_raw.csv": expansion["raw_sha256"],
    })

    def collect_local(node):
        if isinstance(node, dict):
            if "local_path" in node and "sha256" in node:
                paths[node["local_path"]] = node["sha256"]
            for value in node.values():
                collect_local(value)
        elif isinstance(node, list):
            for value in node:
                collect_local(value)

    collect_local(source_manifest)
    paths["src/run_confirmatory.py"] = file_hash(Path(__file__))
    recovery = None
    if mode in {"recovery-smoke", "recovery-full"}:
        recovery, recovery_paths = load_recovery_amendment(root)
        paths.update(recovery_paths)
        environment = verify_environment(cfg)
        paths.update(verify_recovery_clearance(
            root, mode, decision_006_sha, decision_007_sha, recovery, dict(paths), environment,
        ))
    else:
        paths.update(verify_clearance(root, mode, decision_004_sha, decision_005_sha))
    observed = {}
    for relative, expected in paths.items():
        path = root / relative
        require(path.is_file(), f"missing frozen input: {relative}")
        observed[relative] = file_hash(path)
        require(observed[relative] == expected, f"hash mismatch for {relative}")
    require(cfg["status"] == "FROZEN_PRE_RESULT" and cfg["formal_run_authorized"] is True, "base run not authorized")
    require(cfg["protocol_sha256"] == BASE_PROTOCOL_SHA, "base protocol identity mismatch")
    require(overlay["base_contract"]["config_sha256"] == BASE_CONFIG_SHA, "overlay base mismatch")
    require(overlay["source_manifest_amendment"]["sha256"] == SOURCE_AMENDMENT_SHA, "source overlay mismatch")
    require(overlay["protocol_amendment"]["sha256"] == PROTOCOL_AMENDMENT_SHA, "protocol overlay mismatch")
    require(source_amendment["base_manifest"]["sha256"] == SOURCE_MANIFEST_SHA, "source amendment base mismatch")
    require(
        {row["vector"] for row in source_amendment["corrections"]} == {
            "ordered_internal_id_tab_smiles", "ordered_cluster", "ordered_ds", "ordered_parent_OCNT_ID",
            "ordered_availability", "ordered_formal_hit",
        }
        and len(source_amendment["corrections"]) == 6,
        "unexpected source amendment correction set",
    )
    amendment_audit_text = (root / "artifacts/confirmatory/confirm_001/amendment_001_independent_audit.md").read_text(encoding="utf-8")
    require("Verdict: `CLEAR_TO_IMPLEMENT`" in amendment_audit_text, "amendment audit did not clear implementation")
    require(overlay["formal_input_precedence"] == [
        "artifacts/confirmatory/confirm_001/protocol.md",
        "artifacts/confirmatory/confirm_001/config.json",
        "data/external/source_freeze_20260911/source_manifest.json",
        "data/external/source_freeze_20260911/source_manifest_amendment_001.json",
        "artifacts/confirmatory/confirm_001/protocol_amendment_001.md",
        "artifacts/confirmatory/confirm_001/config_amendment_001.json",
    ], "formal input precedence mismatch")
    if recovery is not None:
        require(recovery["formal_input_precedence"] == [
            "artifacts/confirmatory/confirm_001/protocol.md",
            "artifacts/confirmatory/confirm_001/config.json",
            "data/external/source_freeze_20260911/source_manifest.json",
            "data/external/source_freeze_20260911/source_manifest_amendment_001.json",
            "artifacts/confirmatory/confirm_001/protocol_amendment_001.md",
            "artifacts/confirmatory/confirm_001/config_amendment_001.json",
            "artifacts/confirmatory/confirm_001/failure_001.md",
            "artifacts/confirmatory/confirm_001/protocol_amendment_002.md",
            "artifacts/confirmatory/confirm_001/config_amendment_002.json",
        ], "recovery formal input precedence mismatch")
    if mode not in {"recovery-smoke", "recovery-full"}:
        environment = verify_environment(cfg)
    return cfg, source_amendment, observed, environment, recovery


def descriptor_registry():
    registry = sorted(Descriptors._descList, key=lambda item: item[0])
    names = [name for name, _ in registry]
    require(len(names) == 217 and len(set(names)) == 217, "RDKit descriptor registry mismatch")
    require(text_hash(names) == DESCRIPTOR_NAMES_SHA, "RDKit descriptor-name hash mismatch")
    return registry


def build_features(smiles, registry):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    fps = np.empty((len(smiles), 2048), dtype=np.float32)
    desc = np.empty((len(smiles), len(registry)), dtype=np.float64)
    bitvectors = []
    molecules = []
    canonical = []
    for index, value in enumerate(smiles):
        mol = Chem.MolFromSmiles(str(value))
        require(mol is not None, f"SMILES parse failure at row {index}")
        fp = generator.GetFingerprint(mol)
        array = np.zeros(2048, dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fp, array)
        fps[index] = array
        desc[index] = [function(mol) for _, function in registry]
        bitvectors.append(fp)
        molecules.append(mol)
        canonical.append(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True))
    require(np.isfinite(desc).all(), "nonfinite RDKit descriptor")
    return {
        "fps": fps,
        "desc": desc,
        "bitvectors": bitvectors,
        "molecules": molecules,
        "canonical": np.asarray(canonical, dtype=object),
    }


def max_similarity(query, reference):
    return np.asarray([max(DataStructs.BulkTanimotoSimilarity(fp, reference)) for fp in query], dtype=float)


def scaffold_groups(molecules, source_rows):
    groups = []
    for mol, source_row in zip(molecules, source_rows):
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold.GetNumAtoms():
            value = Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=True)
            groups.append(f"SCAFFOLD::{value}")
        else:
            groups.append(f"ACYCLIC_SOURCE_ROW_{int(source_row)}")
    return np.asarray(groups, dtype=object)


def fold_contract(groups):
    split = list(GroupKFold(n_splits=3).split(np.zeros(len(groups)), groups=groups))
    fold = np.full(len(groups), -1, dtype=np.int64)
    for fold_id, (_, valid) in enumerate(split):
        fold[valid] = fold_id
    require(np.array_equal(np.unique(fold), np.arange(3)), "invalid fold assignment")
    return split, fold


def scaled_folds(X, split):
    prepared = []
    for train, valid in split:
        scaler = StandardScaler()
        prepared.append((train, valid, scaler.fit_transform(X[train]), scaler.transform(X[valid])))
    return prepared


def logistic_model():
    return LogisticRegression(
        penalty="l2",
        C=1.0,
        solver="liblinear",
        class_weight=None,
        max_iter=1000,
        random_state=20260910,
    )


def checked_fit(model, x, y, label):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning, module=r"sklearn\.linear_model\._logistic")
        warnings.simplefilter("error", ConvergenceWarning)
        fitted = model.fit(x, y)
    if isinstance(fitted, LogisticRegression):
        require(np.asarray(fitted.n_iter_).size > 0, f"missing iteration count: {label}")
        require(np.all(np.asarray(fitted.n_iter_) < fitted.max_iter), f"nuisance fit reached max_iter: {label}")
    return fitted


def fold_support(r, split):
    values = {}
    for fold_id, (train, _) in enumerate(split):
        values[f"fold_{fold_id}_train_observed"] = int(np.sum(r[train] == 1))
        values[f"fold_{fold_id}_train_unavailable"] = int(np.sum(r[train] == 0))
    values["fold_class_support"] = bool(
        all(values[f"fold_{fold_id}_train_{state}"] >= 20 for fold_id in range(3) for state in ["observed", "unavailable"])
    )
    return values


def fit_propensity(prepared, r):
    r = np.asarray(r, dtype=int)
    require(np.isin(r, [0, 1]).all(), "nonbinary availability indicator")
    p_raw = np.empty(len(r), dtype=float)
    for train, valid, x_train, x_valid in prepared:
        require(len(np.unique(r[train])) == 2, "propensity training fold has one class")
        fitted = checked_fit(logistic_model(), x_train, r[train], "propensity")
        p_raw[valid] = fitted.predict_proba(x_valid)[:, 1]
    p = np.clip(p_raw, 0.05, 0.95)
    require(np.isfinite(p_raw).all(), "nonfinite fitted propensity")
    return {
        "p_raw": p_raw,
        "p": p,
        "clip_fraction": float(np.mean((p_raw < 0.05) | (p_raw > 0.95))),
    }


def fit_loss(prepared, r, losses):
    r = np.asarray(r, dtype=int)
    mhat = np.empty_like(losses, dtype=float)
    for train, valid, x_train, x_valid in prepared:
        observed = r[train].astype(bool)
        require(observed.any(), "loss-nuisance training fold has no observed rows")
        model = checked_fit(
            Ridge(alpha=10.0, solver="lsqr", tol=1e-6),
            x_train[observed],
            losses[train][observed],
            "loss nuisance",
        )
        prediction = np.asarray(model.predict(x_valid), dtype=float)
        mhat[valid] = prediction[:, None] if prediction.ndim == 1 else prediction
    require(np.isfinite(mhat).all(), "nonfinite loss-nuisance prediction")
    return mhat


def imbalance(r, p, audit_features, feature_names):
    weights = np.asarray(r, dtype=float) / np.asarray(p, dtype=float)
    denominator = float(weights.sum())
    require(denominator > 0, "zero observed weight")
    values = {}
    zero_variance = []
    for column, name in enumerate(feature_names):
        x = audit_features[:, column]
        sd = float(np.std(x, ddof=0))
        if sd == 0:
            zero_variance.append(name)
        else:
            weighted = float(np.sum(weights * x) / denominator)
            values[name] = abs(weighted - float(np.mean(x))) / sd
    return {
        "smd": values,
        "zero_variance": zero_variance,
        "max_abs_smd": float(max(values.values())) if values else 0.0,
        "max_normalized_weight": float(np.max(weights[r.astype(bool)] / denominator)),
    }


def estimator(r, loss, mhat, p, groups, audit_features, audit_names):
    r = np.asarray(r, dtype=float)
    loss = np.asarray(loss, dtype=float)
    mhat = np.asarray(mhat, dtype=float)
    p = np.asarray(p, dtype=float)
    require(
        np.isfinite(np.column_stack([r, loss, mhat, p])).all()
        and np.isin(r, [0.0, 1.0]).all()
        and ((p > 0) & (p <= 1)).all(),
        "invalid estimator input",
    )
    weights = r / p
    weight_sum = float(weights.sum())
    require(weight_sum > 0, "no observed estimator rows")
    ipw = float(np.sum(weights * loss) / weight_sum)
    augmentation = r * (loss - mhat) / p
    pseudo = mhat + augmentation
    aipw = float(np.mean(pseudo))
    influence = pseudo - aipw
    codes = pd.factorize(np.asarray(groups, dtype=object), sort=False)[0]
    cluster_sums = np.bincount(codes, weights=influence)
    count = len(cluster_sums)
    require(count > 1, "fewer than two scaffold groups")
    se = float(np.sqrt((count / (count - 1)) * np.square(cluster_sums).sum()) / len(r))
    ess = float(weight_sum**2 / np.square(weights).sum())
    audit = imbalance(r, p, audit_features, audit_names)
    values = [ipw, aipw, se, ess, audit["max_abs_smd"], audit["max_normalized_weight"]]
    require(np.isfinite(values).all(), "nonfinite estimator output")
    return {
        "ipw": ipw,
        "aipw": aipw,
        "se": se,
        "lcl": aipw - 1.96 * se,
        "ucl": aipw + 1.96 * se,
        "ess": ess,
        "ess_fraction": ess / int(r.sum()),
        "augmentation": augmentation,
        "pseudo": pseudo,
        "influence": influence,
        **audit,
    }


def prepare_biogen(root, cfg, source_amendment):
    derivative_path = root / cfg["source_freeze"]["biogen"]["analysis_derivative"]["path"]
    official_path = root / cfg["source_freeze"]["biogen"]["official_data"]["local_path"]
    derivative = pd.read_csv(derivative_path)
    official = pd.read_csv(official_path)
    derivative["source_row_index"] = np.arange(len(derivative), dtype=np.int64)
    expected_columns = cfg["source_freeze"]["biogen"]["analysis_derivative"]["column_names"]
    require(list(derivative.columns[:-1]) == expected_columns and len(derivative) == 3521, "Biogen derivative schema mismatch")
    require(len(official) == 3521 and list(official.columns) == cfg["source_freeze"]["biogen"]["official_data"]["column_names"], "Biogen official schema mismatch")
    mapping = cfg["source_freeze"]["biogen"]["source_to_derived_mapping"]
    for source_column, derived_column in mapping.items():
        left = official[source_column]
        right = derivative[derived_column]
        if pd.api.types.is_numeric_dtype(left):
            require(np.array_equal(left.to_numpy(float), right.to_numpy(float), equal_nan=True), f"Biogen lineage mismatch: {derived_column}")
        else:
            require(left.astype(str).tolist() == right.astype(str).tolist(), f"Biogen lineage mismatch: {derived_column}")
    corrected = {row["vector"]: row["corrected_sha256"] for row in source_amendment["corrections"]}
    require(text_hash(official["Internal ID"].astype(str) + "\t" + official["SMILES"].astype(str)) == corrected["ordered_internal_id_tab_smiles"], "Biogen identity-vector hash mismatch")
    require(text_hash(derivative["cluster"]) == corrected["ordered_cluster"], "Biogen cluster-vector hash mismatch")
    require(text_hash(derivative["ds"]) == corrected["ordered_ds"], "Biogen split-vector hash mismatch")
    require(set(derivative["ds"]) == {"train", "test"}, "Biogen split label mismatch")
    require(int((derivative["ds"] == "train").sum()) == 2463 and int((derivative["ds"] == "test").sum()) == 1058, "Biogen split count mismatch")
    train_clusters = set(derivative.loc[derivative["ds"] == "train", "cluster"])
    test_clusters = set(derivative.loc[derivative["ds"] == "test", "cluster"])
    require(not (train_clusters & test_clusters), "Biogen cluster crosses split")
    registry = descriptor_registry()
    features = build_features(derivative["SMILES"].tolist(), registry)
    require(len(set(features["canonical"])) == 3521, "Biogen canonical duplicate")
    train_canonical = set(features["canonical"][derivative["ds"].eq("train").to_numpy()])
    test_canonical = set(features["canonical"][derivative["ds"].eq("test").to_numpy()])
    require(not (train_canonical & test_canonical), "canonical structure crosses split")
    descriptor_names = [name for name, _ in registry]
    descriptor_index = {name: index for index, name in enumerate(descriptor_names)}
    cohort_rows = []
    for split_name in ["all", "train", "test"]:
        split_mask = np.ones(len(derivative), dtype=bool) if split_name == "all" else derivative["ds"].eq(split_name).to_numpy()
        for endpoint in [*ENDPOINTS, *REGULAR_ANCHORS, *STRESS_ANCHORS]:
            finite = np.isfinite(pd.to_numeric(derivative[endpoint], errors="coerce").to_numpy(float)) & split_mask
            cohort_rows.append({
                "dataset": "Biogen_public_ADME", "split": split_name, "endpoint": endpoint,
                "n_rows": int(split_mask.sum()), "n_finite": int(finite.sum()),
                "n_unavailable": int(split_mask.sum() - finite.sum()), "coverage": float(finite.sum() / split_mask.sum()),
            })
    structure = {
        "rows": 3521,
        "smiles_parse_failures": 0,
        "unique_canonical_smiles": 3521,
        "canonical_train_test_overlap": 0,
        "ordered_canonical_smiles_sha256": text_hash(features["canonical"]),
        "descriptor_count": len(descriptor_names),
        "descriptor_names_sha256": text_hash(descriptor_names),
        "descriptor_names": descriptor_names,
        "all_descriptors_finite": True,
    }
    return derivative, features, descriptor_index, pd.DataFrame(cohort_rows, columns=COHORT_COLUMNS), structure


def fit_predictions(derivative, features, descriptor_index, cfg):
    prediction_rows = []
    contexts = {}
    for endpoint in ENDPOINTS:
        values = pd.to_numeric(derivative[endpoint], errors="coerce").to_numpy(float)
        train_idx = np.flatnonzero(derivative["ds"].eq("train").to_numpy() & np.isfinite(values))
        target_idx = np.flatnonzero(derivative["ds"].eq("test").to_numpy() & np.isfinite(values))
        require(np.all(np.diff(target_idx) > 0), "sentinel cohort row order mismatch")
        rf_cfg = cfg["prediction_models"]["ecfp_rf"]
        rf = RandomForestRegressor(
            n_estimators=rf_cfg["n_estimators"], criterion=rf_cfg["criterion"],
            min_samples_leaf=rf_cfg["min_samples_leaf"], max_features=rf_cfg["max_features"],
            n_jobs=rf_cfg["n_jobs"], random_state=rf_cfg["random_state"],
        ).fit(features["fps"][train_idx], values[train_idx])
        x_train = np.hstack([features["fps"][train_idx], features["desc"][train_idx]]).astype(np.float32)
        x_target = np.hstack([features["fps"][target_idx], features["desc"][target_idx]]).astype(np.float32)
        xgb_cfg = cfg["prediction_models"]["ecfp_desc_xgb"]
        xgb_model = XGBRegressor(**{key: value for key, value in xgb_cfg.items() if key != "features"}).fit(x_train, values[train_idx])
        predictions = {
            "ecfp_rf": np.asarray(rf.predict(features["fps"][target_idx]), dtype=float),
            "ecfp_desc_xgb": np.asarray(xgb_model.predict(x_target), dtype=float),
        }
        require(all(np.isfinite(value).all() and np.ptp(value) > 0 for value in predictions.values()), f"invalid prediction: {endpoint}")
        similarity = max_similarity([features["bitvectors"][i] for i in target_idx], [features["bitvectors"][i] for i in train_idx])
        low_unscaled = np.column_stack([
            features["desc"][target_idx][:, [descriptor_index[name] for name in LOW_DESCRIPTOR_NAMES]], similarity
        ])
        representations = {
            "lowdim_11": low_unscaled.astype(np.float64),
            "full_2266": np.hstack([
                features["fps"][target_idx], features["desc"][target_idx], similarity[:, None]
            ]).astype(np.float64),
        }
        groups = scaffold_groups([features["molecules"][i] for i in target_idx], target_idx)
        split, fold = fold_contract(groups)
        losses = np.column_stack([np.abs(values[target_idx] - predictions[model]) for model in MODELS])
        context = {
            "target_idx": target_idx,
            "source_row_index": target_idx.astype(np.int64),
            "cohort_row_index": np.arange(len(target_idx), dtype=np.int64),
            "molecule_id": derivative.iloc[target_idx]["Name"].astype(str).to_numpy(),
            "groups": groups,
            "scaffold_sha256": np.asarray([string_hash(value) for value in groups]),
            "fold": fold,
            "split": split,
            "representations": representations,
            "prepared": {name: scaled_folds(value, split) for name, value in representations.items()},
            "audit_features": low_unscaled,
            "audit_names": [*LOW_DESCRIPTOR_NAMES, "max_train_similarity"],
            "losses": losses,
            "predictions": predictions,
            "reference": values[target_idx],
            "availability": {
                anchor: np.isfinite(pd.to_numeric(derivative.iloc[target_idx][anchor], errors="coerce").to_numpy(float)).astype(int)
                for anchor in [*REGULAR_ANCHORS, *STRESS_ANCHORS]
            },
            "hashes": {
                "source_row_index_sha256": int_hash(target_idx),
                "scaffold_group_vector_sha256": text_hash(groups),
                "fold_sha256": int_hash(fold),
            },
        }
        contexts[endpoint] = context
        for model_index, model in enumerate(MODELS):
            for cohort_index, source_index in enumerate(target_idx):
                prediction_rows.append({
                    "endpoint": endpoint, "model": model, "source_row_index": int(source_index),
                    "cohort_row_index": cohort_index, "molecule_id": derivative.iloc[source_index]["Name"],
                    "reference": values[source_index], "prediction": predictions[model][cohort_index],
                    "loss": losses[cohort_index, model_index],
                })
    require(len(contexts["LOG_HLM"]["target_idx"]) == 929 and len(contexts["LOG_RLM"]["target_idx"]) == 922, "sentinel cohort count mismatch")
    return contexts, pd.DataFrame(prediction_rows, columns=PREDICTION_COLUMNS)


def generator_contract(context, anchor):
    r = context["availability"][anchor]
    coverage = float(np.mean(r))
    low_scaled = StandardScaler().fit_transform(context["representations"]["lowdim_11"])
    model = checked_fit(logistic_model(), low_scaled, r, "MAR generator")
    score = np.asarray(model.decision_function(low_scaled), dtype=float)
    lower, upper = -50.0, 50.0
    for _ in range(200):
        midpoint = (lower + upper) / 2
        mean_probability = float(np.mean(np.clip(expit(score + midpoint), 0.10, 0.90)))
        if mean_probability < coverage:
            lower = midpoint
        else:
            upper = midpoint
    shift = (lower + upper) / 2
    mar_probability = np.clip(expit(score + shift), 0.10, 0.90)
    error = abs(float(np.mean(mar_probability)) - coverage)
    require(error <= 1e-10 and np.isfinite(mar_probability).all(), "MAR generator coverage mismatch")
    return {
        "mcar": {"p": np.full(len(r), coverage, dtype=float), "score": np.zeros(len(r), dtype=float), "shift": 0.0, "error": 0.0},
        "mar": {"p": mar_probability, "score": score, "shift": float(shift), "error": error},
    }


def build_generator_vectors(contexts):
    generators = {}
    rows = []
    availability_rows = []
    for endpoint, context in contexts.items():
        for anchor in [*REGULAR_ANCHORS, *STRESS_ANCHORS]:
            r = context["availability"][anchor]
            availability_rows.append({
                "endpoint": endpoint,
                "anchor": anchor,
                "n_total": int(len(r)),
                "n_available": int(r.sum()),
                "n_unavailable": int(len(r) - r.sum()),
                "coverage": float(r.mean()),
                "availability_sha256": binary_hash(r),
            })
        for anchor in REGULAR_ANCHORS:
            generators[(endpoint, anchor)] = generator_contract(context, anchor)
            for mechanism in MECHANISMS:
                generator = generators[(endpoint, anchor)][mechanism]
                p_hash = float_hash(generator["p"])
                for row_index in range(len(context["source_row_index"])):
                    rows.append({
                        "endpoint": endpoint,
                        "anchor": anchor,
                        "mechanism": mechanism,
                        "source_row_index": int(context["source_row_index"][row_index]),
                        "cohort_row_index": row_index,
                        "direct_coverage": float(context["availability"][anchor].mean()),
                        "generator_score": generator["score"][row_index],
                        "intercept_shift": generator["shift"],
                        "p_oracle": generator["p"][row_index],
                        "p_oracle_sha256": p_hash,
                        "coverage_error": generator["error"],
                        "all_finite": bool(np.isfinite(generator["p"][row_index])),
                    })
    availability = pd.DataFrame(availability_rows, columns=AVAILABILITY_COLUMNS)
    vectors = pd.DataFrame(rows, columns=GENERATOR_COLUMNS)
    require(len(availability) == 8 and len(vectors) == 7404, "availability/generator grid mismatch")
    return generators, availability, vectors


def scenario_base(context, anchor, mechanism, seed, r, p_oracle, shift, support):
    return {
        "endpoint": None,
        "anchor": anchor,
        "mechanism": mechanism,
        "seed": int(seed),
        "n_total": int(len(r)),
        "n_observed": int(r.sum()),
        "observed_fraction": float(r.mean()),
        "direct_anchor_coverage": float(context["availability"][anchor].mean()),
        "generator_intercept_shift": shift,
        **context["hashes"],
        "R_sha256": binary_hash(r),
        "p_oracle_sha256": float_hash(p_oracle),
        **support,
    }


def semisynthetic_grid(contexts, generators, cfg):
    rows = []
    replay_frames = []
    replay_seed = int(cfg["integrity"]["row_replay_seed"])
    for endpoint, context in contexts.items():
        for anchor in REGULAR_ANCHORS:
            for mechanism in MECHANISMS:
                p_oracle = generators[(endpoint, anchor)][mechanism]["p"]
                shift = generators[(endpoint, anchor)][mechanism]["shift"]
                for seed in range(41001, 41101):
                    uniform = np.random.default_rng(seed).random(len(p_oracle))
                    r = (uniform < p_oracle).astype(int)
                    support = fold_support(r, context["split"])
                    propensity = {
                        representation: fit_propensity(context["prepared"][representation], r)
                        for representation in REPRESENTATIONS
                    }
                    mhat = {
                        representation: fit_loss(context["prepared"][representation], r, context["losses"])
                        for representation in REPRESENTATIONS
                    }
                    oracle_cache = {}
                    for loss_representation in REPRESENTATIONS:
                        for model_index, model in enumerate(MODELS):
                            oracle_cache[(loss_representation, model)] = estimator(
                                r,
                                context["losses"][:, model_index],
                                mhat[loss_representation][:, model_index],
                                p_oracle,
                                context["groups"],
                                context["audit_features"],
                                context["audit_names"],
                            )
                    base = scenario_base(context, anchor, mechanism, seed, r, p_oracle, shift, support)
                    base["endpoint"] = endpoint
                    for propensity_representation in REPRESENTATIONS:
                        prop = propensity[propensity_representation]
                        prop_audit = imbalance(r, prop["p"], context["audit_features"], context["audit_names"])
                        for loss_representation in REPRESENTATIONS:
                            for model_index, model in enumerate(MODELS):
                                fitted = estimator(
                                    r,
                                    context["losses"][:, model_index],
                                    mhat[loss_representation][:, model_index],
                                    prop["p"],
                                    context["groups"],
                                    context["audit_features"],
                                    context["audit_names"],
                                )
                                oracle = oracle_cache[(loss_representation, model)]
                                truth = float(context["losses"][:, model_index].mean())
                                row = {
                                    **base,
                                    "model": model,
                                    "propensity_representation": propensity_representation,
                                    "loss_representation": loss_representation,
                                    "p_raw_sha256": float_hash(prop["p_raw"]),
                                    "p_fitted_sha256": float_hash(prop["p"]),
                                    "mhat_sha256": float_hash(mhat[loss_representation][:, model_index]),
                                    "full_truth": truth,
                                    "cc": float(context["losses"][r.astype(bool), model_index].mean()),
                                    "or_mean": float(mhat[loss_representation][:, model_index].mean()),
                                    "fitted_ipw": fitted["ipw"],
                                    "fitted_aipw": fitted["aipw"],
                                    "fitted_se": fitted["se"],
                                    "fitted_lcl": fitted["lcl"],
                                    "fitted_ucl": fitted["ucl"],
                                    "fitted_covered": bool(fitted["lcl"] <= truth <= fitted["ucl"]),
                                    "fitted_ess": fitted["ess"],
                                    "fitted_ess_fraction": fitted["ess_fraction"],
                                    "fitted_clip_fraction": prop["clip_fraction"],
                                    "fitted_max_abs_smd": prop_audit["max_abs_smd"],
                                    "fitted_max_normalized_weight": prop_audit["max_normalized_weight"],
                                    "fitted_smd_json": json.dumps(prop_audit["smd"], sort_keys=True),
                                    "oracle_ipw": oracle["ipw"],
                                    "oracle_aipw": oracle["aipw"],
                                    "oracle_se": oracle["se"],
                                    "oracle_lcl": oracle["lcl"],
                                    "oracle_ucl": oracle["ucl"],
                                    "oracle_covered": bool(oracle["lcl"] <= truth <= oracle["ucl"]),
                                    "oracle_ess": oracle["ess"],
                                    "oracle_ess_fraction": oracle["ess_fraction"],
                                    "oracle_max_abs_smd": oracle["max_abs_smd"],
                                    "oracle_max_normalized_weight": oracle["max_normalized_weight"],
                                    "oracle_smd_json": json.dumps(oracle["smd"], sort_keys=True),
                                    "fold_class_support": support["fold_class_support"],
                                    "all_finite": True,
                                    "cc_abs_error": abs(float(context["losses"][r.astype(bool), model_index].mean()) - truth),
                                    "fitted_aipw_abs_error": abs(fitted["aipw"] - truth),
                                    "oracle_ipw_abs_error": abs(oracle["ipw"] - truth),
                                    "oracle_aipw_abs_error": abs(oracle["aipw"] - truth),
                                }
                                rows.append(row)
                                if seed == replay_seed:
                                    replay = pd.DataFrame({
                                        "endpoint": endpoint,
                                        "anchor": anchor,
                                        "mechanism": mechanism,
                                        "seed": seed,
                                        "model": model,
                                        "propensity_representation": propensity_representation,
                                        "loss_representation": loss_representation,
                                        "source_row_index": context["source_row_index"],
                                        "cohort_row_index": context["cohort_row_index"],
                                        "molecule_id": context["molecule_id"],
                                        "scaffold_group": context["groups"],
                                        "scaffold_sha256": context["scaffold_sha256"],
                                        "fold": context["fold"],
                                        "R": r,
                                        "loss": context["losses"][:, model_index],
                                        "mhat": mhat[loss_representation][:, model_index],
                                        "p_raw": prop["p_raw"],
                                        "p_fitted": prop["p"],
                                        "p_oracle": p_oracle,
                                    })
                                    for audit_index, audit_name in enumerate(context["audit_names"]):
                                        replay[f"audit_{audit_name}"] = context["audit_features"][:, audit_index]
                                    replay["augmentation_fitted"] = fitted["augmentation"]
                                    replay["pseudo_outcome_fitted"] = fitted["pseudo"]
                                    replay["influence_fitted"] = fitted["influence"]
                                    replay["augmentation_oracle"] = oracle["augmentation"]
                                    replay["pseudo_outcome_oracle"] = oracle["pseudo"]
                                    replay["influence_oracle"] = oracle["influence"]
                                    replay_frames.append(replay[REPLAY_COLUMNS])
                print(f"SEMISYNTHETIC {endpoint} {anchor} {mechanism}", flush=True)
    table = pd.DataFrame(rows, columns=SEMISYNTHETIC_COLUMNS)
    replay = pd.concat(replay_frames, ignore_index=True)
    require(len(table) == 6400 and len(replay) == 59232, "semisynthetic/replay row grid mismatch")
    return table, replay


def overlap_eligible(support, prop, prop_audit):
    return bool(
        support["fold_class_support"]
        and prop["clip_fraction"] <= 0.10
        and prop_audit["ess_fraction"] >= 0.50
        and prop_audit["max_abs_smd"] <= 0.10
        and prop_audit["max_normalized_weight"] <= 0.05
    )


def actual_mask_grid(contexts):
    rows = []
    for endpoint, context in contexts.items():
        for anchor in [*REGULAR_ANCHORS, *STRESS_ANCHORS]:
            r = context["availability"][anchor]
            support = fold_support(r, context["split"])
            propensity = {
                representation: fit_propensity(context["prepared"][representation], r)
                for representation in REPRESENTATIONS
            }
            mhat = {
                representation: fit_loss(context["prepared"][representation], r, context["losses"])
                for representation in REPRESENTATIONS
            }
            for propensity_representation in REPRESENTATIONS:
                prop = propensity[propensity_representation]
                audit = imbalance(r, prop["p"], context["audit_features"], context["audit_names"])
                weights = r / prop["p"]
                ess = float(weights.sum() ** 2 / np.square(weights).sum())
                audit["ess_fraction"] = ess / int(r.sum())
                eligible = overlap_eligible(support, prop, audit)
                for loss_representation in REPRESENTATIONS:
                    for model_index, model in enumerate(MODELS):
                        fitted = estimator(
                            r,
                            context["losses"][:, model_index],
                            mhat[loss_representation][:, model_index],
                            prop["p"],
                            context["groups"],
                            context["audit_features"],
                            context["audit_names"],
                        )
                        truth = float(context["losses"][:, model_index].mean())
                        cc = float(context["losses"][r.astype(bool), model_index].mean())
                        rows.append({
                            "endpoint": endpoint,
                            "anchor": anchor,
                            "mask_type": "direct_empirical",
                            "model": model,
                            "propensity_representation": propensity_representation,
                            "loss_representation": loss_representation,
                            "n_total": int(len(r)),
                            "n_observed": int(r.sum()),
                            "observed_fraction": float(r.mean()),
                            **context["hashes"],
                            "R_sha256": binary_hash(r),
                            "p_raw_sha256": float_hash(prop["p_raw"]),
                            "p_fitted_sha256": float_hash(prop["p"]),
                            "mhat_sha256": float_hash(mhat[loss_representation][:, model_index]),
                            "full_truth": truth,
                            "cc": cc,
                            "or_mean": float(mhat[loss_representation][:, model_index].mean()),
                            "fitted_ipw": fitted["ipw"],
                            "fitted_aipw": fitted["aipw"],
                            "fitted_se": fitted["se"],
                            "fitted_lcl": fitted["lcl"],
                            "fitted_ucl": fitted["ucl"],
                            "fitted_covered": bool(fitted["lcl"] <= truth <= fitted["ucl"]),
                            "fitted_ess": fitted["ess"],
                            "fitted_ess_fraction": fitted["ess_fraction"],
                            "fitted_clip_fraction": prop["clip_fraction"],
                            "fitted_max_abs_smd": audit["max_abs_smd"],
                            "fitted_max_normalized_weight": audit["max_normalized_weight"],
                            "fitted_smd_json": json.dumps(audit["smd"], sort_keys=True),
                            **{column: support[column] for column in FOLD_SUPPORT_COLUMNS},
                            "fold_class_support": support["fold_class_support"],
                            "all_finite": True,
                            "overlap_eligible": eligible,
                            "overlap_label": "ELIGIBLE_OVERLAP" if eligible else "INELIGIBLE_OVERLAP",
                            "cc_abs_error": abs(cc - truth),
                            "fitted_aipw_abs_error": abs(fitted["aipw"] - truth),
                        })
            print(f"DIRECT {endpoint} {anchor}", flush=True)
    table = pd.DataFrame(rows, columns=ACTUAL_COLUMNS)
    require(len(table) == 64, "actual-mask row grid mismatch")
    return table


def prepare_pxr(root, cfg, source_amendment):
    pxr_cfg = cfg["source_freeze"]["pxr"]
    single_path = root / pxr_cfg["files"]["single_concentration_train"]["local_path"]
    dose_path = root / pxr_cfg["files"]["dose_response_train"]["local_path"]
    single = pd.read_csv(single_path)
    dose = pd.read_csv(dose_path)
    require(list(single.columns) == pxr_cfg["files"]["single_concentration_train"]["column_names"] and len(single) == 21003, "PXR screen schema mismatch")
    require(list(dose.columns) == pxr_cfg["files"]["dose_response_train"]["column_names"] and len(dose) == 4139, "PXR dose-response schema mismatch")
    require(not dose["OCNT_ID"].duplicated().any(), "duplicate PXR dose-response ID")
    single["source_row_index"] = np.arange(len(single), dtype=np.int64)
    parent = single[
        single["compound_class"].eq("Library")
        & single["concentration_M"].eq(8.251e-06)
    ].sort_values("source_row_index").copy()
    parent["cohort_row_index"] = np.arange(len(parent), dtype=np.int64)
    require(len(parent) == 10747 and not parent["OCNT_ID"].duplicated().any(), "PXR parent mismatch")
    followup_ids = set(dose["OCNT_ID"].astype(str))
    r = parent["OCNT_ID"].astype(str).isin(followup_ids).to_numpy(int)
    formal_hit = ((parent["log2_fc_estimate"] > 1.0) & (parent["fdr_bh"] < 0.05)).to_numpy(int)
    corrected = {row["vector"]: row["corrected_sha256"] for row in source_amendment["corrections"]}
    require(text_hash(parent["OCNT_ID"].astype(str)) == corrected["ordered_parent_OCNT_ID"], "PXR parent-ID vector mismatch")
    require(binary_hash(r) == corrected["ordered_availability"], "PXR availability vector mismatch")
    require(binary_hash(formal_hit) == corrected["ordered_formal_hit"], "PXR formal-hit vector mismatch")
    require(int(r.sum()) == 2722 and int(formal_hit.sum()) == 873 and int(r[formal_hit.astype(bool)].sum()) == 869, "PXR linkage/stratum count mismatch")
    registry = descriptor_registry()
    registry_map = dict(registry)
    features = build_features(parent["SMILES"].tolist(), [(name, registry_map[name]) for name in LOW_DESCRIPTOR_NAMES])
    descriptor_10 = features["desc"]
    require(len(set(features["canonical"])) == 10747 and np.isfinite(descriptor_10).all(), "PXR structure audit mismatch")
    source_rows = parent["source_row_index"].to_numpy(np.int64)
    groups = scaffold_groups(features["molecules"], source_rows)
    split, fold = fold_contract(groups)
    fdr_score = -np.log10(np.maximum(parent["fdr_bh"].to_numpy(float), 1e-12))
    require(np.isfinite(fdr_score).all() and fdr_score.min() >= 0 and fdr_score.max() <= 12, "PXR fdr_score mismatch")
    specifications = {
        "structure_10": descriptor_10,
        "structure_screen_12": np.column_stack([descriptor_10, parent["log2_fc_estimate"].to_numpy(float), fdr_score]),
    }
    propensity = {
        name: fit_propensity(scaled_folds(X, split), r) for name, X in specifications.items()
    }
    dose_by_id = dose.set_index("OCNT_ID")
    dose_smiles = parent["OCNT_ID"].map(dose_by_id["SMILES"])
    dose_batch = parent["OCNT_ID"].map(dose_by_id["OCNT Batch"])
    linked = r.astype(bool)
    exact_match = pd.Series(pd.NA, index=parent.index, dtype="boolean")
    batch_match = pd.Series(pd.NA, index=parent.index, dtype="boolean")
    exact_match.loc[linked] = parent.loc[linked, "SMILES"].to_numpy() == dose_smiles.loc[linked].to_numpy()
    batch_match.loc[linked] = parent.loc[linked, "OCNT Batch"].to_numpy() == dose_batch.loc[linked].to_numpy()
    require(int(exact_match.fillna(False).sum()) == 2375, "PXR exact-SMILES join count mismatch")
    require(int(batch_match.fillna(False).sum()) == 2700, "PXR batch join count mismatch")
    combined = pd.concat([
        single[["OCNT_ID", "SMILES", "OCNT Batch"]],
        dose[["OCNT_ID", "SMILES", "OCNT Batch"]],
    ], ignore_index=True)
    cardinality = combined.groupby("OCNT_ID").agg(
        cross_file_exact_smiles_nunique=("SMILES", "nunique"),
        cross_file_batch_nunique=("OCNT Batch", "nunique"),
    )
    require(int((cardinality["cross_file_exact_smiles_nunique"] > 1).sum()) == 351 and int((cardinality["cross_file_batch_nunique"] > 1).sum()) == 23, "PXR cross-file cardinality mismatch")
    linkage = pd.DataFrame({
        "source_row_index": source_rows,
        "cohort_row_index": parent["cohort_row_index"].to_numpy(np.int64),
        "OCNT_ID": parent["OCNT_ID"].astype(str).to_numpy(),
        "R_followup": r,
        "formal_hit": formal_hit,
        "smiles_parse_ok": True,
        "descriptor_10_all_finite": np.isfinite(descriptor_10).all(axis=1),
        "canonical_smiles_sha256": [string_hash(value) for value in features["canonical"]],
        "murcko_group_sha256": [string_hash(value) for value in groups],
        "fold": fold,
    })
    for index, name in enumerate(LOW_DESCRIPTOR_NAMES):
        linkage[name] = descriptor_10[:, index]
    linkage["log2_fc_estimate"] = parent["log2_fc_estimate"].to_numpy(float)
    linkage["fdr_bh"] = parent["fdr_bh"].to_numpy(float)
    linkage["fdr_score"] = fdr_score
    for name, result in propensity.items():
        linkage[f"p_raw_{name}"] = result["p_raw"]
        linkage[f"p_clipped_{name}"] = result["p"]
        linkage[f"clipped_{name}"] = (result["p_raw"] < 0.05) | (result["p_raw"] > 0.95)
    linkage["exact_smiles_match_if_linked"] = exact_match.to_numpy()
    linkage["batch_match_if_linked"] = batch_match.to_numpy()
    linkage["cross_file_exact_smiles_nunique"] = parent["OCNT_ID"].map(cardinality["cross_file_exact_smiles_nunique"]).to_numpy(int)
    linkage["cross_file_batch_nunique"] = parent["OCNT_ID"].map(cardinality["cross_file_batch_nunique"]).to_numpy(int)
    linkage = linkage[PXR_LINKAGE_COLUMNS]
    require(len(linkage) == 10747, "PXR linkage row count mismatch")
    hashes = {
        "ordered_source_row_index": int_hash(source_rows),
        "ordered_murcko_group": text_hash(groups),
        "ordered_fold": int_hash(fold),
        "ordered_canonical_smiles": text_hash(features["canonical"]),
        "ordered_parent_OCNT_ID": text_hash(parent["OCNT_ID"].astype(str)),
        "ordered_availability": binary_hash(r),
        "ordered_formal_hit": binary_hash(formal_hit),
    }
    structure = {
        "rows": 10747,
        "smiles_parse_failures": 0,
        "unique_canonical_smiles": 10747,
        "descriptor_10_all_finite_rows": 10747,
        "ordered_canonical_smiles_sha256": hashes["ordered_canonical_smiles"],
        "ordered_murcko_group_sha256": hashes["ordered_murcko_group"],
    }
    return linkage, specifications, split, hashes, structure


def summarize_pxr(linkage_path, cfg, hashes):
    linkage = read_result_csv(linkage_path)
    require(list(linkage.columns) == PXR_LINKAGE_COLUMNS and len(linkage) == 10747, "reloaded PXR linkage schema mismatch")
    r = linkage["R_followup"].to_numpy(int)
    formal_hit = linkage["formal_hit"].to_numpy(int)
    fold = linkage["fold"].to_numpy(int)
    split = [(np.flatnonzero(fold != value), np.flatnonzero(fold == value)) for value in range(3)]
    strata_rows = []
    for value in [1, 0]:
        selected = formal_hit == value
        total = int(selected.sum())
        available = int(r[selected].sum())
        strata_rows.append({
            "formal_hit": value,
            "n_total": total,
            "n_followup_available": available,
            "n_followup_unavailable": total - available,
            "advancement_probability": available / total,
        })
    strata = pd.DataFrame(strata_rows, columns=PXR_STRATA_COLUMNS)
    require(
        strata.set_index("formal_hit").loc[1, ["n_total", "n_followup_available", "n_followup_unavailable"]].tolist() == [873, 869, 4]
        and strata.set_index("formal_hit").loc[0, ["n_total", "n_followup_available", "n_followup_unavailable"]].tolist() == [9874, 1853, 8021],
        "PXR selection-stratum mismatch",
    )
    rows = []
    thresholds = cfg["pxr_boundary"]["abstention_triggers_any"]
    specs = {
        "structure_10": LOW_DESCRIPTOR_NAMES,
        "structure_screen_12": [*LOW_DESCRIPTOR_NAMES, "log2_fc_estimate", "fdr_score"],
    }
    support = fold_support(r, split)
    for specification, names in specs.items():
        p_raw = linkage[f"p_raw_{specification}"].to_numpy(float)
        p = linkage[f"p_clipped_{specification}"].to_numpy(float)
        X = linkage[names].to_numpy(float)
        audit = imbalance(r, p, X, names)
        weights = r / p
        ess = float(weights.sum() ** 2 / np.square(weights).sum())
        clip_fraction = float(np.mean((p_raw < 0.05) | (p_raw > 0.95)))
        trigger_clip = clip_fraction > thresholds["maximum_clip_fraction"]
        trigger_ess = ess / int(r.sum()) < thresholds["minimum_ess_over_observed"]
        trigger_smd = audit["max_abs_smd"] > thresholds["maximum_max_abs_smd"]
        rows.append({
            "specification": f"pxr_{specification}",
            "n_total": len(r),
            "n_observed": int(r.sum()),
            "n_unavailable": int(len(r) - r.sum()),
            "coverage": float(r.mean()),
            "clip_fraction": clip_fraction,
            "ess": ess,
            "ess_fraction": ess / int(r.sum()),
            "max_abs_smd": audit["max_abs_smd"],
            "max_normalized_weight": audit["max_normalized_weight"],
            "smd_json": json.dumps(audit["smd"], sort_keys=True),
            "zero_variance_features_json": json.dumps(audit["zero_variance"]),
            **{column: support[column] for column in FOLD_SUPPORT_COLUMNS},
            "fold_class_support": support["fold_class_support"],
            "all_finite": bool(np.isfinite(np.column_stack([p_raw, p, X])).all()),
            "source_row_index_sha256": int_hash(linkage["source_row_index"].to_numpy(int)),
            "murcko_group_vector_sha256": hashes["ordered_murcko_group"],
            "fold_sha256": int_hash(fold),
            "p_raw_sha256": float_hash(p_raw),
            "p_clipped_sha256": float_hash(p),
            "availability_sha256": binary_hash(r),
            "formal_hit_sha256": binary_hash(formal_hit),
            "trigger_clip": trigger_clip,
            "trigger_ess": trigger_ess,
            "trigger_smd": trigger_smd,
            "overlap_pass": bool(not (trigger_clip or trigger_ess or trigger_smd) and support["fold_class_support"]),
        })
    overlap = pd.DataFrame(rows, columns=PXR_OVERLAP_COLUMNS)
    minimum_n = cfg["pxr_boundary"]["abstention_triggers_any"]["selection_stratum_minimum_n"]
    lower, upper = cfg["pxr_boundary"]["abstention_triggers_any"]["selection_stratum_extreme_probability"]
    extreme = bool(((strata["n_total"] >= minimum_n) & ((strata["advancement_probability"] <= lower) | (strata["advancement_probability"] >= upper))).any())
    propensity_triggers = [
        {
            "specification": row["specification"],
            "clip": bool(row["trigger_clip"]),
            "ess": bool(row["trigger_ess"]),
            "smd": bool(row["trigger_smd"]),
            "any": bool(row["trigger_clip"] or row["trigger_ess"] or row["trigger_smd"]),
        }
        for row in overlap.to_dict("records")
    ]
    boundary = {
        "study_id": "confirm_001",
        "decision": "ABSTAIN_NONIDENTIFIABLE_FULL_LIBRARY_DOSE_RESPONSE",
        "selection_strata_extreme": extreme,
        "propensity_triggers": propensity_triggers,
        "documentary_triggers": {
            "missing_documented_promotion_branch_or_selection_variable": True,
            "public_files_do_not_reconstruct_official_parent_and_promoted_populations": True,
            "missing_nonadvanced_dose_response_truth_without_randomized_audit_or_identification_design": True,
        },
        "full_library_point_estimate_allowed": False,
        "allowed_output": cfg["pxr_boundary"]["allowed_output"],
    }
    return strata, overlap, boundary


def factorial_contrasts(scenarios):
    rows = []
    medians = scenarios.groupby(
        ["endpoint", "anchor", "mechanism", "model", "propensity_representation", "loss_representation"]
    )["fitted_aipw_abs_error"].median()
    for endpoint in ENDPOINTS:
        for anchor in REGULAR_ANCHORS:
            for mechanism in MECHANISMS:
                for model in MODELS:
                    value = lambda p, m: float(medians.loc[(endpoint, anchor, mechanism, model, p, m)])
                    low_low = value("lowdim_11", "lowdim_11")
                    full_low = value("full_2266", "lowdim_11")
                    low_full = value("lowdim_11", "full_2266")
                    full_full = value("full_2266", "full_2266")
                    contrasts = {
                        "P": 0.5 * ((full_low - low_low) + (full_full - low_full)),
                        "M": 0.5 * ((low_full - low_low) + (full_full - full_low)),
                        "P_x_M": (full_full - low_full) - (full_low - low_low),
                    }
                    for name, estimate in contrasts.items():
                        rows.append({
                            "row_type": "factorial_contrast", "endpoint": endpoint, "anchor": anchor,
                            "mechanism": mechanism, "model": model, "loss_representation": None,
                            "contrast": name, "estimate": estimate, "numerator": None, "denominator": None,
                            "ratio": None, "rows": 400, "status": "exploratory_only",
                        })
    selected = scenarios[scenarios["propensity_representation"].eq("full_2266")]
    for endpoint in ENDPOINTS:
        for mechanism in MECHANISMS:
            for loss_representation in REPRESENTATIONS:
                part = selected[
                    selected["endpoint"].eq(endpoint)
                    & selected["mechanism"].eq(mechanism)
                    & selected["loss_representation"].eq(loss_representation)
                ]
                require(len(part) == 400, "H3 aggregation grid mismatch")
                numerator = float(part["oracle_aipw_abs_error"].median())
                denominator = float(part["oracle_ipw_abs_error"].median())
                rows.append({
                    "row_type": "augmentation_ratio", "endpoint": endpoint, "anchor": None,
                    "mechanism": mechanism, "model": None, "loss_representation": loss_representation,
                    "contrast": None, "estimate": None, "numerator": numerator, "denominator": denominator,
                    "ratio": numerator / denominator if denominator > 0 else None, "rows": len(part),
                    "status": "exploratory_only_no_gate",
                })
    table = pd.DataFrame(rows, columns=CONTRAST_COLUMNS)
    require(len(table) == 56, "factorial contrast row count mismatch")
    return table


def seed_value_matrix(frame, value_column):
    ordered = frame.sort_values(["seed", "anchor", "model"])
    require(
        len(ordered) == 400
        and ordered.groupby("seed").size().eq(4).all()
        and ordered["seed"].drop_duplicates().tolist() == list(range(41001, 41101)),
        "bootstrap seed grid mismatch",
    )
    return ordered[value_column].to_numpy(float).reshape(100, 4)


def bootstrap_intervals(scenarios, cfg):
    boot_cfg = cfg["confirmatory_hypotheses"]["bootstrap"]
    seed_indices = np.random.default_rng(52001).integers(
        0, 100, size=(10000, 100), endpoint=False, dtype=np.int64
    )
    seed_hash = int_hash(seed_indices)
    require(seed_hash == boot_cfg["seed_index_matrix_sha256_numpy_1_26_4"], "bootstrap seed-index hash mismatch")
    rows = []
    for endpoint in ENDPOINTS:
        for mechanism in MECHANISMS:
            base = scenarios[
                scenarios["endpoint"].eq(endpoint)
                & scenarios["mechanism"].eq(mechanism)
                & scenarios["loss_representation"].eq("lowdim_11")
            ]
            full = base[base["propensity_representation"].eq("full_2266")]
            low = base[base["propensity_representation"].eq("lowdim_11")]
            support_valid = bool(full["fold_class_support"].all() and low["fold_class_support"].all())
            finite_valid = bool(full["all_finite"].all() and low["all_finite"].all())
            definitions = [
                (
                    "H1_PROPENSITY_INSTABILITY_CONFIRMED",
                    "oracle_to_fitted_aipw_error_ratio",
                    seed_value_matrix(full, "oracle_aipw_abs_error"),
                    seed_value_matrix(full, "fitted_aipw_abs_error"),
                    0.50,
                    "less_or_equal",
                    "ratio",
                ),
                (
                    "H2_HIGH_DIMENSION_INSTABILITY_CONFIRMED",
                    "low_to_full_fitted_aipw_error_ratio",
                    seed_value_matrix(low, "fitted_aipw_abs_error"),
                    seed_value_matrix(full, "fitted_aipw_abs_error"),
                    0.50,
                    "less_or_equal",
                    "ratio",
                ),
                (
                    "H2_HIGH_DIMENSION_INSTABILITY_CONFIRMED",
                    "full_minus_low_clip_fraction",
                    seed_value_matrix(full, "fitted_clip_fraction"),
                    seed_value_matrix(low, "fitted_clip_fraction"),
                    0.20,
                    "greater_or_equal",
                    "difference",
                ),
            ]
            for hypothesis, statistic, numerator, denominator, threshold, direction, operation in definitions:
                numerator_point = float(np.median(numerator))
                denominator_point = float(np.median(denominator))
                numerator_boot = np.median(numerator[seed_indices], axis=(1, 2))
                denominator_boot = np.median(denominator[seed_indices], axis=(1, 2))
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
                    point = numerator_point - denominator_point
                    require(denominator_valid and np.isfinite(values).all(), "nonfinite bootstrap difference")
                    lower, upper = np.quantile(values, [0.025, 0.975], method="linear")
                    status = "OK"
                if not finite_valid:
                    status = "NONADJUDICATIVE_NONFINITE"
                elif not support_valid:
                    status = "NONADJUDICATIVE_FOLD_SUPPORT"
                point_pass = bool(
                    denominator_valid and support_valid and finite_valid
                    and (point <= threshold if direction == "less_or_equal" else point >= threshold)
                )
                rows.append({
                    "hypothesis": hypothesis, "statistic": statistic, "endpoint": endpoint,
                    "mechanism": mechanism, "point_estimate": point, "lower": float(lower), "upper": float(upper),
                    "threshold": threshold, "direction": direction, "point_pass": point_pass,
                    "denominator_valid": denominator_valid, "resamples": 10000,
                    "seed_index_matrix_sha256": seed_hash,
                    "status": status,
                })
    table = pd.DataFrame(rows, columns=BOOTSTRAP_COLUMNS)
    require(len(table) == 12, "bootstrap row count mismatch")
    return table, seed_hash


def admission_and_translation(scenarios, actual, cfg):
    candidate = scenarios[
        scenarios["propensity_representation"].eq("lowdim_11")
        & scenarios["loss_representation"].eq("lowdim_11")
    ]
    admission_cells = []
    criteria = cfg["admission_gate"]["criteria"]
    for endpoint in ENDPOINTS:
        for anchor in REGULAR_ANCHORS:
            for model in MODELS:
                for mechanism in MECHANISMS:
                    part = candidate[
                        candidate["endpoint"].eq(endpoint)
                        & candidate["anchor"].eq(anchor)
                        & candidate["model"].eq(model)
                        & candidate["mechanism"].eq(mechanism)
                    ]
                    require(len(part) == 100, "admission cell grid mismatch")
                    b_cc = float(part["cc_abs_error"].median())
                    b_aipw = float(part["fitted_aipw_abs_error"].median())
                    denominator_valid = bool(np.isfinite(b_cc) and b_cc > 0)
                    ratio = b_aipw / b_cc if denominator_valid else None
                    threshold = (
                        criteria["maximum_median_aipw_to_complete_case_absolute_estimation_error_ratio_mcar"]
                        if mechanism == "mcar"
                        else criteria["maximum_median_aipw_to_complete_case_absolute_estimation_error_ratio_mar"]
                    )
                    coverage = float(part["fitted_covered"].mean())
                    ess_q05 = float(np.quantile(part["fitted_ess_fraction"], 0.05, method="linear"))
                    clip_median = float(part["fitted_clip_fraction"].median())
                    clip_q95 = float(np.quantile(part["fitted_clip_fraction"], 0.95, method="linear"))
                    smd_q95 = float(np.quantile(part["fitted_max_abs_smd"], 0.95, method="linear"))
                    support = bool(part["fold_class_support"].all())
                    finite = bool(part["all_finite"].all() and np.isfinite(part[[
                        "full_truth", "cc", "fitted_aipw", "fitted_se", "fitted_ess_fraction",
                        "fitted_clip_fraction", "fitted_max_abs_smd",
                    ]].to_numpy(float)).all())
                    passed = bool(
                        denominator_valid and ratio <= threshold
                        and coverage >= criteria["minimum_empirical_95pct_coverage"]
                        and ess_q05 >= criteria["minimum_fifth_percentile_ess_over_observed"]
                        and clip_median <= criteria["maximum_median_clip_fraction"]
                        and clip_q95 <= criteria["maximum_95th_percentile_clip_fraction"]
                        and smd_q95 <= criteria["maximum_95th_percentile_max_abs_smd"]
                        and support and finite
                    )
                    admission_cells.append({
                        "endpoint": endpoint, "anchor": anchor, "model": model, "mechanism": mechanism,
                        "rows": len(part), "B_CC": b_cc, "B_AIPW": b_aipw, "error_ratio": ratio,
                        "error_threshold": threshold, "denominator_valid": denominator_valid,
                        "coverage": coverage, "ess_fraction_q05": ess_q05, "clip_fraction_median": clip_median,
                        "clip_fraction_q95": clip_q95, "max_abs_smd_q95": smd_q95,
                        "fold_class_support": support, "all_finite": finite, "pass": passed,
                    })
    admission_pass = bool(len(admission_cells) == 16 and all(cell["pass"] for cell in admission_cells))
    direct_candidate = actual[
        actual["anchor"].isin(REGULAR_ANCHORS)
        & actual["propensity_representation"].eq("lowdim_11")
        & actual["loss_representation"].eq("lowdim_11")
    ]
    overlap_cells = []
    for row in direct_candidate.sort_values(["endpoint", "anchor", "model"]).to_dict("records"):
        overlap_cells.append({
            "endpoint": row["endpoint"], "anchor": row["anchor"], "model": row["model"],
            "fold_class_support": bool(row["fold_class_support"]), "clip_fraction": row["fitted_clip_fraction"],
            "ess_fraction": row["fitted_ess_fraction"], "max_abs_smd": row["fitted_max_abs_smd"],
            "max_normalized_weight": row["fitted_max_normalized_weight"], "pass": bool(row["overlap_eligible"]),
            "label": row["overlap_label"],
        })
    require(len(overlap_cells) == 8, "direct overlap cell grid mismatch")
    overlap_pass = bool(all(cell["pass"] for cell in overlap_cells))
    cc_median = float(direct_candidate["cc_abs_error"].median())
    aipw_median = float(direct_candidate["fitted_aipw_abs_error"].median())
    denominator_valid = bool(np.isfinite(cc_median) and cc_median > 0)
    reduction = 1 - aipw_median / cc_median if denominator_valid else None
    better = int((direct_candidate["fitted_aipw_abs_error"] < direct_candidate["cc_abs_error"]).sum())
    if not admission_pass:
        translation_status = "NOT_ADJUDICATED_ADMISSION_FAILED"
        translation_pass = False
    elif not overlap_pass:
        translation_status = "NOT_ADJUDICATED_OVERLAP_FAILED"
        translation_pass = False
    elif not denominator_valid:
        translation_status = "NONADJUDICATIVE_FAIL"
        translation_pass = False
    else:
        translation_pass = bool(reduction >= 0.25 and better >= 6)
        translation_status = "PASS" if translation_pass else "FAIL"
    return {
        "study_id": "confirm_001",
        "candidate_specification": {"propensity_representation": "lowdim_11", "loss_representation": "lowdim_11"},
        "known_truth_admission": {"pass": admission_pass, "required_cells": 16, "cells": admission_cells},
        "direct_mask_overlap": {"pass": overlap_pass, "required_cells": 8, "cells": overlap_cells},
        "direct_mask_translation": {
            "status": translation_status, "pass": translation_pass, "required_cells": 8,
            "cc_median_abs_error": cc_median, "aipw_median_abs_error": aipw_median,
            "denominator_valid": denominator_valid, "error_reduction": reduction,
            "cells_aipw_better": better, "required_aipw_better": 6,
        },
    }


def fold_table(contexts):
    frames = []
    for endpoint, context in contexts.items():
        frame = pd.DataFrame({
            "endpoint": endpoint,
            "source_row_index": context["source_row_index"],
            "cohort_row_index": context["cohort_row_index"],
            "molecule_id": context["molecule_id"],
            "scaffold_group": context["groups"],
            "scaffold_sha256": context["scaffold_sha256"],
            "fold": context["fold"],
            "source_row_index_sha256": context["hashes"]["source_row_index_sha256"],
            "scaffold_group_vector_sha256": context["hashes"]["scaffold_group_vector_sha256"],
            "fold_sha256": context["hashes"]["fold_sha256"],
        })
        frames.append(frame[FOLD_COLUMNS])
    table = pd.concat(frames, ignore_index=True)
    require(len(table) == 1851, "fold-assignment row count mismatch")
    return table


def close(left, right, tolerance, label):
    left = float(left)
    right = float(right)
    require(np.isfinite(left) and np.isfinite(right), f"nonfinite {label}")
    require(abs(left - right) <= tolerance, f"{label} mismatch: {left} != {right}")


def vector_close(left, right, tolerance, label):
    left = np.asarray(left)
    right = np.asarray(right)
    require(left.shape == right.shape, f"{label} shape mismatch")
    if np.issubdtype(left.dtype, np.number) and np.issubdtype(right.dtype, np.number):
        require(np.isfinite(left.astype(float)).all() and np.isfinite(right.astype(float)).all(), f"nonfinite {label}")
        require(np.allclose(left.astype(float), right.astype(float), rtol=0, atol=tolerance), f"{label} mismatch")
    else:
        require(np.array_equal(left.astype(str), right.astype(str)), f"{label} mismatch")


def require_schema(frame, columns, rows, key, label):
    require(list(frame.columns) == list(columns), f"{label} schema mismatch")
    require(len(frame) == rows, f"{label} row count mismatch")
    require(not frame.duplicated(list(key)).any(), f"duplicate {label} key")


def require_numeric_finite(frame, columns, label):
    values = frame[list(columns)].to_numpy(float)
    require(np.isfinite(values).all(), f"nonfinite {label}")


def safe_ratio(numerator, denominator):
    numerator = float(numerator)
    denominator = float(denominator)
    if not np.isfinite(numerator) or not np.isfinite(denominator) or denominator == 0:
        return None
    value = numerator / denominator
    return float(value) if np.isfinite(value) else None


def hypothesis_results(scenarios, cfg):
    h1_cells = []
    h2_cells = []
    h1_limit = cfg["confirmatory_hypotheses"]["H1_PROPENSITY_INSTABILITY_CONFIRMED"]["maximum_ratio"]
    h2_cfg = cfg["confirmatory_hypotheses"]["H2_HIGH_DIMENSION_INSTABILITY_CONFIRMED"]
    for endpoint, mechanism in product(ENDPOINTS, MECHANISMS):
        base = scenarios[
            scenarios["endpoint"].eq(endpoint)
            & scenarios["mechanism"].eq(mechanism)
            & scenarios["loss_representation"].eq("lowdim_11")
        ]
        full = base[base["propensity_representation"].eq("full_2266")]
        low = base[base["propensity_representation"].eq("lowdim_11")]
        require(len(full) == len(low) == 400, "confirmatory aggregation row mismatch")
        b_oracle = float(full["oracle_aipw_abs_error"].median())
        b_full = float(full["fitted_aipw_abs_error"].median())
        b_low = float(low["fitted_aipw_abs_error"].median())
        h1_ratio = safe_ratio(b_oracle, b_full)
        h2_ratio = safe_ratio(b_low, b_full)
        support_valid = bool(full["fold_class_support"].all() and low["fold_class_support"].all())
        finite_valid = bool(full["all_finite"].all() and low["all_finite"].all())
        h1_adjudicative = bool(support_valid and finite_valid and h1_ratio is not None)
        h2_adjudicative = bool(support_valid and finite_valid and h2_ratio is not None)
        clip_full = float(full["fitted_clip_fraction"].median())
        clip_low = float(low["fitted_clip_fraction"].median())
        clip_reduction = clip_full - clip_low
        h1_cells.append({
            "endpoint": endpoint,
            "mechanism": mechanism,
            "rows": 400,
            "B_oracle_aipw_fullP_lowM": b_oracle,
            "B_fitted_aipw_fullP_lowM": b_full,
            "ratio": h1_ratio,
            "threshold": float(h1_limit),
            "denominator_valid": h1_ratio is not None,
            "fold_class_support": support_valid,
            "all_finite": finite_valid,
            "adjudicative": h1_adjudicative,
            "status": "OK" if h1_adjudicative else (
                "NONADJUDICATIVE_FOLD_SUPPORT" if not support_valid else
                "NONADJUDICATIVE_NONFINITE" if not finite_valid else
                "NONADJUDICATIVE_DENOMINATOR"
            ),
            "pass": bool(h1_adjudicative and b_oracle <= h1_limit * b_full),
        })
        h2_cells.append({
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
        })
    return {
        "H1": {"pass": bool(all(row["pass"] for row in h1_cells)), "required_cells": 4, "cells": h1_cells},
        "H2": {"pass": bool(all(row["pass"] for row in h2_cells)), "required_cells": 4, "cells": h2_cells},
    }


def replay_integrity(replay_path, scenario_path, contexts, generators, tolerance):
    replay = pd.read_parquet(replay_path)
    scenarios = read_result_csv(scenario_path)
    require_schema(
        replay, REPLAY_COLUMNS, 59232,
        ["endpoint", "anchor", "mechanism", "seed", "model", "propensity_representation", "loss_representation", "source_row_index"],
        "row replay",
    )
    require_numeric_finite(replay, [
        "source_row_index", "cohort_row_index", "fold", "R", "loss", "mhat", "p_raw", "p_fitted", "p_oracle",
        *[f"audit_{name}" for name in LOW_DESCRIPTOR_NAMES], "audit_max_train_similarity",
        "augmentation_fitted", "pseudo_outcome_fitted", "influence_fitted",
        "augmentation_oracle", "pseudo_outcome_oracle", "influence_oracle",
    ], "row replay")
    require(np.isin(replay["R"].to_numpy(int), [0, 1]).all(), "nonbinary replay mask")
    require(((replay[["p_raw", "p_fitted", "p_oracle"]].to_numpy(float) > 0)
             & (replay[["p_raw", "p_fitted", "p_oracle"]].to_numpy(float) <= 1)).all(), "invalid replay probability")
    keys = ["endpoint", "anchor", "mechanism", "seed", "model", "propensity_representation", "loss_representation"]
    expected_cells = set(product(ENDPOINTS, REGULAR_ANCHORS, MECHANISMS, [41001], MODELS, REPRESENTATIONS, REPRESENTATIONS))
    require(set(replay[keys].drop_duplicates().itertuples(index=False, name=None)) == expected_cells, "row replay cell grid mismatch")
    scenario_index = scenarios[scenarios["seed"].eq(41001)].set_index(keys)
    require(len(scenario_index) == len(expected_cells) and not scenario_index.index.duplicated().any(), "seed-41001 scenario grid mismatch")
    maxima = {}

    def track(label, left, right):
        difference = np.abs(np.asarray(left, dtype=float) - np.asarray(right, dtype=float))
        maximum = float(np.max(difference)) if difference.size else 0.0
        maxima[label] = max(maxima.get(label, 0.0), maximum)
        require(np.isfinite(difference).all() and maximum <= tolerance, f"row replay formula mismatch: {label}")

    for key, part in replay.groupby(keys, sort=False):
        endpoint, anchor, mechanism, seed, model, propensity_representation, loss_representation = key
        context = contexts[endpoint]
        part = part.sort_values("source_row_index").reset_index(drop=True)
        require(len(part) == len(context["source_row_index"]), f"row replay cell length mismatch: {key}")
        vector_close(part["source_row_index"], context["source_row_index"], tolerance, "replay source rows")
        vector_close(part["cohort_row_index"], context["cohort_row_index"], tolerance, "replay cohort rows")
        vector_close(part["molecule_id"].astype(str), context["molecule_id"].astype(str), tolerance, "replay molecule IDs")
        vector_close(part["scaffold_group"].astype(str), context["groups"].astype(str), tolerance, "replay scaffold groups")
        vector_close(part["scaffold_sha256"].astype(str), context["scaffold_sha256"].astype(str), tolerance, "replay scaffold hashes")
        vector_close(part["fold"], context["fold"], tolerance, "replay folds")
        p_oracle = generators[(endpoint, anchor)][mechanism]["p"]
        r = (np.random.default_rng(int(seed)).random(len(p_oracle)) < p_oracle).astype(int)
        vector_close(part["R"], r, tolerance, "replay mask")
        vector_close(part["p_oracle"], p_oracle, tolerance, "replay oracle propensity")
        track("p_fitted_clip", part["p_fitted"], np.clip(part["p_raw"], 0.05, 0.95))
        model_index = MODELS.index(model)
        vector_close(part["loss"], context["losses"][:, model_index], tolerance, "replay loss")
        audit = part[[*[f"audit_{name}" for name in LOW_DESCRIPTOR_NAMES], "audit_max_train_similarity"]].to_numpy(float)
        vector_close(audit, context["audit_features"], tolerance, "replay audit features")
        fitted = estimator(
            r, part["loss"].to_numpy(float), part["mhat"].to_numpy(float), part["p_fitted"].to_numpy(float),
            context["groups"], audit, context["audit_names"],
        )
        oracle = estimator(
            r, part["loss"].to_numpy(float), part["mhat"].to_numpy(float), p_oracle,
            context["groups"], audit, context["audit_names"],
        )
        stored = scenario_index.loc[key]
        require(stored["R_sha256"] == binary_hash(r), "replay mask hash mismatch")
        require(stored["p_oracle_sha256"] == float_hash(p_oracle), "replay oracle-propensity hash mismatch")
        require(stored["p_raw_sha256"] == float_hash(part["p_raw"]), "replay raw-propensity hash mismatch")
        require(stored["p_fitted_sha256"] == float_hash(part["p_fitted"]), "replay fitted-propensity hash mismatch")
        require(stored["mhat_sha256"] == float_hash(part["mhat"]), "replay mhat hash mismatch")
        track("truth", stored["full_truth"], part["loss"].mean())
        track("complete_case", stored["cc"], part.loc[part["R"].eq(1), "loss"].mean())
        track("outcome_regression", stored["or_mean"], part["mhat"].mean())
        track("clip_fraction", stored["fitted_clip_fraction"], np.mean((part["p_raw"] < 0.05) | (part["p_raw"] > 0.95)))
        for prefix, result in [("fitted", fitted), ("oracle", oracle)]:
            for stored_name, result_name in [
                ("ipw", "ipw"), ("aipw", "aipw"), ("se", "se"), ("lcl", "lcl"), ("ucl", "ucl"),
                ("ess", "ess"), ("ess_fraction", "ess_fraction"), ("max_abs_smd", "max_abs_smd"),
                ("max_normalized_weight", "max_normalized_weight"),
            ]:
                track(f"{prefix}_{stored_name}", stored[f"{prefix}_{stored_name}"], result[result_name])
            require(bool(stored[f"{prefix}_covered"]) == bool(result["lcl"] <= stored["full_truth"] <= result["ucl"]), f"{prefix} coverage mismatch")
            stored_smd = json.loads(stored[f"{prefix}_smd_json"])
            require(set(stored_smd) == set(context["audit_names"]), f"{prefix} SMD feature mismatch")
            for name in context["audit_names"]:
                track(f"{prefix}_smd_{name}", stored_smd[name], result["smd"][name])
            for vector_name in ["augmentation", "pseudo", "influence"]:
                replay_name = "pseudo_outcome" if vector_name == "pseudo" else vector_name
                track(f"{prefix}_{replay_name}", part[f"{replay_name}_{prefix}"], result[vector_name])
        support = fold_support(r, context["split"])
        for column in FOLD_SUPPORT_COLUMNS:
            require(int(stored[column]) == support[column], f"replay fold support mismatch: {column}")
        require(bool(stored["fold_class_support"]) == support["fold_class_support"], "replay fold class-support mismatch")
        require(bool(stored["all_finite"]), "replay all-finite flag failed")
    require(
        replay.groupby(["endpoint", "anchor", "mechanism", "seed", "propensity_representation", "source_row_index"])[["p_raw", "p_fitted"]].nunique().max().max() == 1,
        "replay propensity not shared across model/loss cells",
    )
    require(
        replay.groupby(["endpoint", "anchor", "mechanism", "seed", "model", "loss_representation", "source_row_index"])["mhat"].nunique().max() == 1,
        "replay loss nuisance not shared across propensity cells",
    )
    require(
        replay.groupby(["endpoint", "anchor", "mechanism", "seed", "source_row_index"])[["R", "p_oracle"]].nunique().max().max() == 1,
        "replay mask/oracle propensity not shared across factorial cells",
    )
    return {"pass": True, "rows": len(replay), "cells": len(expected_cells), "tolerance": tolerance, "maximum_absolute_difference": maxima}


def compare_frames(expected, observed, columns, keys, tolerance, label):
    require(list(expected.columns) == list(columns) and list(observed.columns) == list(columns), f"{label} schema mismatch")
    require(len(expected) == len(observed), f"{label} row count mismatch")
    left = expected.copy()
    right = observed.copy()
    left["__key"] = left[list(keys)].fillna("<NULL>").astype(str).agg("\x1f".join, axis=1)
    right["__key"] = right[list(keys)].fillna("<NULL>").astype(str).agg("\x1f".join, axis=1)
    require(not left["__key"].duplicated().any() and not right["__key"].duplicated().any(), f"duplicate {label} key")
    left = left.sort_values("__key").reset_index(drop=True)
    right = right.sort_values("__key").reset_index(drop=True)
    require(np.array_equal(left["__key"], right["__key"]), f"{label} key grid mismatch")
    for column in columns:
        if pd.api.types.is_numeric_dtype(left[column]) and pd.api.types.is_numeric_dtype(right[column]):
            require(
                np.allclose(left[column].to_numpy(float), right[column].to_numpy(float), rtol=0, atol=tolerance, equal_nan=True),
                f"{label} numeric mismatch: {column}",
            )
        else:
            require(
                np.array_equal(left[column].fillna("<NULL>").astype(str), right[column].fillna("<NULL>").astype(str)),
                f"{label} value mismatch: {column}",
            )


def compare_json(expected, observed, tolerance, label="root"):
    if isinstance(expected, dict):
        require(isinstance(observed, dict) and set(expected) == set(observed), f"{label} JSON keys mismatch")
        for key in expected:
            compare_json(expected[key], observed[key], tolerance, f"{label}.{key}")
    elif isinstance(expected, list):
        require(isinstance(observed, list) and len(expected) == len(observed), f"{label} JSON list mismatch")
        for index, (left, right) in enumerate(zip(expected, observed)):
            compare_json(left, right, tolerance, f"{label}[{index}]")
    elif expected is None or isinstance(expected, (str, bool)):
        require(expected == observed, f"{label} JSON value mismatch")
    elif isinstance(expected, (int, float, np.integer, np.floating)):
        close(expected, observed, tolerance, label)
    else:
        require(expected == observed, f"{label} JSON value mismatch")


def durable_integrity(output, contexts, generators, expected_tables, expected_structure, cfg):
    tolerance = float(cfg["integrity"]["replay_tolerance"])
    paths = {name: output / name for name in expected_tables}
    tables = {
        name: (pd.read_parquet(path) if path.suffix == ".parquet" else read_result_csv(path))
        for name, path in paths.items()
    }
    specifications = {
        "cohort_audit.csv": (COHORT_COLUMNS, 18, ["dataset", "split", "endpoint"]),
        "availability_map.csv": (AVAILABILITY_COLUMNS, 8, ["endpoint", "anchor"]),
        "predictions.parquet": (PREDICTION_COLUMNS, 3702, ["endpoint", "model", "source_row_index"]),
        "nuisance_fold_assignments.parquet": (FOLD_COLUMNS, 1851, ["endpoint", "source_row_index"]),
        "generator_vectors.parquet": (GENERATOR_COLUMNS, 7404, ["endpoint", "anchor", "mechanism", "source_row_index"]),
        "semisynthetic_scenarios.csv": (SEMISYNTHETIC_COLUMNS, 6400, ["endpoint", "anchor", "mechanism", "seed", "model", "propensity_representation", "loss_representation"]),
        "actual_mask_scenarios.csv": (ACTUAL_COLUMNS, 64, ["endpoint", "anchor", "model", "propensity_representation", "loss_representation"]),
        "factorial_contrasts.csv": (CONTRAST_COLUMNS, 56, ["row_type", "endpoint", "anchor", "mechanism", "model", "loss_representation", "contrast"]),
        "bootstrap_intervals.csv": (BOOTSTRAP_COLUMNS, 12, ["hypothesis", "statistic", "endpoint", "mechanism"]),
        "pxr_linkage_audit.csv": (PXR_LINKAGE_COLUMNS, 10747, ["source_row_index"]),
        "pxr_selection_strata.csv": (PXR_STRATA_COLUMNS, 2, ["formal_hit"]),
        "pxr_overlap_summary.csv": (PXR_OVERLAP_COLUMNS, 2, ["specification"]),
    }
    for name, (columns, rows, key) in specifications.items():
        require_schema(tables[name], columns, rows, key, name)
    require_schema(
        pd.read_parquet(output / "row_replay_seed41001.parquet"), REPLAY_COLUMNS, 59232,
        ["endpoint", "anchor", "mechanism", "seed", "model", "propensity_representation", "loss_representation", "source_row_index"],
        "row_replay_seed41001.parquet",
    )
    scenario_keys = ["endpoint", "anchor", "mechanism", "seed", "model", "propensity_representation", "loss_representation"]
    require(
        set(tables["semisynthetic_scenarios.csv"][scenario_keys].itertuples(index=False, name=None))
        == set(product(ENDPOINTS, REGULAR_ANCHORS, MECHANISMS, range(41001, 41101), MODELS, REPRESENTATIONS, REPRESENTATIONS)),
        "semisynthetic Cartesian grid mismatch",
    )
    actual_keys = ["endpoint", "anchor", "model", "propensity_representation", "loss_representation"]
    require(
        set(tables["actual_mask_scenarios.csv"][actual_keys].itertuples(index=False, name=None))
        == set(product(ENDPOINTS, [*REGULAR_ANCHORS, *STRESS_ANCHORS], MODELS, REPRESENTATIONS, REPRESENTATIONS)),
        "actual-mask Cartesian grid mismatch",
    )
    for name in [
        "cohort_audit.csv", "availability_map.csv", "predictions.parquet", "nuisance_fold_assignments.parquet",
        "generator_vectors.parquet", "semisynthetic_scenarios.csv", "actual_mask_scenarios.csv",
        "bootstrap_intervals.csv", "pxr_linkage_audit.csv", "pxr_selection_strata.csv", "pxr_overlap_summary.csv",
    ]:
        numeric = tables[name].select_dtypes(include=[np.number]).columns
        require_numeric_finite(tables[name], numeric, name)
    require(tables["generator_vectors.parquet"]["all_finite"].astype(bool).all(), "generator finite flag failed")
    require(tables["semisynthetic_scenarios.csv"]["all_finite"].astype(bool).all(), "scenario finite flag failed")
    require(tables["actual_mask_scenarios.csv"]["all_finite"].astype(bool).all(), "actual finite flag failed")
    compare_frames(
        expected_tables["cohort_audit.csv"], tables["cohort_audit.csv"], COHORT_COLUMNS,
        ["dataset", "split", "endpoint"], tolerance, "cohort audit",
    )
    compare_frames(
        expected_tables["availability_map.csv"], tables["availability_map.csv"], AVAILABILITY_COLUMNS,
        ["endpoint", "anchor"], tolerance, "availability map",
    )
    compare_frames(
        expected_tables["predictions.parquet"], tables["predictions.parquet"], PREDICTION_COLUMNS,
        ["endpoint", "model", "source_row_index"], tolerance, "predictions",
    )
    compare_frames(
        expected_tables["nuisance_fold_assignments.parquet"], tables["nuisance_fold_assignments.parquet"], FOLD_COLUMNS,
        ["endpoint", "source_row_index"], tolerance, "fold assignments",
    )
    compare_frames(
        expected_tables["generator_vectors.parquet"], tables["generator_vectors.parquet"], GENERATOR_COLUMNS,
        ["endpoint", "anchor", "mechanism", "source_row_index"], tolerance, "generator vectors",
    )
    scenario = tables["semisynthetic_scenarios.csv"]
    require(
        scenario.groupby(["endpoint", "anchor", "mechanism", "seed"])[
            ["R_sha256", "p_oracle_sha256"]
        ].nunique().max().max() == 1,
        "semisynthetic mask/oracle propensity not shared across factorial cells",
    )
    require(
        scenario.groupby(["endpoint", "anchor", "mechanism", "seed", "propensity_representation"])[
            ["p_raw_sha256", "p_fitted_sha256", "fitted_ess", "fitted_clip_fraction"]
        ].nunique().max().max() == 1,
        "semisynthetic fitted propensity not shared across model/loss cells",
    )
    require(
        scenario.groupby(["endpoint", "anchor", "mechanism", "seed", "model", "loss_representation"])["mhat_sha256"].nunique().max() == 1,
        "semisynthetic loss nuisance not shared across propensity cells",
    )
    require(scenario.groupby("endpoint")[["source_row_index_sha256", "scaffold_group_vector_sha256", "fold_sha256"]].nunique().max().max() == 1, "semisynthetic endpoint vectors not shared")
    actual = tables["actual_mask_scenarios.csv"]
    require(
        actual.groupby(["endpoint", "anchor"])["R_sha256"].nunique().max() == 1,
        "actual mask not shared across factorial cells",
    )
    require(
        actual.groupby(["endpoint", "anchor", "propensity_representation"])[
            ["p_raw_sha256", "p_fitted_sha256", "fitted_ess", "fitted_clip_fraction"]
        ].nunique().max().max() == 1,
        "actual-mask propensity not shared across model/loss cells",
    )
    require(
        actual.groupby(["endpoint", "anchor", "model", "loss_representation"])["mhat_sha256"].nunique().max() == 1,
        "actual-mask loss nuisance not shared across propensity cells",
    )
    require(actual.groupby("endpoint")[["source_row_index_sha256", "scaffold_group_vector_sha256", "fold_sha256"]].nunique().max().max() == 1, "actual endpoint vectors not shared")
    expected_contrasts = factorial_contrasts(scenario)
    compare_frames(
        expected_contrasts, tables["factorial_contrasts.csv"], CONTRAST_COLUMNS,
        ["row_type", "endpoint", "anchor", "mechanism", "model", "loss_representation", "contrast"],
        tolerance, "factorial contrasts",
    )
    expected_bootstrap, _ = bootstrap_intervals(scenario, cfg)
    compare_frames(
        expected_bootstrap, tables["bootstrap_intervals.csv"], BOOTSTRAP_COLUMNS,
        ["hypothesis", "statistic", "endpoint", "mechanism"], tolerance, "bootstrap intervals",
    )
    stored_structure = json.loads((output / "structure_audit.json").read_text(encoding="utf-8"))
    compare_json(expected_structure, stored_structure, tolerance, "structure_audit")
    replay = replay_integrity(
        output / "row_replay_seed41001.parquet", output / "semisynthetic_scenarios.csv",
        contexts, generators, tolerance,
    )
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
        "replay": replay,
    }, tables


def output_hashes(output):
    return {
        str(path.relative_to(output)).replace("\\", "/"): file_hash(path)
        for path in sorted(output.iterdir())
        if path.is_file() and path.name != "run_manifest.json"
    }


def full_decision(integrity, hypotheses, admission, pxr_boundary):
    if not integrity["pass"]:
        return "STOP_IMPLEMENTATION_OR_DATA"
    if not hypotheses["H1"]["pass"] or not hypotheses["H2"]["pass"]:
        return "STOP_NO_EXTERNAL_CONFIRMATION"
    boundary_complete = bool(
        pxr_boundary["decision"] == "ABSTAIN_NONIDENTIFIABLE_FULL_LIBRARY_DOSE_RESPONSE"
        and pxr_boundary["full_library_point_estimate_allowed"] is False
    )
    if (
        admission["known_truth_admission"]["pass"]
        and admission["direct_mask_overlap"]["pass"]
        and admission["direct_mask_translation"]["pass"]
        and boundary_complete
    ):
        return "GO_WRITE_CORRECTION_BOUNDED"
    require(boundary_complete, "PXR boundary incomplete")
    return "GO_WRITE_FALSIFICATION_BOUNDARY"


def run_full(
    root, output, cfg, source_amendment, input_hashes, environment, started, started_utc,
    execution_mode, input_precedence, clearance_paths, config_paths,
):
    derivative, features, descriptor_index, cohort, structure = prepare_biogen(root, cfg, source_amendment)
    contexts, predictions = fit_predictions(derivative, features, descriptor_index, cfg)
    folds = fold_table(contexts)
    generators, availability, generator_vectors = build_generator_vectors(contexts)
    scenarios, replay = semisynthetic_grid(contexts, generators, cfg)
    actual = actual_mask_grid(contexts)
    pxr_linkage, _, _, pxr_hashes, pxr_structure = prepare_pxr(root, cfg, source_amendment)
    expected_tables = {
        "cohort_audit.csv": cohort,
        "availability_map.csv": availability,
        "predictions.parquet": predictions,
        "nuisance_fold_assignments.parquet": folds,
        "generator_vectors.parquet": generator_vectors,
    }
    write_result_csv(cohort, output / "cohort_audit.csv")
    write_json(output / "structure_audit.json", structure)
    write_result_csv(availability, output / "availability_map.csv")
    predictions.to_parquet(output / "predictions.parquet", index=False)
    folds.to_parquet(output / "nuisance_fold_assignments.parquet", index=False)
    generator_vectors.to_parquet(output / "generator_vectors.parquet", index=False)
    write_result_csv(scenarios, output / "semisynthetic_scenarios.csv")
    write_result_csv(actual, output / "actual_mask_scenarios.csv")
    replay.to_parquet(output / "row_replay_seed41001.parquet", index=False)
    write_result_csv(pxr_linkage, output / "pxr_linkage_audit.csv")
    reloaded_scenarios = read_result_csv(output / "semisynthetic_scenarios.csv")
    reloaded_actual = read_result_csv(output / "actual_mask_scenarios.csv")
    contrasts = factorial_contrasts(reloaded_scenarios)
    bootstrap, bootstrap_hash = bootstrap_intervals(reloaded_scenarios, cfg)
    admission = admission_and_translation(reloaded_scenarios, reloaded_actual, cfg)
    write_result_csv(contrasts, output / "factorial_contrasts.csv")
    write_result_csv(bootstrap, output / "bootstrap_intervals.csv")
    write_json(output / "admission_gate.json", admission)
    pxr_strata, pxr_overlap, pxr_boundary = summarize_pxr(output / "pxr_linkage_audit.csv", cfg, pxr_hashes)
    write_result_csv(pxr_strata, output / "pxr_selection_strata.csv")
    write_result_csv(pxr_overlap, output / "pxr_overlap_summary.csv")
    write_json(output / "pxr_boundary_decision.json", pxr_boundary)
    expected_tables.update({
        "semisynthetic_scenarios.csv": scenarios,
        "actual_mask_scenarios.csv": actual,
        "factorial_contrasts.csv": contrasts,
        "bootstrap_intervals.csv": bootstrap,
        "pxr_linkage_audit.csv": pxr_linkage,
        "pxr_selection_strata.csv": pxr_strata,
        "pxr_overlap_summary.csv": pxr_overlap,
    })
    integrity, durable = durable_integrity(output, contexts, generators, expected_tables, structure, cfg)
    durable_admission = admission_and_translation(
        durable["semisynthetic_scenarios.csv"], durable["actual_mask_scenarios.csv"], cfg
    )
    compare_json(
        durable_admission,
        json.loads((output / "admission_gate.json").read_text(encoding="utf-8")),
        cfg["integrity"]["replay_tolerance"],
        "admission_gate",
    )
    _, expected_overlap, expected_boundary = summarize_pxr(output / "pxr_linkage_audit.csv", cfg, pxr_hashes)
    compare_frames(
        expected_overlap, durable["pxr_overlap_summary.csv"], PXR_OVERLAP_COLUMNS, ["specification"],
        cfg["integrity"]["replay_tolerance"], "PXR overlap summary",
    )
    compare_json(
        expected_boundary,
        json.loads((output / "pxr_boundary_decision.json").read_text(encoding="utf-8")),
        cfg["integrity"]["replay_tolerance"],
        "PXR boundary decision",
    )
    hypotheses = hypothesis_results(durable["semisynthetic_scenarios.csv"], cfg)
    conjunctive = bool(hypotheses["H1"]["pass"] and hypotheses["H2"]["pass"])
    decision = full_decision(integrity, hypotheses, durable_admission, expected_boundary)
    summary = {
        "study_id": "confirm_001",
        "pilot_decision_unchanged": "STOP_PILOT_NEGATIVE",
        "forensic_role": "development_only_no_confirmatory_aggregation",
        "integrity": integrity,
        "H1": hypotheses["H1"],
        "H2": hypotheses["H2"],
        "conjunctive_confirmation": conjunctive,
        "admission": durable_admission["known_truth_admission"],
        "direct_overlap": durable_admission["direct_mask_overlap"],
        "direct_translation": durable_admission["direct_mask_translation"],
        "pxr_boundary": expected_boundary,
        "decision": decision,
        "anti_claims": cfg["anti_claims"],
    }
    write_json(output / "summary.json", summary)
    verification = (
        "# confirm_001 runner verification\n\n"
        f"- Internal durable-grid and seed-41001 formula replay: `{'PASS' if integrity['pass'] else 'FAIL'}`.\n"
        f"- H1 all four cells: `{'PASS' if hypotheses['H1']['pass'] else 'FAIL'}`.\n"
        f"- H2 all four cells: `{'PASS' if hypotheses['H2']['pass'] else 'FAIL'}`.\n"
        f"- Frozen decision: `{decision}`.\n\n"
        "Scope: this file records the runner's reload-and-replay check. It is not the independent full verifier, "
        "does not alter `STOP_PILOT_NEGATIVE`, and does not identify full-library PXR dose-response performance.\n"
    )
    (output / "verification.md").write_text(verification, encoding="utf-8")
    required_outputs = {
        "cohort_audit.csv", "structure_audit.json", "availability_map.csv", "predictions.parquet",
        "nuisance_fold_assignments.parquet", "generator_vectors.parquet", "semisynthetic_scenarios.csv",
        "actual_mask_scenarios.csv", "row_replay_seed41001.parquet", "factorial_contrasts.csv",
        "bootstrap_intervals.csv", "admission_gate.json", "pxr_linkage_audit.csv",
        "pxr_selection_strata.csv", "pxr_overlap_summary.csv", "pxr_boundary_decision.json",
        "summary.json", "verification.md",
    }
    require(set(output_hashes(output)) == required_outputs, "formal output inventory mismatch before manifest")
    finished_utc = datetime.now(timezone.utc).isoformat()
    verifier_sha = input_hashes["src/verify_confirmatory.py"]
    manifest = {
        "study_id": "confirm_001",
        "mode": execution_mode,
        "status": "COMPLETED_PENDING_INDEPENDENT_FULL_VERIFICATION",
        "started_utc": started_utc,
        "finished_utc": finished_utc,
        "duration_seconds": float(time.time() - started),
        "command": [sys.executable, *sys.argv],
        "input_precedence": input_precedence,
        "input_hashes": input_hashes,
        "runner_sha256": file_hash(Path(__file__)),
        "verifier_sha256": verifier_sha,
        "clearance_hashes": {path: input_hashes[path] for path in clearance_paths},
        "environment": environment,
        "config_paths": config_paths,
        "seeds": {
            "mask_seeds": {"first": 41001, "last": 41100, "count": 100},
            "row_replay_seed": 41001,
            "bootstrap_seed": 52001,
            "bootstrap_resamples": 10000,
            "bootstrap_seed_index_matrix_sha256": bootstrap_hash,
        },
        "feature_contract": {
            "descriptor_registry_count": 217,
            "descriptor_names": structure["descriptor_names"],
            "descriptor_names_sha256": structure["descriptor_names_sha256"],
            "nuisance_representations": {
                endpoint: {
                    name: {"shape": list(contexts[endpoint]["representations"][name].shape), "dtype": str(contexts[endpoint]["representations"][name].dtype)}
                    for name in REPRESENTATIONS
                }
                for endpoint in ENDPOINTS
            },
            "lowdim_names": [*LOW_DESCRIPTOR_NAMES, "max_train_similarity"],
            "mcar_generator_neutral_storage": {"generator_score": 0.0, "intercept_shift": 0.0},
        },
        "grid_rows": {
            "cohort_audit.csv": 18,
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
        },
        "vector_hashes": {
            "biogen": {
                endpoint: {
                    **contexts[endpoint]["hashes"],
                    "generator_p_oracle_sha256": {
                        f"{anchor}_{mechanism}": float_hash(generators[(endpoint, anchor)][mechanism]["p"])
                        for anchor in REGULAR_ANCHORS for mechanism in MECHANISMS
                    },
                }
                for endpoint in ENDPOINTS
            },
            "pxr": pxr_hashes,
            "bootstrap_seed_index_matrix_sha256": bootstrap_hash,
        },
        "pxr_structure_audit": pxr_structure,
        "integrity": integrity,
        "decision": decision,
        "output_hashes_before_manifest": output_hashes(output),
    }
    write_json(output / "run_manifest.json", manifest)


def run_smoke(
    output, input_hashes, environment, started, started_utc,
    execution_mode, clearance_paths, recovery,
):
    rng = np.random.default_rng(20260911)
    low = rng.normal(size=(72, 11)).astype(np.float64)
    full = np.column_stack([low, rng.normal(size=(72, 2255))]).astype(np.float64)
    groups = np.repeat([f"SMOKE_GROUP_{index}" for index in range(24)], 3)
    split, fold = fold_contract(groups)
    r = np.tile([0, 1], 36).astype(int)
    losses = np.abs(rng.normal(size=(72, 2)))
    representations = {"lowdim_11": low, "full_2266": full}
    prepared = {name: scaled_folds(value, split) for name, value in representations.items()}
    propensity = {name: fit_propensity(prepared[name], r) for name in REPRESENTATIONS}
    mhat = {name: fit_loss(prepared[name], r, losses) for name in REPRESENTATIONS}
    estimates = []
    for propensity_representation, loss_representation, model_index in product(REPRESENTATIONS, REPRESENTATIONS, range(2)):
        result = estimator(
            r, losses[:, model_index], mhat[loss_representation][:, model_index], propensity[propensity_representation]["p"],
            groups, low, [f"x{index}" for index in range(11)],
        )
        estimates.append([result[name] for name in ["ipw", "aipw", "se", "ess", "ess_fraction"]])
    numerator = 408.40724836686263
    denominator = 0.004053686965995372
    fixture = pd.DataFrame({
        "fixture_id": ["synthetic_large_ratio"],
        "numerator": [numerator],
        "denominator": [denominator],
        "ratio": [numerator / denominator],
    })
    fixture_path = output / "roundtrip_fixture.csv"
    write_result_csv(fixture, fixture_path)
    reloaded_fixture = read_result_csv(fixture_path)
    numeric_columns = ["numerator", "denominator", "ratio"]
    component_bit_exact = bool(np.array_equal(
        fixture[numeric_columns].to_numpy(dtype="<f8").view("<u8"),
        reloaded_fixture[numeric_columns].to_numpy(dtype="<f8").view("<u8"),
    ))
    stored_ratio = float(reloaded_fixture.loc[0, "ratio"])
    recomputed_ratio = float(reloaded_fixture.loc[0, "numerator"] / reloaded_fixture.loc[0, "denominator"])
    ratio_difference = abs(stored_ratio - recomputed_ratio)
    checks = {
        "folds_complete": bool(np.array_equal(np.unique(fold), np.arange(3))),
        "fold_class_support": bool(fold_support(r, split)["fold_class_support"]),
        "raw_propensity_finite": bool(all(np.isfinite(value["p_raw"]).all() for value in propensity.values())),
        "clipped_propensity_finite": bool(all(np.isfinite(value["p"]).all() for value in propensity.values())),
        "loss_nuisance_finite": bool(all(np.isfinite(value).all() for value in mhat.values())),
        "estimates_finite": bool(np.isfinite(np.asarray(estimates, dtype=float)).all()),
        "nuisance_dtype_float64": bool(all(value.dtype == np.float64 for value in representations.values())),
        "csv_roundtrip_component_bit_exact": component_bit_exact,
        "csv_roundtrip_ratio_within_absolute_tolerance": bool(ratio_difference <= 1e-12),
        "csv_roundtrip_large_ratio": bool(abs(stored_ratio) >= recovery["csv_roundtrip_contract"]["recovery_smoke_minimum_ratio_magnitude"]),
    }
    require(all(checks.values()), "synthetic smoke failed")
    smoke = {
        "study_id": "confirm_001",
        "mode": execution_mode,
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
        "fold_sha256": int_hash(fold),
        "input_hashes": input_hashes,
        "runner_sha256": file_hash(Path(__file__)),
        "verifier_sha256": input_hashes["src/verify_confirmatory.py"],
    }
    write_json(output / "smoke.json", smoke)
    manifest = {
        "study_id": "confirm_001",
        "mode": execution_mode,
        "status": "PASS_NON_ADJUDICATIVE",
        "started_utc": started_utc,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "duration_seconds": float(time.time() - started),
        "command": [sys.executable, *sys.argv],
        "input_hashes": input_hashes,
        "runner_sha256": file_hash(Path(__file__)),
        "verifier_sha256": input_hashes["src/verify_confirmatory.py"],
        "clearance_hashes": {path: input_hashes[path] for path in clearance_paths},
        "environment": environment,
        "synthetic_only": True,
        "output_hashes_before_manifest": output_hashes(output),
    }
    write_json(output / "run_manifest.json", manifest)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["recovery-smoke", "recovery-full"], required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--decision-004-sha")
    parser.add_argument("--decision-005-sha")
    parser.add_argument("--decision-006-sha", required=True)
    parser.add_argument("--decision-007-sha")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    cfg, source_amendment, input_hashes, environment, recovery = verify_inputs(
        root, args.mode, args.decision_004_sha, args.decision_005_sha,
        args.decision_006_sha, args.decision_007_sha,
    )
    recovery_smoke = (root / recovery["recovery"]["smoke_output_path"]).resolve()
    recovery_full = (root / recovery["recovery"]["full_output_path"]).resolve()
    if args.mode == "recovery-full":
        require(args.decision_007_sha is not None, "recovery-full requires Decision 007")
        require(args.output is None or args.output.resolve() == recovery_full, "recovery-full output path mismatch")
        require(recovery_smoke.is_dir(), "authorized recovery smoke is missing")
        require((root / RECOVERY_SMOKE_VERIFICATION_PATH).is_file(), "independent recovery smoke verification is missing")
        require(not recovery_full.exists(), "recovery full result directory must be absent before launch")
        output = recovery_full
    else:
        require(args.decision_007_sha is None, "Decision 007 is not accepted for recovery-smoke")
        require(args.output is not None, "recovery-smoke requires an explicit --output")
        output = args.output.resolve()
        require(output == recovery_smoke, "recovery-smoke output path mismatch")
        require(not (root / DECISION_007_PATH).exists(), "Decision 007 must be absent before recovery smoke")
        require(not (root / RECOVERY_SMOKE_VERIFICATION_PATH).exists(), "recovery smoke verification record must be absent before launch")
        require(not recovery_smoke.exists(), "recovery smoke output directory must be absent before launch")
        require(not recovery_full.exists(), "recovery full output exists before smoke")
    output.mkdir(parents=True, exist_ok=False)
    started = time.time()
    started_utc = datetime.now(timezone.utc).isoformat()
    print(f"START mode={args.mode} runner={file_hash(Path(__file__))}", flush=True)
    if args.mode == "recovery-smoke":
        run_smoke(
            output, input_hashes, environment, started, started_utc,
            args.mode, [DECISION_006_PATH], recovery,
        )
    else:
        run_full(
            root, output, cfg, source_amendment, input_hashes, environment, started, started_utc,
            args.mode,
            recovery["formal_input_precedence"],
            [DECISION_006_PATH, DECISION_007_PATH],
            {
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
        )
    print(f"DONE mode={args.mode} seconds={time.time() - started:.3f}", flush=True)


if __name__ == "__main__":
    main()
