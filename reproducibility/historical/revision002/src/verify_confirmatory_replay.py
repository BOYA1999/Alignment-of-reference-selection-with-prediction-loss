import argparse
import hashlib
import importlib.util
import json
import re
import sys
from itertools import product
from pathlib import Path

import numpy as np


FAILURE_SHA256 = "3a6ac0a67988d3037831f9f21a084b814f78e701ec464c5e224b09018e74c8a1"
PROTOCOL_AMENDMENT_003_SHA256 = "458d0217da500b89a89c28018e07f2df18b3d741f1fb13e068451cef0a447adb"
CONFIG_AMENDMENT_003_SHA256 = "690875bca0afcd95a44946c57aa552c3287c1c80efb690cb412318565999ab60"
AMENDMENT_003_AUDIT_SHA256 = "61a67e578d57ce790182646c61ba8da18ce8f49d7fb9f6853f4033b21156fe75"
RUNNER_SHA256 = "650a8b8f6b9bf65db69efc7ba60340aafa99c3c4b8b5bba957a8e6c6fa800c5c"
ORIGINAL_VERIFIER_SHA256 = "740e31590ec40f40f7783e59889d265c3d5d20296901c3968d4b593986238715"
DECISION_006_SHA256 = "96cd052f18c9442ebe9b3b3524da7b88539c0896f6541dbc675cf9cb094cee0b"
DECISION_007_SHA256 = "faea63b749016aaaf570888ca5cbc99fda856088c53ac02f489fca3df7f45fe7"
FULL_MANIFEST_SHA256 = "8da9b9c1d0425fff118d983a8491aceb57ae8b6085300eacf18cd8e0d7e6eac4"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def field(text, label):
    pattern = re.compile(rf"^\s*(?:[-*]\s*)?(?:\*\*)?{re.escape(label)}(?:\*\*)?\s*:\s*(.+?)\s*$")
    values = [match.group(1).strip().strip("`") for line in text.splitlines() if (match := pattern.match(line))]
    require(len(values) == 1 and values[0], f"missing or duplicate decision field: {label}")
    return values[0]


def load_original_verifier(root):
    path = root / "src/verify_confirmatory.py"
    require(sha256(path) == ORIGINAL_VERIFIER_SHA256, "original verifier hash mismatch")
    spec = importlib.util.spec_from_file_location("confirm001_original_verifier", path)
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    return verifier


def verify_authorization(root, decision_008_sha):
    paths = {
        "artifacts/confirmatory/confirm_001/independent_verification_failure_001.md": FAILURE_SHA256,
        "artifacts/confirmatory/confirm_001/protocol_amendment_003.md": PROTOCOL_AMENDMENT_003_SHA256,
        "artifacts/confirmatory/confirm_001/config_amendment_003.json": CONFIG_AMENDMENT_003_SHA256,
        "artifacts/confirmatory/confirm_001/amendment_003_independent_audit.md": AMENDMENT_003_AUDIT_SHA256,
        "src/run_confirmatory.py": RUNNER_SHA256,
        "src/verify_confirmatory.py": ORIGINAL_VERIFIER_SHA256,
        "artifacts/decisions/006_confirm_001_clear_to_recovery_smoke.md": DECISION_006_SHA256,
        "artifacts/decisions/007_confirm_001_clear_to_recovery_full.md": DECISION_007_SHA256,
        "artifacts/confirmatory/confirm_001/full_recovery_001/run_manifest.json": FULL_MANIFEST_SHA256,
    }
    for relative, expected in paths.items():
        path = root / relative
        require(path.is_file() and sha256(path) == expected, f"authorization input mismatch: {relative}")

    amendment_audit = (root / "artifacts/confirmatory/confirm_001/amendment_003_independent_audit.md").read_text(encoding="utf-8")
    require(field(amendment_audit, "Verdict") == "CLEAR_TO_IMPLEMENT_VERIFIER_REPLAY", "Amendment 003 audit did not clear implementation")
    require(field(amendment_audit, "Failure record SHA-256") == FAILURE_SHA256, "Amendment 003 audit failure-record mismatch")
    require(field(amendment_audit, "Protocol Amendment 003 SHA-256") == PROTOCOL_AMENDMENT_003_SHA256, "Amendment 003 audit protocol mismatch")
    require(field(amendment_audit, "Config Amendment 003 SHA-256") == CONFIG_AMENDMENT_003_SHA256, "Amendment 003 audit config mismatch")

    config_path = root / "artifacts/confirmatory/confirm_001/config_amendment_003.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    require(config["status"] == "PENDING_INDEPENDENT_VERIFICATION_RECOVERY_AUDIT", "Amendment 003 config state mismatch")
    require(config["runner_execution_authorized"] is False and config["scientific_interpretation_allowed_before_pass"] is False, "Amendment 003 boundary mismatch")

    decision_path = root / config["fixed_paths"]["decision_008"]
    require(re.fullmatch(r"[0-9a-f]{64}", decision_008_sha or "") is not None, "trusted Decision 008 SHA-256 is required")
    require(decision_path.is_file() and sha256(decision_path) == decision_008_sha, "Decision 008 hash mismatch")
    text = decision_path.read_text(encoding="utf-8")
    own_sha = sha256(Path(__file__))
    static_path = root / config["fixed_paths"]["static_audit"]
    require(static_path.is_file(), "verification-replay static audit is missing")
    static_sha = sha256(static_path)
    required = {
        "Decision": "CLEAR_TO_VERIFICATION_REPLAY",
        "Mode": "verification-replay",
        "Supplemental verifier SHA-256": own_sha,
        "Original runner SHA-256": RUNNER_SHA256,
        "Original verifier SHA-256": ORIGINAL_VERIFIER_SHA256,
        "Decision 006 SHA-256": DECISION_006_SHA256,
        "Consumed Decision 007 SHA-256": DECISION_007_SHA256,
        "Failure record SHA-256": FAILURE_SHA256,
        "Protocol Amendment 003 SHA-256": PROTOCOL_AMENDMENT_003_SHA256,
        "Config Amendment 003 SHA-256": CONFIG_AMENDMENT_003_SHA256,
        "Amendment 003 audit SHA-256": AMENDMENT_003_AUDIT_SHA256,
        "Verification-replay static audit SHA-256": static_sha,
        "Verification-replay static audit verdict": "CLEAR_TO_VERIFICATION_REPLAY",
        "Full output path": config["immutable_full_output"]["path"],
        "Full output file count": str(config["immutable_full_output"]["file_count"]),
        "Full run manifest SHA-256": FULL_MANIFEST_SHA256,
        "Verification output path": config["fixed_paths"]["verification_output"],
    }
    for label, expected in required.items():
        require(field(text, label) == expected, f"Decision 008 field mismatch: {label}")
    static_text = static_path.read_text(encoding="utf-8")
    require(field(static_text, "Verdict") == "CLEAR_TO_VERIFICATION_REPLAY", "static audit did not clear replay")
    require(field(static_text, "Supplemental verifier SHA-256") == own_sha, "static audit supplemental-verifier mismatch")
    require(field(static_text, "Original verifier SHA-256") == ORIGINAL_VERIFIER_SHA256, "static audit original-verifier mismatch")

    full = root / config["immutable_full_output"]["path"]
    entries = list(full.iterdir()) if full.is_dir() else []
    expected_files = set(config["immutable_full_output"]["expected_files"])
    require({path.name for path in entries} == expected_files and all(path.is_file() for path in entries), "immutable full-output inventory mismatch")
    require(sha256(full / "run_manifest.json") == FULL_MANIFEST_SHA256, "immutable full manifest changed")
    verification_output = root / config["fixed_paths"]["verification_output"]
    require(not verification_output.exists(), "independent verification output must be absent before replay")
    return config, full, verification_output, own_sha


def verify_stored_generator_durability(verifier, full):
    path = full / "generator_vectors.parquet"
    frame = verifier.read_table(path)
    verifier.require_schema(frame, "generator_vectors.parquet")
    stored_hash_cells = 0
    formula_bit_exact_cells = 0
    for endpoint, anchor, mechanism in product(verifier.ENDPOINTS, verifier.REGULAR_ANCHORS, verifier.MECHANISMS):
        part = frame[
            frame["endpoint"].eq(endpoint)
            & frame["anchor"].eq(anchor)
            & frame["mechanism"].eq(mechanism)
        ].sort_values("source_row_index")
        hashes = part["p_oracle_sha256"].astype(str).unique()
        probability = part["p_oracle"].to_numpy(dtype="<f8")
        require(len(hashes) == 1 and verifier.vector_hash(probability, "float") == hashes[0], f"stored generator probability rehash mismatch: {endpoint}/{anchor}/{mechanism}")
        stored_hash_cells += 1
        if mechanism == "mar":
            shifts = part["intercept_shift"].to_numpy(dtype="<f8")
            require(np.all(shifts.view("<u8") == shifts[0:1].view("<u8")), f"stored MAR shift is not constant: {endpoint}/{anchor}")
            replay = np.clip(verifier.expit(part["generator_score"].to_numpy(dtype="<f8") + shifts[0]), 0.10, 0.90)
            require(np.array_equal(replay.view("<u8"), probability.view("<u8")), f"stored MAR score-shift formula mismatch: {endpoint}/{anchor}")
            formula_bit_exact_cells += 1
    return {
        "generator_vectors_sha256": sha256(path),
        "stored_probability_hash_cells_bit_exact": stored_hash_cells,
        "stored_mar_score_shift_formula_cells_bit_exact": formula_bit_exact_cells,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--decision-006-sha", required=True)
    parser.add_argument("--decision-007-sha", required=True)
    parser.add_argument("--decision-008-sha", required=True)
    args = parser.parse_args()
    require(args.decision_006_sha == DECISION_006_SHA256, "Decision 006 argument mismatch")
    require(args.decision_007_sha == DECISION_007_SHA256, "Decision 007 argument mismatch")
    root = Path(__file__).resolve().parents[1]
    config, full, verification_output, own_sha = verify_authorization(root, args.decision_008_sha)
    verifier = load_original_verifier(root)
    layout_evidence = {}
    original_loader = verifier.load_biogen_contexts

    def f_order_loader(*loader_args, **loader_kwargs):
        data, contexts, structure = original_loader(*loader_args, **loader_kwargs)
        for endpoint, context in contexts.items():
            original = context["lowdim"]
            require(original.dtype == np.float64 and original.flags.c_contiguous and not original.flags.f_contiguous, f"unexpected preconversion layout: {endpoint}")
            require(list(original.strides) == config["root_cause"]["verifier_strides_bytes"][endpoint], f"unexpected preconversion strides: {endpoint}")
            converted = np.asfortranarray(original, dtype=np.float64)
            require(converted.shape == original.shape and converted.dtype == np.float64 and converted.flags.f_contiguous, f"invalid F-order conversion: {endpoint}")
            require(list(converted.strides) == config["root_cause"]["runner_strides_bytes"][endpoint], f"unexpected postconversion strides: {endpoint}")
            require(np.array_equal(converted.view("<u8"), original.view("<u8")), f"F-order conversion changed values: {endpoint}")
            require(np.array_equal(context["representations"]["lowdim_11"].view("<u8"), original.view("<u8")), f"lowdim context aliases differ before conversion: {endpoint}")
            context["lowdim"] = converted
            context["representations"]["lowdim_11"] = converted
            layout_evidence[endpoint] = {
                "shape": list(converted.shape),
                "dtype": str(converted.dtype),
                "preconversion_c_contiguous": True,
                "preconversion_strides_bytes": list(original.strides),
                "postconversion_f_contiguous": True,
                "postconversion_strides_bytes": list(converted.strides),
                "elementwise_bit_identity": True,
            }
        return data, contexts, structure

    verifier.load_biogen_contexts = f_order_loader
    generator_storage = verify_stored_generator_durability(verifier, full)
    base_config_path = root / "artifacts/confirmatory/confirm_001/config.json"
    protocol_path = root / "artifacts/confirmatory/confirm_001/protocol.md"
    require(verifier.sha256(protocol_path) == verifier.EXPECTED_PROTOCOL_SHA256, "base protocol hash mismatch")
    require(verifier.sha256(base_config_path) == verifier.EXPECTED_CONFIG_SHA256, "base config hash mismatch")
    base_config = json.loads(base_config_path.read_text(encoding="utf-8"))
    corrections, config_overlay = verifier.amended_hashes(root)
    recovery = verifier.load_recovery_amendment(root)
    environment = verifier.check_environment(base_config)
    inputs, runner_sha, original_verifier_sha = verifier.expected_input_hashes(
        root, base_config, base_config_path, config_overlay, recovery, environment, "recovery-full",
        None, None, args.decision_006_sha, args.decision_007_sha,
    )
    require(runner_sha == RUNNER_SHA256 and original_verifier_sha == ORIGINAL_VERIFIER_SHA256, "original code identity mismatch")
    decision, nuisance_difference = verifier.verify_full(
        root, full, base_config, corrections, config_overlay, recovery, inputs, runner_sha,
        original_verifier_sha, environment, args.decision_006_sha, args.decision_007_sha,
    )
    require(set(layout_evidence) == set(verifier.ENDPOINTS), "layout conversion evidence incomplete")
    result = {
        "status": "PASS",
        "mode": "verification-replay",
        "acceptance_state": "PASS_PENDING_POST_VERIFICATION_AUDIT",
        "scientific_decision": decision,
        "full_output_path": config["immutable_full_output"]["path"],
        "full_run_manifest_sha256": FULL_MANIFEST_SHA256,
        "decision_006_sha256": args.decision_006_sha,
        "decision_007_sha256": args.decision_007_sha,
        "decision_008_sha256": args.decision_008_sha,
        "original_runner_sha256": runner_sha,
        "original_verifier_sha256": original_verifier_sha,
        "supplemental_verifier_sha256": own_sha,
        "failure_record_sha256": FAILURE_SHA256,
        "protocol_amendment_003_sha256": PROTOCOL_AMENDMENT_003_SHA256,
        "config_amendment_003_sha256": CONFIG_AMENDMENT_003_SHA256,
        "amendment_003_audit_sha256": AMENDMENT_003_AUDIT_SHA256,
        "layout_evidence": layout_evidence,
        "generator_storage_checks": generator_storage,
        "seed41001_independent_nuisance_refit_maximum_absolute_difference": nuisance_difference,
        "scope": {
            "independent_source_feature_prediction_generator_and_seed41001_formula_reconstruction": True,
            "independent_seed41001_crossfit_propensity_and_ridge_refit": True,
            "independent_remaining_99_seed_and_direct_or_pxr_nuisance_refit": False,
            "remaining_formal_nuisance_vectors_checked_by_hash_and_grid_invariants": True,
            "full_output_write_performed": False,
        },
    }
    payload = json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n"
    with verification_output.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(payload)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(json.dumps({"status": "FAIL", "error": str(error)}, ensure_ascii=False), file=sys.stderr)
        raise
