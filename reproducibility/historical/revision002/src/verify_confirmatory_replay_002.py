import argparse
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np


RUNNER_SHA256 = "650a8b8f6b9bf65db69efc7ba60340aafa99c3c4b8b5bba957a8e6c6fa800c5c"
ORIGINAL_VERIFIER_SHA256 = "740e31590ec40f40f7783e59889d265c3d5d20296901c3968d4b593986238715"
SUPPLEMENTAL_VERIFIER_001_SHA256 = "4e95a20d5a3c74d5c96122b80fc37cf5d44892c4438e846b781899fbf2212a57"
FAILURE_001_SHA256 = "3a6ac0a67988d3037831f9f21a084b814f78e701ec464c5e224b09018e74c8a1"
FAILURE_002_SHA256 = "455c8e3928099fb31452801b5cec1c12d3ab716eaba1a03dd4ca8ee59a94e634"
PROTOCOL_AMENDMENT_002_SHA256 = "d34c1053da408a3e349c2ecdc20d1e49f8aa5cda9dd491746146d621c5988d9f"
CONFIG_AMENDMENT_002_SHA256 = "144af9c376929f24699452051abed2ddb5243ab8cdf47c8f4024e219ac278984"
AMENDMENT_002_AUDIT_SHA256 = "7a266a7614729ccf5bb2ff13890eabb6d2e6e32cc0ebc66ea21a34c8f1a821c3"
PROTOCOL_AMENDMENT_003_SHA256 = "458d0217da500b89a89c28018e07f2df18b3d741f1fb13e068451cef0a447adb"
CONFIG_AMENDMENT_003_SHA256 = "690875bca0afcd95a44946c57aa552c3287c1c80efb690cb412318565999ab60"
AMENDMENT_003_AUDIT_SHA256 = "61a67e578d57ce790182646c61ba8da18ce8f49d7fb9f6853f4033b21156fe75"
PROTOCOL_AMENDMENT_004_SHA256 = "00482a584e320e7b8188fe29a6db9f1f5bf1349874062c473e1d407273e8f3e3"
CONFIG_AMENDMENT_004_SHA256 = "5be83a824a04f26729e1d51fe1aa547ebf750430ad1469a9dee7270b364a65df"
AMENDMENT_004_AUDIT_SHA256 = "104484171d0e567463bf8584c49fcf84ea2bd8d36e0f9bf29a4b413f0e040373"
RECOVERY_STATIC_AUDIT_SHA256 = "a36f7202eb6fae23229cfa790af8166d3b607e458ece3280b42d6eed26a430cf"
VERIFICATION_REPLAY_001_STATIC_AUDIT_SHA256 = "df4c666f97390258e164d7c65cda527d388d7e82f7fefbfc9ba44cd6010ffe35"
DECISION_006_SHA256 = "96cd052f18c9442ebe9b3b3524da7b88539c0896f6541dbc675cf9cb094cee0b"
DECISION_007_SHA256 = "faea63b749016aaaf570888ca5cbc99fda856088c53ac02f489fca3df7f45fe7"
DECISION_008_SHA256 = "3cc6426bdebeb739d31d8d0cee1092d7de7ac8c4625b6b2004011390fec43216"
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


def load_module(path, expected_sha, name):
    require(Path(path).is_file() and sha256(path) == expected_sha, f"module identity mismatch: {Path(path).name}")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_authorization(root, decision_009_sha):
    paths = {
        "src/run_confirmatory.py": RUNNER_SHA256,
        "src/verify_confirmatory.py": ORIGINAL_VERIFIER_SHA256,
        "src/verify_confirmatory_replay.py": SUPPLEMENTAL_VERIFIER_001_SHA256,
        "artifacts/confirmatory/confirm_001/independent_verification_failure_001.md": FAILURE_001_SHA256,
        "artifacts/confirmatory/confirm_001/independent_verification_failure_002.md": FAILURE_002_SHA256,
        "artifacts/confirmatory/confirm_001/protocol_amendment_002.md": PROTOCOL_AMENDMENT_002_SHA256,
        "artifacts/confirmatory/confirm_001/config_amendment_002.json": CONFIG_AMENDMENT_002_SHA256,
        "artifacts/confirmatory/confirm_001/amendment_002_independent_audit.md": AMENDMENT_002_AUDIT_SHA256,
        "artifacts/confirmatory/confirm_001/protocol_amendment_003.md": PROTOCOL_AMENDMENT_003_SHA256,
        "artifacts/confirmatory/confirm_001/config_amendment_003.json": CONFIG_AMENDMENT_003_SHA256,
        "artifacts/confirmatory/confirm_001/amendment_003_independent_audit.md": AMENDMENT_003_AUDIT_SHA256,
        "artifacts/confirmatory/confirm_001/protocol_amendment_004.md": PROTOCOL_AMENDMENT_004_SHA256,
        "artifacts/confirmatory/confirm_001/config_amendment_004.json": CONFIG_AMENDMENT_004_SHA256,
        "artifacts/confirmatory/confirm_001/amendment_004_independent_audit.md": AMENDMENT_004_AUDIT_SHA256,
        "artifacts/confirmatory/confirm_001/recovery_001_static_code_audit.md": RECOVERY_STATIC_AUDIT_SHA256,
        "artifacts/confirmatory/confirm_001/verification_replay_001_static_audit.md": VERIFICATION_REPLAY_001_STATIC_AUDIT_SHA256,
        "artifacts/decisions/006_confirm_001_clear_to_recovery_smoke.md": DECISION_006_SHA256,
        "artifacts/decisions/007_confirm_001_clear_to_recovery_full.md": DECISION_007_SHA256,
        "artifacts/decisions/008_confirm_001_clear_to_verification_replay.md": DECISION_008_SHA256,
        "artifacts/confirmatory/confirm_001/full_recovery_001/run_manifest.json": FULL_MANIFEST_SHA256,
    }
    for relative, expected in paths.items():
        path = root / relative
        require(path.is_file() and sha256(path) == expected, f"authorization input mismatch: {relative}")

    prior = load_module(root / "src/verify_confirmatory_replay.py", SUPPLEMENTAL_VERIFIER_001_SHA256, "confirm001_replay_001")
    config_003, prior_full, failed_output, _ = prior.verify_authorization(root, DECISION_008_SHA256)
    require(not failed_output.exists(), "failed verification output must remain absent")

    audit_path = root / "artifacts/confirmatory/confirm_001/amendment_004_independent_audit.md"
    audit_text = audit_path.read_text(encoding="utf-8")
    audit_fields = {
        "Verdict": "CLEAR_TO_IMPLEMENT_VERIFICATION_REPLAY_002",
        "Failure record 002 SHA-256": FAILURE_002_SHA256,
        "Protocol Amendment 004 SHA-256": PROTOCOL_AMENDMENT_004_SHA256,
        "Config Amendment 004 SHA-256": CONFIG_AMENDMENT_004_SHA256,
        "Consumed Decision 008 SHA-256": DECISION_008_SHA256,
        "Supplemental verifier 001 SHA-256": SUPPLEMENTAL_VERIFIER_001_SHA256,
        "Original verifier SHA-256": ORIGINAL_VERIFIER_SHA256,
        "Original runner SHA-256": RUNNER_SHA256,
        "Full run manifest SHA-256": FULL_MANIFEST_SHA256,
    }
    for label, expected in audit_fields.items():
        require(field(audit_text, label) == expected, f"Amendment 004 audit field mismatch: {label}")

    config_path = root / "artifacts/confirmatory/confirm_001/config_amendment_004.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    require(config["status"] == "PENDING_INDEPENDENT_VERIFICATION_RECOVERY_AUDIT", "Amendment 004 config state mismatch")
    require(config["runner_execution_authorized"] is False, "runner execution boundary mismatch")
    require(config["scientific_interpretation_allowed_before_pass_and_post_audit"] is False, "scientific interpretation boundary mismatch")
    require(config["root_cause"]["source_and_linkage_values_bit_exact"] is True, "PXR root-cause identity mismatch")
    require(config["root_cause"]["pxr_trigger_or_boundary_changes"] == 0, "PXR scientific-boundary mismatch")
    require(config["preserved_checks"]["original_smd_json_exact_comparison"] is True, "exact SMD JSON gate was not preserved")
    require(config["preserved_checks"]["tolerance_relaxed"] is False, "tolerance relaxation is prohibited")
    require(config["preserved_checks"]["compare_frames_patch_allowed"] is False, "compare_frames patch is prohibited")
    require(config["preserved_checks"]["smd_map_patch_allowed"] is False, "smd_map patch is prohibited")
    require(config["preserved_checks"]["verify_pxr_linkage_patch_allowed"] is False, "verify_pxr_linkage patch is prohibited")

    decision_path = root / config["fixed_paths"]["decision_009"]
    require(re.fullmatch(r"[0-9a-f]{64}", decision_009_sha or "") is not None, "trusted Decision 009 SHA-256 is required")
    require(decision_path.is_file() and sha256(decision_path) == decision_009_sha, "Decision 009 hash mismatch")
    decision_text = decision_path.read_text(encoding="utf-8")
    own_sha = sha256(Path(__file__))
    static_path = root / config["fixed_paths"]["static_audit"]
    require(static_path.is_file(), "verification-replay-002 static audit is missing")
    static_sha = sha256(static_path)
    required_decision_fields = {
        "Decision": "CLEAR_TO_VERIFICATION_REPLAY_002",
        "Mode": "verification-replay-002",
        "Supplemental verifier 002 SHA-256": own_sha,
        "Supplemental verifier 001 SHA-256": SUPPLEMENTAL_VERIFIER_001_SHA256,
        "Original runner SHA-256": RUNNER_SHA256,
        "Original verifier SHA-256": ORIGINAL_VERIFIER_SHA256,
        "Decision 006 SHA-256": DECISION_006_SHA256,
        "Consumed Decision 007 SHA-256": DECISION_007_SHA256,
        "Consumed Decision 008 SHA-256": DECISION_008_SHA256,
        "Failure record 001 SHA-256": FAILURE_001_SHA256,
        "Failure record 002 SHA-256": FAILURE_002_SHA256,
        "Protocol Amendment 003 SHA-256": PROTOCOL_AMENDMENT_003_SHA256,
        "Config Amendment 003 SHA-256": CONFIG_AMENDMENT_003_SHA256,
        "Amendment 003 audit SHA-256": AMENDMENT_003_AUDIT_SHA256,
        "Protocol Amendment 004 SHA-256": PROTOCOL_AMENDMENT_004_SHA256,
        "Config Amendment 004 SHA-256": CONFIG_AMENDMENT_004_SHA256,
        "Amendment 004 audit SHA-256": AMENDMENT_004_AUDIT_SHA256,
        "Verification-replay 001 static audit SHA-256": VERIFICATION_REPLAY_001_STATIC_AUDIT_SHA256,
        "Verification-replay 002 static audit SHA-256": static_sha,
        "Verification-replay 002 static audit verdict": "CLEAR_TO_VERIFICATION_REPLAY_002",
        "Full output path": config["immutable_full_output"]["path"],
        "Full output file count": str(config["immutable_full_output"]["file_count"]),
        "Full run manifest SHA-256": FULL_MANIFEST_SHA256,
        "Failed verification output path": config["fixed_paths"]["failed_verification_output"],
        "Failed verification output status": "ABSENT",
        "Verification output path": config["fixed_paths"]["verification_output"],
    }
    for label, expected in required_decision_fields.items():
        require(field(decision_text, label) == expected, f"Decision 009 field mismatch: {label}")

    static_text = static_path.read_text(encoding="utf-8")
    required_static_fields = {
        "Verdict": "CLEAR_TO_VERIFICATION_REPLAY_002",
        "Supplemental verifier 002 SHA-256": own_sha,
        "Supplemental verifier 001 SHA-256": SUPPLEMENTAL_VERIFIER_001_SHA256,
        "Original verifier SHA-256": ORIGINAL_VERIFIER_SHA256,
        "Original runner SHA-256": RUNNER_SHA256,
        "Failure record 002 SHA-256": FAILURE_002_SHA256,
        "Protocol Amendment 004 SHA-256": PROTOCOL_AMENDMENT_004_SHA256,
        "Config Amendment 004 SHA-256": CONFIG_AMENDMENT_004_SHA256,
        "Amendment 004 audit SHA-256": AMENDMENT_004_AUDIT_SHA256,
        "Full run manifest SHA-256": FULL_MANIFEST_SHA256,
    }
    for label, expected in required_static_fields.items():
        require(field(static_text, label) == expected, f"verification-replay-002 static audit field mismatch: {label}")

    full = root / config["immutable_full_output"]["path"]
    entries = list(full.iterdir()) if full.is_dir() else []
    expected_files = set(config["immutable_full_output"]["expected_files"])
    require({path.name for path in entries} == expected_files and all(path.is_file() for path in entries), "immutable full-output inventory mismatch")
    require(full.resolve() == prior_full.resolve(), "immutable full-output path mismatch")
    require(sha256(full / "run_manifest.json") == FULL_MANIFEST_SHA256, "immutable full manifest changed")
    verification_output = root / config["fixed_paths"]["verification_output"]
    require(not failed_output.exists(), "failed verification output must remain absent before replay 002")
    require(not verification_output.exists(), "verification output 002 must be absent before replay")
    return prior, config_003, config, full, failed_output, verification_output, own_sha


def bit_preserving_f_order(original, expected_shape, c_strides, f_strides, label):
    require(original.dtype == np.float64, f"unexpected dtype before F-order conversion: {label}")
    require(list(original.shape) == list(expected_shape), f"unexpected shape before F-order conversion: {label}")
    require(original.flags.c_contiguous and not original.flags.f_contiguous, f"unexpected layout before F-order conversion: {label}")
    require(list(original.strides) == list(c_strides), f"unexpected strides before F-order conversion: {label}")
    converted = np.asfortranarray(original, dtype=np.float64)
    require(converted.shape == original.shape and converted.dtype == np.float64, f"shape or dtype changed during F-order conversion: {label}")
    require(converted.flags.f_contiguous and not converted.flags.c_contiguous, f"invalid layout after F-order conversion: {label}")
    require(list(converted.strides) == list(f_strides), f"unexpected strides after F-order conversion: {label}")
    require(np.array_equal(converted.view("<u8"), original.view("<u8")), f"F-order conversion changed values: {label}")
    return converted


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--decision-006-sha", required=True)
    parser.add_argument("--decision-007-sha", required=True)
    parser.add_argument("--decision-008-sha", required=True)
    parser.add_argument("--decision-009-sha", required=True)
    args = parser.parse_args()
    require(args.decision_006_sha == DECISION_006_SHA256, "Decision 006 argument mismatch")
    require(args.decision_007_sha == DECISION_007_SHA256, "Decision 007 argument mismatch")
    require(args.decision_008_sha == DECISION_008_SHA256, "Decision 008 argument mismatch")
    root = Path(__file__).resolve().parents[1]
    prior, config_003, config, full, failed_output, verification_output, own_sha = verify_authorization(root, args.decision_009_sha)
    verifier = prior.load_original_verifier(root)
    lowdim_layout_evidence = {}
    pxr_layout_evidence = {}
    original_lowdim_loader = verifier.load_biogen_contexts
    original_pxr_loader = verifier.load_pxr_context

    def lowdim_loader(*loader_args, **loader_kwargs):
        data, contexts, structure = original_lowdim_loader(*loader_args, **loader_kwargs)
        for endpoint, context in contexts.items():
            original = context["lowdim"]
            converted = bit_preserving_f_order(
                original,
                original.shape,
                config_003["root_cause"]["verifier_strides_bytes"][endpoint],
                config_003["root_cause"]["runner_strides_bytes"][endpoint],
                f"Biogen/{endpoint}/lowdim_11",
            )
            require(np.array_equal(context["representations"]["lowdim_11"].view("<u8"), original.view("<u8")), f"lowdim context aliases differ before conversion: {endpoint}")
            context["lowdim"] = converted
            context["representations"]["lowdim_11"] = converted
            lowdim_layout_evidence[endpoint] = {
                "shape": list(converted.shape),
                "dtype": str(converted.dtype),
                "preconversion_c_contiguous": True,
                "preconversion_strides_bytes": list(original.strides),
                "postconversion_f_contiguous": True,
                "postconversion_strides_bytes": list(converted.strides),
                "elementwise_bit_identity": True,
            }
        return data, contexts, structure

    def pxr_loader(*loader_args, **loader_kwargs):
        context = original_pxr_loader(*loader_args, **loader_kwargs)
        correction = config["sole_new_code_correction"]
        expected_keys = {"pxr_structure_10", "pxr_structure_screen_12"}
        require(set(context["features"]) == expected_keys, "unexpected PXR feature-context keys")
        for key in sorted(expected_keys):
            original = context["features"][key]
            record = config["root_cause"][key]
            names = correction["source_identity_keys"][f"features.{key}"]
            expected_source = context["expected_linkage"][names].to_numpy(dtype="<f8")
            require(expected_source.shape == original.shape and expected_source.dtype == np.float64, f"invalid source-derived PXR identity matrix: {key}")
            require(np.array_equal(expected_source.view("<u8"), original.view("<u8")), f"PXR source-derived feature identity mismatch: {key}")
            converted = bit_preserving_f_order(
                original,
                record["shape"],
                record["source_c_strides_bytes"],
                record["runner_f_strides_bytes"],
                f"PXR/{key}",
            )
            context["features"][key] = converted
            pxr_layout_evidence[key] = {
                "shape": list(converted.shape),
                "dtype": str(converted.dtype),
                "preconversion_c_contiguous": True,
                "preconversion_strides_bytes": list(original.strides),
                "postconversion_f_contiguous": True,
                "postconversion_strides_bytes": list(converted.strides),
                "elementwise_bit_identity": True,
                "source_expected_linkage_bit_identity": True,
            }
        return context

    verifier.load_biogen_contexts = lowdim_loader
    verifier.load_pxr_context = pxr_loader
    generator_storage = prior.verify_stored_generator_durability(verifier, full)
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
    require(set(lowdim_layout_evidence) == set(verifier.ENDPOINTS), "Biogen layout evidence incomplete")
    require(set(pxr_layout_evidence) == {"pxr_structure_10", "pxr_structure_screen_12"}, "PXR layout evidence incomplete")
    require(not failed_output.exists(), "failed verification output appeared during replay 002")
    result = {
        "status": "PASS",
        "mode": "verification-replay-002",
        "acceptance_state": "PASS_PENDING_POST_VERIFICATION_AUDIT",
        "scientific_decision": decision,
        "full_output_path": config["immutable_full_output"]["path"],
        "full_run_manifest_sha256": FULL_MANIFEST_SHA256,
        "decision_006_sha256": args.decision_006_sha,
        "decision_007_sha256": args.decision_007_sha,
        "decision_008_sha256": args.decision_008_sha,
        "decision_009_sha256": args.decision_009_sha,
        "original_runner_sha256": runner_sha,
        "original_verifier_sha256": original_verifier_sha,
        "supplemental_verifier_001_sha256": SUPPLEMENTAL_VERIFIER_001_SHA256,
        "supplemental_verifier_002_sha256": own_sha,
        "failure_record_001_sha256": FAILURE_001_SHA256,
        "failure_record_002_sha256": FAILURE_002_SHA256,
        "protocol_amendment_003_sha256": PROTOCOL_AMENDMENT_003_SHA256,
        "config_amendment_003_sha256": CONFIG_AMENDMENT_003_SHA256,
        "amendment_003_audit_sha256": AMENDMENT_003_AUDIT_SHA256,
        "protocol_amendment_004_sha256": PROTOCOL_AMENDMENT_004_SHA256,
        "config_amendment_004_sha256": CONFIG_AMENDMENT_004_SHA256,
        "amendment_004_audit_sha256": AMENDMENT_004_AUDIT_SHA256,
        "failed_verification_output_absent": True,
        "lowdim_layout_evidence": lowdim_layout_evidence,
        "pxr_layout_evidence": pxr_layout_evidence,
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
