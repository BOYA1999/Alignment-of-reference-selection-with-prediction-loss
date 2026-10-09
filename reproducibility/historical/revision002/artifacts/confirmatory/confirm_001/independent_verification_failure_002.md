# `confirm_001` independent verification failure 002

Recorded at: 2026-09-11T09:09:46+08:00
State: INDEPENDENT_VERIFICATION_REPLAY_001_FAILED_BEFORE_RESULT_ACCEPTANCE
Gate outcome: STOP_IMPLEMENTATION_OR_DATA
Attempt mode: verification-replay
Exit code: 1
First failing check: value mismatch: pxr_overlap_summary.csv/smd_json
Verification output created: false
Consumed Decision 008 SHA-256: 3cc6426bdebeb739d31d8d0cee1092d7de7ac8c4625b6b2004011390fec43216
Verification-replay static audit SHA-256: df4c666f97390258e164d7c65cda527d388d7e82f7fefbfc9ba44cd6010ffe35
Supplemental verifier SHA-256: 4e95a20d5a3c74d5c96122b80fc37cf5d44892c4438e846b781899fbf2212a57
Original verifier SHA-256: 740e31590ec40f40f7783e59889d265c3d5d20296901c3968d4b593986238715
Original runner SHA-256: 650a8b8f6b9bf65db69efc7ba60340aafa99c3c4b8b5bba957a8e6c6fa800c5c
Decision 006 SHA-256: 96cd052f18c9442ebe9b3b3524da7b88539c0896f6541dbc675cf9cb094cee0b
Consumed Decision 007 SHA-256: faea63b749016aaaf570888ca5cbc99fda856088c53ac02f489fca3df7f45fe7
Failure record 001 SHA-256: 3a6ac0a67988d3037831f9f21a084b814f78e701ec464c5e224b09018e74c8a1
Protocol Amendment 003 SHA-256: 458d0217da500b89a89c28018e07f2df18b3d741f1fb13e068451cef0a447adb
Config Amendment 003 SHA-256: 690875bca0afcd95a44946c57aa552c3287c1c80efb690cb412318565999ab60
Amendment 003 audit SHA-256: 61a67e578d57ce790182646c61ba8da18ce8f49d7fb9f6853f4033b21156fe75
Full run manifest SHA-256: 8da9b9c1d0425fff118d983a8491aceb57ae8b6085300eacf18cd8e0d7e6eac4
Full run summary SHA-256: 784132d1df5ab6a99a934bf2c53d17c7d11f650beb4b733be222011f060986ae
Full run verification narrative SHA-256: 7ef5b2cdb2d1885b074250402efd313aa9758b337b4f486147ddc10e0042af99
Generator vectors SHA-256: a6b817b9a3db207dd84f1a816a7374a5a21b186c853b73c2f2c422a1bfc830f7
Immutable full-output file count: 19
Immutable full-output directory count: 0
Frozen output hashes matching manifest: 18/18
Failed verification output path: artifacts/confirmatory/confirm_001/full_recovery_001_independent_verification.json

## Authorized command and failure

The Decision-008-authorized command was executed exactly once:

```text
python src/verify_confirmatory_replay.py --decision-006-sha 96cd052f18c9442ebe9b3b3524da7b88539c0896f6541dbc675cf9cb094cee0b --decision-007-sha faea63b749016aaaf570888ca5cbc99fda856088c53ac02f489fca3df7f45fe7 --decision-008-sha 3cc6426bdebeb739d31d8d0cee1092d7de7ac8c4625b6b2004011390fec43216
```

The supplemental verifier completed its authorization, stored-generator durability, F-order low-dimensional reconstruction, original environment, and frozen-input checks. The unchanged original `verify_full` flow then stopped at `compare_frames(expected_overlap, pxr_overlap, "pxr_overlap_summary.csv", ("specification",), tolerance)`. The exact emitted error was:

```text
{"status": "FAIL", "error": "value mismatch: pxr_overlap_summary.csv/smd_json"}
AssertionError: value mismatch: pxr_overlap_summary.csv/smd_json
```

The exception propagated, the process exited nonzero, and exclusive creation of the external PASS record was never reached. The failed verification-output path remains absent. Decision 008 is consumed and may not authorize another execution.

## Read-only root-cause isolation

The two independently source-reconstructed PXR feature matrices and the corresponding columns reloaded from the frozen round-trip linkage CSV are elementwise `float64` bit-identical. The source matrices are C-contiguous, with shapes and byte strides `(10747, 10) / [80, 8]` and `(10747, 12) / [96, 8]`; the runner's DataFrame-selected matrices are F-contiguous, with byte strides `[8, 85976]` for both specifications.

The runner computes each SMD with separate one-dimensional `std`, `sum`, and `mean` reductions. The original verifier computes all columns through C-order `axis=0` vectorized reductions. Different floating-point reduction order changes all 22 JSON-embedded SMD values only in their terminal bits. The maximum absolute discrepancy is `5.009881398621019e-15` for `pxr_structure_10` and `2.2752633110911802e-14` for `pxr_structure_screen_12`, both below the frozen `rtol=0, atol=1e-12` numerical tolerance. The stored values are reproduced bit-for-bit either by the runner's columnwise formula on the independent source matrix or by the unchanged verifier formula after a bit-preserving F-order conversion. The current verifier nevertheless rejects them because generic frame comparison treats `smd_json` as an exact string.

This is a layout-triggered verifier reduction-order defect. It is not CSV round-trip damage, a feature-value difference, an output-hash failure, or a runner calculation defect. The structure-only specification remains `trigger_clip=false`, `trigger_ess=false`, `trigger_smd=false`, `overlap_pass=true`; the structure-plus-screen specification remains `trigger_clip=true`, `trigger_ess=true`, `trigger_smd=true`, `overlap_pass=false`. The formal-hit count and hash, selection strata, and PXR boundary object are unchanged. These diagnostic facts do not accept or interpret the quarantined overall scientific result.

## Frozen aftermath and recovery boundary

After failure, `full_recovery_001` still contains exactly the original 19 files and no directories; all 18 pre-manifest output hashes match the unchanged manifest. The original runner, original verifier, supplemental verifier, all prior decisions and amendments, and all full-run outputs remain unchanged.

No replay, code change, result interpretation, or PASS record is authorized by this failure record. A further attempt requires a new verifier-only amendment, independent amendment and static-code audits, a new hash-bound supplemental verifier at a new path, a new external write-once verification path, and a separately hashed Decision 009. The original failed verification path must remain absent.
