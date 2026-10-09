# `confirm_001` verifier-only recovery amendment 004

Created at: 2026-09-11T09:12:00+08:00

State: PENDING_INDEPENDENT_VERIFICATION_RECOVERY_AUDIT
Amendment class: IMPLEMENTATION_ONLY_VERIFIER_REPLAY_002

## Trigger and quarantined state

The Decision-008-authorized supplemental verifier replay exited with `value mismatch: pxr_overlap_summary.csv/smd_json`. The immutable record is `artifacts/confirmatory/confirm_001/independent_verification_failure_002.md`, SHA-256 `455c8e3928099fb31452801b5cec1c12d3ab716eaba1a03dd4ca8ee59a94e634`. Decision 008 is consumed. Its intended external verification path remains absent and must remain permanently absent.

The completed full directory remains frozen at exactly 19 files and no directories, with run-manifest SHA-256 `8da9b9c1d0425fff118d983a8491aceb57ae8b6085300eacf18cd8e0d7e6eac4` and all 18 pre-manifest output hashes matching. It may not be edited, appended to, overwritten, resumed, renamed, deleted, rerun, or used as an output destination. The runner decision remains quarantined and scientifically uninterpretable.

## Proven implementation defect

For each PXR specification, the verifier's independently source-reconstructed feature matrix and the matching feature columns reloaded from the frozen round-trip linkage CSV are elementwise `float64` bit-identical. The independently reconstructed matrices are C-contiguous: `(10747, 10)` with byte strides `[80, 8]` and `(10747, 12)` with byte strides `[96, 8]`. The runner's DataFrame-selected matrices are F-contiguous with byte strides `[8, 85976]` in both cases.

The runner computes SMDs through separate one-dimensional column reductions. The original verifier computes the same formula through `axis=0` vectorized reductions on the C-order source matrix. Layout-triggered reduction order changes all 22 JSON-embedded SMD values only in terminal bits: the maximum absolute discrepancy is `5.009881398621019e-15` for `pxr_structure_10` and `2.2752633110911802e-14` for `pxr_structure_screen_12`. These are below the unchanged frozen `rtol=0, atol=1e-12` tolerance, but generic frame comparison treats `smd_json` as an exact string and therefore fails.

This is not CSV round-trip damage or a feature-value difference. The stored SMD values are reproduced bit-for-bit by the frozen runner columnwise formula applied to either bit-identical matrix and by the unchanged verifier formula applied after F-order conversion. The PXR clip, ESS and SMD triggers, both overlap-pass values, formal-hit vector, selection strata and boundary object remain unchanged. This diagnosis concerns only the first failing verifier check and does not predict that a new complete replay will pass.

## Sole new permitted correction

A new supplemental verifier at `src/verify_confirmatory_replay_002.py` may import the unchanged, hash-bound original verifier and preserve the Amendment-003 Biogen low-dimensional F-order correction. In addition, immediately after the original `load_pxr_context` independently reconstructs its source context, it may replace only these two arrays:

- `context["features"]["pxr_structure_10"]`
- `context["features"]["pxr_structure_screen_12"]`

Each replacement must be produced solely as `np.asfortranarray(original, dtype=np.float64)`. Before replacement, the verifier must require the frozen shape, `float64` dtype, C-contiguity, non-F-contiguity and exact C-order strides. It must independently require bit identity between the original array and the corresponding columns of the source-derived `expected_linkage` object. After replacement, it must require unchanged shape and dtype, F-contiguity, exact `[8, 85976]` strides and elementwise `uint64` bit identity with the original. The two converted arrays must replace only the two named keys. No stored SMD or other durable full-run result may supply an expected feature value or be used to construct a replacement.

The new verifier must retain the Amendment-003 low-dimensional conversion and its stored-generator durability checks unchanged. It must run the unchanged original verifier's complete full flow, including the original exact `smd_json` comparison, every exact hash, and every `rtol=0, atol=1e-12` numerical check. It may not patch or replace `compare_frames`, `smd_map`, `verify_pxr_linkage`, any threshold, or any scientific decision function.

## Binding, write boundary and failure behavior

The new verifier must bind by SHA-256 the unchanged runner and original verifier; the first supplemental verifier; D006, D007 and consumed D008; failure records 001 and 002; Amendments 002, 003 and 004 with their configs and audits; both earlier static audits; the complete full-output manifest; its own code and static audit; and a new Decision 009.

It may not import or invoke the runner. Before computation it must require both the failed first verification path and the new verification path to be absent, and it must recheck the exact 19-file full-output inventory and manifest hash. Only after all supplemental checks and the unchanged original `verify_full` return successfully may it create exactly one external record at `artifacts/confirmatory/confirm_001/full_recovery_001_independent_verification_002.json` using exclusive creation. Any exception must propagate as failure and must not create a PASS record.

If the replay fails, D009 is consumed, both the frozen full directory and all prior evidence remain unchanged, no PASS record may be written, and another additive failure record and authorization chain are required. No retry under D009 is permitted.

## Authorization sequence and fixed paths

The following new paths are initially absent and write-once:

- `artifacts/confirmatory/confirm_001/protocol_amendment_004.md`
- `artifacts/confirmatory/confirm_001/config_amendment_004.json`
- `artifacts/confirmatory/confirm_001/amendment_004_independent_audit.md`
- `src/verify_confirmatory_replay_002.py`
- `artifacts/confirmatory/confirm_001/verification_replay_002_static_audit.md`
- `artifacts/decisions/009_confirm_001_clear_to_verification_replay_002.md`
- `artifacts/confirmatory/confirm_001/full_recovery_001_independent_verification_002.json`

Implementation requires an independent Amendment-004 audit with verdict `CLEAR_TO_IMPLEMENT_VERIFICATION_REPLAY_002`. Execution requires a separate static-code audit with verdict `CLEAR_TO_VERIFICATION_REPLAY_002` and a separately hashed Decision 009. Decision 009 may authorize exactly one read-only verification replay and no runner execution.

## Unchanged scientific contract

All source bytes and readers, cohorts, endpoints, anchors, feature names and values, models and parameters, folds, masks, seeds, RNG scopes, estimands, formulas, thresholds, bootstrap indices, admission criteria, overlap criteria, decision precedence, anti-claims, output bytes, and negative/null-result retention remain unchanged. No tolerance, exact hash, JSON comparison, trigger or decision rule is removed or relaxed. This amendment authorizes neither implementation nor execution and cannot accept or interpret any scientific result.
