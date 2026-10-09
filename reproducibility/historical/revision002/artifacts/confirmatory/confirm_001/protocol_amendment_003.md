# `confirm_001` verifier-only recovery amendment 003

Created at: 2026-09-11T08:45:00+08:00

State: PENDING_INDEPENDENT_VERIFICATION_RECOVERY_AUDIT
Amendment class: IMPLEMENTATION_ONLY_VERIFIER_REPLAY

## Trigger and quarantined state

The Decision-007-authorized runner completed `artifacts/confirmatory/confirm_001/full_recovery_001`, but the independent verifier exited with `generator probability hash mismatch`. The immutable failure record is `independent_verification_failure_001.md`. The runner manifest remains `COMPLETED_PENDING_INDEPENDENT_FULL_VERIFICATION`; no scientific decision in the runner summary is accepted or interpretable until an independent verifier passes.

The completed full directory contains exactly 19 files and is frozen by run-manifest SHA-256 `8da9b9c1d0425fff118d983a8491aceb57ae8b6085300eacf18cd8e0d7e6eac4`. It may not be edited, appended to, overwritten, resumed, rerun, or used as an output destination. Decision 007 is consumed.

## Proven implementation defect

The runner and verifier independently reconstruct identical `float64` values for the 11-column generator input, but the runner's row-then-column indexing produces F-contiguous storage whereas the verifier's column-then-row indexing produces C-contiguous storage. `StandardScaler` reduction order and liblinear amplify that layout-only difference. Explicit F-contiguous conversion of the verifier's bit-identical source-derived matrix reproduces the runner's four MAR score vectors, intercept shifts, oracle probability vectors, and SHA-256 hashes bit-for-bit. Four MCAR vectors already match exactly. Across all four MAR cells, 100 frozen seeds, and 370,200 threshold comparisons, the two layouts caused zero availability-mask flips.

The failure remains binding because two MAR vector differences exceed the frozen `1e-12` absolute equality tolerance and all four exact hashes differ. Hash checks and tolerances may not be removed or relaxed.

## Sole permitted correction

A new supplemental verifier may import the unchanged, hash-bound `src/verify_confirmatory.py` and replace only its source-derived Biogen low-dimensional array at the generator and nuisance entry point with:

`np.asfortranarray(lowdim, dtype=np.float64)`

Before replacement it must prove that the source-derived matrix is C-contiguous; after replacement it must prove F-contiguity, `float64` dtype, unchanged shape, and elementwise bit identity. The same converted matrix must replace `context["lowdim"]` and `context["representations"]["lowdim_11"]`; no other context value may change.

The supplemental verifier must additionally:

1. bind the unchanged runner, original verifier, D006, D007, failure record, this amendment and its config/audits, the complete full-output manifest, its own code, and a new Decision 008;
2. independently rehash every stored generator probability vector from the Parquet float64 bytes;
3. reconstruct stored MAR probabilities from stored score plus shift through the frozen `expit` and `[0.10, 0.90]` clip and require bit identity;
4. run the original verifier's complete independent full verification, preserving every exact hash and `rtol=0, atol=1e-12` check;
5. create exactly one external verification record at `artifacts/confirmatory/confirm_001/full_recovery_001_independent_verification.json` using exclusive creation only after all checks pass.

It may not invoke the runner, write inside `full_recovery_001`, use any durable result as an expected source-derived generator vector, bypass a failed check, catch and convert verification exceptions to PASS, or create a verification record after failure.

## Authorization sequence and fixed paths

The following paths are initially absent and write-once:

- supplemental verifier: `src/verify_confirmatory_replay.py`
- amendment audit: `artifacts/confirmatory/confirm_001/amendment_003_independent_audit.md`
- supplemental static audit: `artifacts/confirmatory/confirm_001/verification_replay_001_static_audit.md`
- authorization: `artifacts/decisions/008_confirm_001_clear_to_verification_replay.md`
- successful verification record: `artifacts/confirmatory/confirm_001/full_recovery_001_independent_verification.json`

Implementation requires an independent amendment audit with verdict `CLEAR_TO_IMPLEMENT_VERIFIER_REPLAY`. Execution requires a static code audit with verdict `CLEAR_TO_VERIFICATION_REPLAY` and a separately hashed Decision 008. Decision 008 may authorize one read-only verification replay of the existing immutable full directory and no runner execution.

If the replay fails after launch, the existing full directory and all prior evidence remain unchanged, no PASS record may be written, and another additive failure record and decision are required.

## Unchanged scientific contract

All source bytes and readers, cohorts, endpoints, anchors, features and feature values, models and parameters, folds, masks, seeds, RNG scopes, estimands, formulas, thresholds, bootstrap indices, admission criteria, overlap criteria, decision precedence, anti-claims, output bytes, and negative/null-result retention remain unchanged. This amendment cannot turn the current quarantined runner decision into accepted evidence; only a complete independent verification PASS can do so.

This document authorizes neither implementation nor execution.
