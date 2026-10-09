# `confirm_001` Post-failure Implementation Recovery Amendment 002

Timestamp: 2026-09-11T07:31:30+08:00  
State: `PENDING_INDEPENDENT_RECOVERY_AMENDMENT_AUDIT`

## Trigger and current scientific state

The first Decision-005-authorized full execution created `artifacts/confirmatory/confirm_001/full/` and stopped with exit code 1 during the runner's durable CSV reload check. The immutable failure record is `failure_001.md`. The failed directory contains 16 partial files but no `summary.json`, `verification.md`, or `run_manifest.json`; integrity PASS, H1/H2 adjudication, and the overall scientific decision were not reached.

The observed exception was `RuntimeError: factorial contrasts numeric mismatch: ratio`. Read-only diagnosis localized the failure to one-unit-in-last-place decimal parsing differences: two exploratory large-ratio rows exceeded the frozen `1e-12` absolute reload tolerance by `1.4551915228366852e-11` and `7.275957614183426e-12`, while their relative differences were approximately `1.43e-16` and `1.28e-16`. Independent diagnosis also found that the same default CSV parser could exceed the absolute tolerance in stored scenario formula fields. This is an implementation-level durable-serialization failure, not a scientific hypothesis result.

## Minimal implementation correction

The correction applies only to result CSV serialization and result CSV reload:

1. Every runner-generated result CSV is written with `float_format="%.17g"`.
2. Every runner or verifier read of a generated result CSV uses `float_precision="round_trip"`.
3. Frozen source and analysis-input CSVs continue to use their original readers and byte identities. Their parsing contract is not changed.
4. Formula replay keeps the same zero-relative, `1e-12` absolute tolerance. No tolerance is relaxed and no row is exempted.
5. The explicit scikit-learn logistic API deprecation warning may be suppressed inside the already fail-closed fit wrapper; `ConvergenceWarning` remains an error and iteration-count checks remain unchanged.

The recovery smoke must add a deterministic, synthetic CSV round-trip fixture containing small denominators and ratios of order `1e5`. It must demonstrate exact restoration of the component floats and a derived-ratio absolute difference no greater than `1e-12`. The fixture is non-adjudicative and may not use the failed scientific results as input.

## Immutable failure evidence and new paths

The failed `artifacts/confirmatory/confirm_001/full/` directory and its 16 files must remain in place and byte-identical to the inventory in `failure_001.md`. It may not be deleted, moved, overwritten, resumed, completed, or used as input to the recovered computation.

Recovery uses two new, initially absent sibling paths:

- synthetic gate: `artifacts/confirmatory/confirm_001/recovery_smoke_001/`
- full computation: `artifacts/confirmatory/confirm_001/full_recovery_001/`

Both paths are write-once. Any failure after creating either path requires another additive failure record and recovery decision; the path may not be reused.

## Two-stage recovery authorization

Decision 006 may authorize only the fixed recovery smoke after an independent audit binds the recovery amendment, immutable failure evidence, old consumed decisions and code snapshots, final repaired runner/verifier hashes, exact recovery paths, and a static code verdict of `CLEAR_TO_RECOVERY_SMOKE`.

Decision 007 may authorize the fixed recovery full run only after the smoke and its independent verifier both PASS and their payload, manifest, code, Decision 006, and verification-record hashes are closed. Decision 005 is consumed evidence and cannot authorize another full execution.

## Unchanged scientific contract

The recovery must recompute from the frozen sources. It may not read any partial failed-result file for scientific computation. All source bytes, source hashes, row orders, cohorts, endpoints, target/anchor roles, features, prediction models, nuisance models, folds, masks, seeds, RNG scopes, estimands, formulas, bootstrap indices, admission criteria, overlap criteria, thresholds, expected grids, decision precedence, anti-claims, and negative/null-result retention remain unchanged.

This amendment does not authorize implementation or execution. Independent amendment audit, repaired-code audit, Decision 006 smoke clearance, smoke PASS plus independent verification, and Decision 007 full clearance are all required in that order.
