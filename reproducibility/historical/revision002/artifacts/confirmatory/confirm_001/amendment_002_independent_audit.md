# Independent audit of `confirm_001` Recovery Amendment 002

Audit date: 2026-09-11  
Audit type: read-only post-failure recovery-contract audit  
Verdict: CLEAR_TO_IMPLEMENT_RECOVERY

Protocol Amendment 002 SHA-256: d34c1053da408a3e349c2ecdc20d1e49f8aa5cda9dd491746146d621c5988d9f  
Config Amendment 002 SHA-256: 144af9c376929f24699452051abed2ddb5243ab8cdf47c8f4024e219ac278984  
Failure record SHA-256: a75a4b853bf9ce5469873870d0375e077f2dd2be7d5e6981febf5abd41667267

This verdict authorizes only implementation and static audit of the bounded recovery correction. It does not authorize a recovery smoke, a recovered full run, scientific adjudication, interpretation, writing, release, or submission.

## Evidence identity and failed-state preservation

- PASS — The base protocol and config remain byte-identical at SHA-256 `1819fa04713ee1a37df66c847825e5de3adcba736fba47c5c40b1c8db7796d7a` and `04fc27fc407a562853dd362d31a4c4cbc01bdb754199d5a18014cdfcbc7a5c23`.
- PASS — Amendment 001, its source-manifest overlay, and its independent audit remain unchanged. Recovery Amendment 002 is additive and preserves their precedence.
- PASS — Decision 004 and consumed Decision 005 remain unchanged at SHA-256 `0e9ab1a55bbe4fa99409c17d1a2c8d9087b9f84f4eedc81ba819edc3b8af39d8` and `753068d66e5ca6dfd78dbcd1e2a56a81be87ea7908da214eb1a7ddd68675c025`.
- PASS — The failed runner and verifier snapshots reproduce the failed code identities exactly: runner `9ce5f27e25e043db622b550446683ca830ec5536cc3e3f54897ed13fd4b22499`; verifier `fc26475e8b99e0c64da5ff6b1073bec82f9a59e1247bb8021223f7f9c04fbe37`.
- PASS — `artifacts/confirmatory/confirm_001/full/` still contains exactly the 16 files recorded in `failure_001.md`; every byte count and SHA-256 matches the machine-readable Amendment 002 inventory. It has no nested directory or reparse point. `summary.json`, `verification.md`, and `run_manifest.json` remain absent.
- PASS — The failure record transparently reports exit code 1, the available traceback tail, the absence of a predesignated raw stdout/stderr log, the two one-ULP ratio discrepancies, and the boundary `FAILED_BEFORE_SCIENTIFIC_DECISION`.

The failed directory is immutable evidence only. It may not be deleted, moved, modified, resumed, completed, or used as scientific input by the recovered computation.

## Recovery correction and scientific invariance

- PASS — The allowed correction is limited to generated-result CSV I/O: `float_format="%.17g"` on writes and `float_precision="round_trip"` on runner/verifier reads. Frozen source/input CSV readers and source bytes remain unchanged.
- PASS — Formula comparison remains `rtol=0` and `atol=1e-12`; no tolerance, row exemption, estimator, seed, gate, or decision rule is relaxed.
- PASS — The required recovery smoke adds a deterministic synthetic large-ratio CSV fixture and remains non-adjudicative. Failed scientific outputs may not supply that fixture.
- PASS — All cohorts, endpoints, row orders, masks, features, prediction and nuisance models, folds, RNG scopes, estimands, formulas, bootstrap indices, thresholds, grids, decision precedence, anti-claims, and negative/null-result retention remain frozen.
- PASS — Optional suppression of the named scikit-learn API `FutureWarning` is bounded to logging hygiene; `ConvergenceWarning` and iteration checks must remain fail-closed.

## New paths and two-stage authorization

At audit time, all three recovery evidence/output paths were absent:

- `artifacts/confirmatory/confirm_001/recovery_smoke_001/`
- `artifacts/confirmatory/confirm_001/recovery_smoke_001_independent_verification.json`
- `artifacts/confirmatory/confirm_001/full_recovery_001/`

They are fixed sibling paths and must be created at most once. The repaired runner and verifier must independently recheck the immutable failed inventory and the absence of the applicable new path before creating or accepting recovery outputs.

Decision 006 may be created only after repaired runner/verifier hashes are final and an independent static audit returns `CLEAR_TO_RECOVERY_SMOKE`. It may authorize only the fixed, synthetic, non-adjudicative recovery smoke and must bind this audit, both Amendment 002 files, the failure record, old D004/D005, failed-code snapshots, final repaired code hashes, and the fixed recovery paths.

Decision 007 may be created only after the recovery smoke and its independent verifier both PASS. It must bind Decision 006, the repaired code and static-audit identities, the recovery smoke payload and manifest hashes, the independent smoke-verification path/hash/status, and the still-absent `full_recovery_001` path. Decision 005 is consumed historical evidence and cannot authorize the recovered full execution.

## Implementation-audit requirements

The next permitted action is the minimal recovery implementation followed by an independent diff-based static audit. That audit must confirm that generated-result CSV handling, recovery path/authorization constants, failed-inventory checks, and the non-adjudicative round-trip smoke fixture are the only material changes; source readers and every scientific computation must remain unchanged. No execution is authorized until the separate Decision 006 gate exists and its trusted SHA-256 is supplied.
