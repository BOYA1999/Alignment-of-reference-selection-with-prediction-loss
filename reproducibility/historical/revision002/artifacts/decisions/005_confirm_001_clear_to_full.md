# Decision 005 — `confirm_001` full-run clearance

Decision date: 2026-09-11

Decision: CLEAR_TO_FULL
Mode: full
Runner SHA-256: 9ce5f27e25e043db622b550446683ca830ec5536cc3e3f54897ed13fd4b22499
Verifier SHA-256: fc26475e8b99e0c64da5ff6b1073bec82f9a59e1247bb8021223f7f9c04fbe37
Decision 004 SHA-256: 0e9ab1a55bbe4fa99409c17d1a2c8d9087b9f84f4eedc81ba819edc3b8af39d8
Smoke status: PASS
Independent verifier smoke status: PASS
Smoke output path: artifacts/confirmatory/confirm_001/smoke_001
Smoke JSON SHA-256: 58a1cca53fc8367737d64422a22e88e46d6292683c1bb713a7aa7a64c453df3a
Smoke run manifest SHA-256: 71747e5fee8b225afb5f9c5a24103feb39222b6e5eabc20b5a5cb7b36c588c3f
Independent smoke verification path: artifacts/confirmatory/confirm_001/smoke_001_independent_verification.json
Independent smoke verification SHA-256: 0c523838f19288d3e01c21edaab70207a7adecdbc3b5e1917ddd0522a05fa557

## Basis

The authorized fixed-path synthetic smoke exited successfully and remained explicitly non-adjudicative. Its seven numerical, fold-support, finiteness, and float64 checks passed. The independent verifier was run separately and reproduced the smoke contract while closing the hashes of Decision 004, both programs, the smoke payload, and the smoke manifest. Two additional read-only audits independently confirmed the evidence chain and the absent formal output directory.

The observed scikit-learn `FutureWarning` concerns a future API deprecation in a version newer than the frozen environment. It is not a convergence warning or a present numerical failure and does not alter this clearance.

## Scope and recovery boundary

This decision authorizes exactly one full `confirm_001` execution at the sole frozen formal directory, using the frozen inputs, environment, estimators, seeds, gates, runner, and verifier identified above. It does not predetermine a favorable scientific result. Negative, null, abstention, or stop outcomes must be retained exactly as computed.

The formal directory must be absent at launch and is write-once. If execution fails after creating that directory, it must not be deleted, overwritten, resumed, or reused under this decision; failure evidence must be recorded and any retry requires a new additive recovery decision.
