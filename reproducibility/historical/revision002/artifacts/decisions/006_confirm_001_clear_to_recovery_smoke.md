# Decision 006 — `confirm_001` recovery-smoke clearance

Decision date: 2026-09-11

Decision: CLEAR_TO_RECOVERY_SMOKE
Mode: recovery-smoke
Runner SHA-256: 650a8b8f6b9bf65db69efc7ba60340aafa99c3c4b8b5bba957a8e6c6fa800c5c
Verifier SHA-256: 740e31590ec40f40f7783e59889d265c3d5d20296901c3968d4b593986238715
Decision 004 SHA-256: 0e9ab1a55bbe4fa99409c17d1a2c8d9087b9f84f4eedc81ba819edc3b8af39d8
Consumed Decision 005 SHA-256: 753068d66e5ca6dfd78dbcd1e2a56a81be87ea7908da214eb1a7ddd68675c025
Failed runner SHA-256: 9ce5f27e25e043db622b550446683ca830ec5536cc3e3f54897ed13fd4b22499
Failed verifier SHA-256: fc26475e8b99e0c64da5ff6b1073bec82f9a59e1247bb8021223f7f9c04fbe37
Failed output path: artifacts/confirmatory/confirm_001/full
Failure record path: artifacts/confirmatory/confirm_001/failure_001.md
Failure record SHA-256: a75a4b853bf9ce5469873870d0375e077f2dd2be7d5e6981febf5abd41667267
Protocol Amendment 002 path: artifacts/confirmatory/confirm_001/protocol_amendment_002.md
Protocol Amendment 002 SHA-256: d34c1053da408a3e349c2ecdc20d1e49f8aa5cda9dd491746146d621c5988d9f
Config Amendment 002 path: artifacts/confirmatory/confirm_001/config_amendment_002.json
Config Amendment 002 SHA-256: 144af9c376929f24699452051abed2ddb5243ab8cdf47c8f4024e219ac278984
Amendment 002 audit path: artifacts/confirmatory/confirm_001/amendment_002_independent_audit.md
Amendment 002 audit SHA-256: 7a266a7614729ccf5bb2ff13890eabb6d2e6e32cc0ebc66ea21a34c8f1a821c3
Amendment 002 audit verdict: CLEAR_TO_IMPLEMENT_RECOVERY
Recovery static audit path: artifacts/confirmatory/confirm_001/recovery_001_static_code_audit.md
Recovery static audit SHA-256: a36f7202eb6fae23229cfa790af8166d3b607e458ece3280b42d6eed26a430cf
Recovery static audit verdict: CLEAR_TO_RECOVERY_SMOKE
Recovery smoke output path: artifacts/confirmatory/confirm_001/recovery_smoke_001
Recovery full output path: artifacts/confirmatory/confirm_001/full_recovery_001

## Basis

Three independent reviews passed on the final hash-bound recovery implementation: static and scientific-contract fidelity, authorization and verifier closure, and numerical round-trip fidelity. The failed full directory remains immutable evidence, Decision 005 is consumed, and every new recovery path was absent before this decision.

## Scope

This decision authorizes exactly one deterministic synthetic, non-adjudicative recovery smoke at the fixed path above, followed by independent verification at the amendment-fixed verification path. It does not authorize the recovery full computation, scientific adjudication, interpretation, writing, release, or submission. Recovery full requires a separate Decision 007 bound to the successful smoke payload, manifest, code identities, Decision 006, and independently written verification record.
