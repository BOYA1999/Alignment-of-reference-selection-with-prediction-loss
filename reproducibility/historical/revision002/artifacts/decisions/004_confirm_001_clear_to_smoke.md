# Decision 004 — `confirm_001` synthetic smoke clearance

Decision date: 2026-09-11

Decision: CLEAR_TO_SMOKE
Mode: smoke
Runner SHA-256: 9ce5f27e25e043db622b550446683ca830ec5536cc3e3f54897ed13fd4b22499
Verifier SHA-256: fc26475e8b99e0c64da5ff6b1073bec82f9a59e1247bb8021223f7f9c04fbe37
Smoke output path: artifacts/confirmatory/confirm_001/smoke_001
Static audit path: artifacts/confirmatory/confirm_001/static_code_audit_001.md
Static audit SHA-256: 84082de50824deba4cd9f77425b1b00e3678916f9fd4fc7e6a62f8c410e8f39c
Static audit verdict: CLEAR_TO_SMOKE

## Basis

The final runner and independent verifier compile successfully and the hash-bound static audit found no blocker in the frozen source, cohort, feature, estimator, gate, authorization, manifest, or write-once contracts. Independent numerical and verifier reviews reached the same clearance conclusion. The audit's two reporting advisories cannot create a favorable scientific decision and remain recorded without changing the frozen analysis.

## Scope

This decision authorizes exactly one synthetic, non-adjudicative smoke run at the fixed output path above, followed by independent smoke verification. It does not authorize the formal computation, scientific interpretation, manuscript drafting, public release, or submission. A separate Decision 005 bound to the successful smoke evidence is required before full execution.
