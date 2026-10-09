# Independent verification of `forensic_001`

Verified: 2026-09-11

## Result

- Status: `PASS`
- Adjudicative integrity: `PASS`
- Frozen diagnostic decision: `MECHANISM_SIGNAL_DETECTED`
- Scenario rows: 640
- Seed-31001 replay rows: 35,520
- F1 frozen-pilot reproduction: `PASS`
- F2 row-formula replay: `PASS`
- Maximum F1 absolute difference: `1.1102230246251565e-16`
- Maximum F2 absolute difference: `9.367506770274758e-17`
- Frozen inputs verified: 12
- Formal outputs verified: 4

## Frozen identities

- Diagnostic runner SHA-256: `d465ba4b162d64a52d807fa520efffd3a944a8dcfcb43d9289b6ce9bcc9044d8`
- Independent verifier SHA-256: `4aa517a2473b5bb37644cbb7350036646c55f765cea62bc8bc9510151e0d6ddc`
- Run manifest SHA-256: `73009b973e7b4de88184f711cd072a75e168bbaa20c811c968dde59a22f3aa55`

## Independent scope

`src/verify_forensic.py` does not import the diagnostic runner, the pilot runner, or their helper functions. It independently reconstructs source-row alignment, prediction losses, RDKit scaffolds, endpoint-specific three-fold GroupKFold assignments, all 40 MCAR/MAR mask vectors, the seed-31001 oracle probabilities, row-level augmentation and pseudo-outcomes, Hajek IPW, AIPW, scaffold-cluster standard errors, F1/F2, pooled mechanism summaries, all three prespecified signatures, the overall decision, and input/output hashes.

The verifier does not independently refit the nuisance models. Nuisance fitting is covered by the pre-run source audit and the stored row-level replay. This verification establishes artifact integrity and exact agreement with the frozen protocol; it is not independent biological validation and does not change `pilot_001 = STOP_PILOT_NEGATIVE`.
