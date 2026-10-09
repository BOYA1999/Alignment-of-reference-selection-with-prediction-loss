# Decision 003 — `GO_CONFIRM_001` after a verified mechanism signal

Date: 2026-09-11  
Decision state: `GO_CONFIRM_001`  
Canonical action: `build_and_run_after_code_audit`  
Formal confirmatory run authorized: **yes, subject to the frozen smoke, static-review, and empty-directory gates**

## Decision

Run the prospectively frozen confirmation study `confirm_001`, using the public Biogen ADME collection as the sole confirmation dataset. Do not reopen or relabel `pilot_001`. External-source readiness is `GO_DATA_READY`, the exact-claim collision addendum is `CLEAR_NARROW_CLAIM`, and the independent protocol audit is `CLEAR_TO_FREEZE`. The protocol and config are therefore frozen before any Biogen prediction, missingness realization, nuisance estimate, risk estimate, or PXR overlap model is generated.

The paper-facing question is narrowed to a domain-specific measurement-process validation problem:

> When ADMET reference measurements are selectively available, can known-truth sentinel controls identify unsafe propensity-based fixed-model risk correction before sparse-endpoint estimates are interpreted?

This is an Original Article route based on new computational experiments. It is not a review, a new AIPW/DML estimator, or a general ADMET validation framework.

## Decisive evidence

`forensic_001` was run after `pilot_001 = STOP_PILOT_NEGATIVE` under a pre-result protocol and passed both integrity checks:

- F1 frozen-pilot reproduction: `PASS`;
- F2 row-level formula replay: `PASS`;
- diagnostic decision: `MECHANISM_SIGNAL_DETECTED`.

The two mechanism signals carried into confirmation are:

- `PROPENSITY_INSTABILITY_SUPPORTED = true`: oracle/fitted AIPW median absolute-estimation-error ratios were `0.1451619170` under MCAR and `0.3680465094` under MAR;
- `HIGH_DIMENSION_INSTABILITY_SUPPORTED = true`: low-/full-dimensional fitted-AIPW absolute-estimation-error ratios were `0.1352882733` and `0.2223513463`, with median clipping reductions `0.5958525346` and `0.6988479263`.

The rejected mechanism is:

- `AUGMENTATION_INSTABILITY_SUPPORTED = false`: oracle-AIPW/oracle-IPW ratios were `1.1276780123` and `1.2489069933`, below the frozen `>2` threshold.

Therefore only propensity estimation and propensity-representation dimensionality are confirmatory. Loss-model representation and all propensity-by-loss interactions remain factorial controls and exploratory analyses. Augmentation instability may not be promoted after seeing confirmation results.

Evidence identities:

- `artifacts/diagnostics/forensic_001/protocol.md` — `e85bb5e05863f5667497fcadf2f8f38312ab82ab85120856e9505fa182a6b2c4`;
- `artifacts/diagnostics/forensic_001/config.json` — `72dfed2fd134fff7b5e6ccc6d6f0300659def6bd0d6d6f510aec4372a929c203`;
- `artifacts/diagnostics/forensic_001/full/summary.json` — `a26cb40c3fd453e9b2b9dc6f3d0108be889f9a3c62115f5a2083fcf23c953b52`;
- `artifacts/diagnostics/forensic_001/full/scenario_decomposition.csv` — `d0a0bd60b9b88f947c01a2a1e90ecce080b6d45ad9eeaf943e38fd0afcfbd49c`;
- `artifacts/diagnostics/forensic_001/full/row_replay_seed31001.parquet` — `af86a61856527e71a3dd293154b1a243bb7e4b0f983eae2f349f172e7dd1a2c0`;
- `artifacts/diagnostics/forensic_001/full/verification.md` — `2c91d5d8abc29a8a337cc2c908e509d89726537db32a234e852ea9778561c2c3`.

The broad live collision audit is `REVISE_AND_RECHECK`, not an unconditional novelty pass. Its frozen identity is `artifacts/scout/reference_availability_collision_audit_20260911.md`, SHA-256 `5c65ac9f5519c774a7ea576c2824a3ba8e61fba2249f5d1258d9b0a967b3f303`. The exact-claim addendum `artifacts/scout/confirm_001_exact_claim_collision_addendum_20260911.md`, SHA-256 `d141c86b381c4febdc050cf47713afe847bad84a01f352f50bfd30ea24149350`, returns `CLEAR_NARROW_CLAIM`: no directly isomorphic study was found, but every statistical component has precedent. The only defensible novelty unit is the linked ADMET measurement-process validation and abstention architecture; estimator, factorial, broad validation-framework, and universal-performance novelty claims remain closed.

## Dataset roles

1. **ExpansionRx — discovery only.** `pilot_001` and `forensic_001` define the hypothesis. No ExpansionRx row enters a confirmatory interval, test, threshold, or external-replication count.
2. **Biogen public ADME — sole confirmation.** The fixed-model risk, known-truth MCAR/MAR controls, real-mask transplantation, and 2 x 2 nuisance factorial are run here only after source freeze.
3. **OpenADMET PXR — boundary only.** It documents staged screen-to-dose-response selection, overlap loss, and non-identifiability. It cannot be used as evidence that full-library dose-response risk was recovered.

## Frozen external sources

The authoritative data-freeze manifest is `data/external/source_freeze_20260911/source_manifest.json`, SHA-256 `be55555aed0b2a5bed095ae8d052500becd527c97c809d6db5b7a14e7548fd78`; the human-readable audit is `artifacts/scout/external_source_freeze_20260911.md`, SHA-256 `79ecee56b1ca62ea38aa3bdfaff96e9e55213a0cc70c87d8dfd400a506d54bc9`. Its decision is `GO_DATA_READY`.

Biogen is pinned to official revision `b00df003de117ce9e5b381afd886095c5f2af2d5`. Official `ADME_public_set_3521.csv` is 388,279 bytes, 3,521 x 10, SHA-256 `2cfabc2667740c224487876c33b23124159ef43294e0f9e4d926cb6276c95a3b`; its pinned MIT `LICENSE` is SHA-256 `06557ca6eeee99c162df4e5e819167c6103f35faffb51f7b0a4833905280a740`. Eight analysis fields reconcile cell-for-cell with zero mismatches to the frozen derivative `evidence/upstream/expansion_ml_comparison/biogen_adme_3521.csv`, SHA-256 `cd310dd3911a9f9763d0581e6202bc20aaf8a810523d01e9ab5e403681a386b4`.

The formal analysis must consume that exact derivative. The BitBIRCH/split chain at construction revision `06582771d48af8096f6b05a5a247b177189a8fa4` is provenance only: its cluster and split vectors independently replayed exactly under the pinned construction environment, with hashes `c85847e68f2d726c325b0255722180201c0cca464fe0e2cbabbdb35734020a79` and `8ea3ebb25b17919946dd12111f4aa9d1f8f6c619b7ff9e28b1752501bbb89a79`. The model and nuisance analysis still uses RDKit `2026.03.3`; it must not regenerate `ds` under the analysis environment.

PXR is pinned to official revision `8041adfc295cac7f0bf85b73a2d87f21c5302ab1`, displayed license CC BY 4.0 with no standalone license file. The two analysis CSV hashes are `fdc2252bb5763ab5742e78cb1ba692bebe131a972cf83aa1bd45e43d8c5bc223` and `e781f8405e919a248b9b832a7cb24e8a615d400ee9dc50beea0d70d8aeaa8c2f`. The exact parent filter is `compound_class == "Library" AND concentration_M == 8.251e-06`, yielding 10,747 unique `OCNT_ID` rows; 2,722 have dose-response follow-up. Availability is a mixed-route, outcome-informed advancement indicator, not an activity label.

## Why the main alternatives lose

- Repeating or extending `pilot_001` on ExpansionRx would be post-result rescue and would not provide independent confirmation.
- Treating augmentation as confirmatory contradicts its frozen negative forensic signature.
- Treating PXR as a correction-success cohort is not identified because most non-advanced compounds lack dose-response truth and the public files do not reconstruct every analogue, diversity, QC, or discretionary promotion branch.
- Adding a model leaderboard, foundation-model comparison, or broad AI review would dilute the measurement-process question and re-enter already crowded literature.

## Frozen transition to `GO_CONFIRM_001`

The scientific and source gates are complete:

1. external data and license freeze: `GO_DATA_READY`;
2. exact-claim collision addendum: `CLEAR_NARROW_CLAIM`;
3. independent protocol/config audit: `CLEAR_TO_FREEZE`;
4. protocol/config status: `FROZEN_PRE_RESULT`, with `formal_run_authorized = true`.

Two operational gates remain mandatory before launch: implement against synthetic fixtures and pass an independent static review/replay check; then verify that the formal result directory is absent or empty. These gates may stop the run but may not change the frozen scientific contract.

Any later source drift, schema mismatch, license contradiction, cross-split duplicate structure, or inability to consume the exact frozen derivative returns the route to data audit. It is not repaired by silently substituting a source, split, endpoint, mask, model, or seed.

## Final pre-result freeze identities

- Frozen protocol SHA-256: `1819fa04713ee1a37df66c847825e5de3adcba736fba47c5c40b1c8db7796d7a`;
- Frozen config SHA-256: `04fc27fc407a562853dd362d31a4c4cbc01bdb754199d5a18014cdfcbc7a5c23`;
- Independent protocol-audit SHA-256: `a4b4c630785c2cd6991872696d2bc2cf23a14235cc19c7331ee1f70731fb5b16`;
- Exact-claim collision-addendum SHA-256: `d141c86b381c4febdc050cf47713afe847bad84a01f352f50bfd30ea24149350`;
- External source-manifest SHA-256: `be55555aed0b2a5bed095ae8d052500becd527c97c809d6db5b7a14e7548fd78`.

The capitalization of hexadecimal hashes is non-semantic; byte identity is authoritative. No `confirm_001` result existed when these identities were frozen.

## Next action

Implement the minimal runner and independent verifier against this frozen contract, pass the operational gates, and run `confirm_001` once. `dingzhen.md` remains unchanged until those artifacts and their QA exist.
