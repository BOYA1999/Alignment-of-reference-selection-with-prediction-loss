# Independent static code audit of `confirm_001`

Audit date: 2026-09-11  
Audit type: read-only scientific-contract and execution-gate review  
Verdict: CLEAR_TO_SMOKE
Runner SHA-256: 9ce5f27e25e043db622b550446683ca830ec5536cc3e3f54897ed13fd4b22499
Verifier SHA-256: fc26475e8b99e0c64da5ff6b1073bec82f9a59e1247bb8021223f7f9c04fbe37

This verdict clears only the frozen synthetic smoke stage. It is not `CLEAR_TO_FULL`, is not a scientific result, and does not replace Decision 004 or Decision 005.

## Audit scope and negative execution statement

The audit read the frozen base protocol/config, the immutable base source manifest, all three Amendment 001 files, the independent amendment audit, and the final runner/verifier source. It checked scientific-contract fidelity, exact grids, fail-closed two-stage authorization, write-once behavior, performance hazards, and verifier independence.

No prediction model, missingness generator, mask realization, nuisance model, bootstrap, PXR propensity model, smoke test, or full analysis was run. No frozen authority file was edited. An independent AST parse of both final Python files passed. At final inspection, the formal result directory, Decision 004, and Decision 005 were absent.

## Frozen identities

| Artifact | SHA-256 | Result |
|---|---|---|
| `artifacts/confirmatory/confirm_001/protocol.md` | `1819fa04713ee1a37df66c847825e5de3adcba736fba47c5c40b1c8db7796d7a` | PASS |
| `artifacts/confirmatory/confirm_001/config.json` | `04fc27fc407a562853dd362d31a4c4cbc01bdb754199d5a18014cdfcbc7a5c23` | PASS |
| `data/external/source_freeze_20260911/source_manifest.json` | `be55555aed0b2a5bed095ae8d052500becd527c97c809d6db5b7a14e7548fd78` | PASS |
| `artifacts/confirmatory/confirm_001/protocol_amendment_001.md` | `989eac8ebd6e26bdcb5cee454ec51f4c00e76ad93bf4b245b92e0a59e9ff4a49` | PASS |
| `artifacts/confirmatory/confirm_001/config_amendment_001.json` | `9008878d6a8c23aae059cf3b4aefcdc813bdd7387cbf1aa84fd6d38db620ec3a` | PASS |
| `data/external/source_freeze_20260911/source_manifest_amendment_001.json` | `5b6b3337a9b7b1a4d359062ab07807527accc37b5cb63798db69e5f23837dfe5` | PASS |
| `artifacts/confirmatory/confirm_001/amendment_001_independent_audit.md` | `cc08015a2a8ea126b90cfd0a53af1a90d548a09f890ce86371490c0bdbc158cf` | PASS |

The six corrected ordered-vector hashes, 217-name RDKit descriptor registry/hash, per-row string-hash rule, and six-file authority precedence are enforced independently in both programs.

## Static checklist

### Source, cohort, features, and prediction

- PASS — Frozen file hashes, schemas, source-to-derivative mapping, physical source-row retention, split counts, cluster isolation, and canonical train/test non-overlap are fail-closed (`src/run_confirmatory.py:284-366`; `src/verify_confirmatory.py:489-653`).
- PASS — Morgan radius 2 / 2,048-bit / chirality features, alphabetically ordered 217 RDKit descriptors, endpoint-specific maximum train similarity, and the exact `lowdim_11` and `full_2266` float64 nuisance representations are implemented. Prediction XGBoost retains the frozen float32 input convention.
- PASS — The two frozen prediction models are fitted once per endpoint with their fixed hyperparameters and no result-dependent tuning or model removal (`src/run_confirmatory.py:626-698`). The verifier independently reconstructs features and refits both prediction models before accepting stored losses.

### Folds, masks, nuisance models, and formulas

- PASS — Deterministic three-fold Bemis-Murcko `GroupKFold`, acyclic source-row singleton labels, physical row order, and one shared endpoint fold vector are enforced.
- PASS — Fold-specific `StandardScaler` is fitted on every nuisance-training fold; propensity uses frozen L2/liblinear settings and strict `[0.05,0.95]` clipping; loss nuisance uses frozen Ridge/LSQR settings and observed training rows only.
- PASS — `ConvergenceWarning` is fail-closed for Logistic and Ridge fits, and Logistic fits must report iterations strictly below `max_iter` (`src/run_confirmatory.py:447-500`; `src/verify_confirmatory.py:178-187`).
- PASS — Fresh `default_rng(seed)` scope, seeds 41001-41100, shared within-sentinel uniform realization, mask reuse across models/factorial cells, MCAR coverage, and 200-update coverage-matched MAR are faithful to the frozen contract.
- PASS — Complete case, outcome regression, Hajek IPW, AIPW augmentation/pseudo-outcome, centered influence, scaffold-cluster SE, confidence interval, ESS, SMD, clipping, and normalized-weight formulas match the protocol.
- PASS — A nuisance fold with fewer than 20 observations in either class makes every affected H1/H2 cell non-adjudicative. H1/H2 now require fold support, finite outputs, and a valid denominator before passing; bootstrap rows carry the same fail-closed status (`src/run_confirmatory.py:1290-1364,1531-1608`; `src/verify_confirmatory.py:972-1022,1146-1230`).

### Exact durable grids and scientific gates

| Durable table | Required rows | Static result |
|---|---:|---|
| `availability_map.csv` | 8 | PASS |
| `predictions.parquet` | 3,702 | PASS |
| `nuisance_fold_assignments.parquet` | 1,851 | PASS |
| `generator_vectors.parquet` | 7,404 | PASS |
| `semisynthetic_scenarios.csv` | 6,400 | PASS |
| `actual_mask_scenarios.csv` | 64 | PASS |
| `row_replay_seed41001.parquet` | 59,232 in exactly 64 cells | PASS |
| `factorial_contrasts.csv` | 56 | PASS |
| `bootstrap_intervals.csv` | 12 | PASS |
| `pxr_linkage_audit.csv` | 10,747 | PASS |
| `pxr_selection_strata.csv` | 2 | PASS |
| `pxr_overlap_summary.csv` | 2 | PASS |

- PASS — H1 and H2 use exactly 400 rows per endpoint-by-mechanism cell with `M=lowdim_11`; H1 selects only the frozen `P=full_2266` copy of oracle AIPW. H3 and the 2-by-2 factorial contrasts remain exploratory.
- PASS — One frozen 10,000-by-100 seed-index matrix is created once from seed 52001, hash-checked, and reused for all 12 seed-cluster bootstrap statistics.
- PASS — The 16-cell known-truth admission gate, all-eight direct overlap gate, direct translation gate, and decision precedence are fail-closed. Stress anchors cannot rescue the direct translation gate.
- PASS — PXR uses the exact parent filter and `OCNT_ID` join, reproduces the 2,375 exact-SMILES and 2,700 batch matches plus 351/23 disagreement cardinalities, implements both frozen propensity specifications, and permits only the abstention boundary rather than a full-library point estimate (`src/run_confirmatory.py:1004-1209`; `src/verify_confirmatory.py:715-957`).

### Verifier independence and durable replay

- PASS — The verifier does not import the runner or its helpers. It independently reconstructs source lineages, canonical structures, descriptors, fingerprints, folds, prediction models, MAR generators, PXR source/linkage expectations, aggregations, gates, decisions, manifests, and exact output inventories.
- PASS — For all 64 seed-41001 replay cells, it regenerates the mask/oracle propensity, independently refits both cross-fitted propensity representations and both multi-target Ridge loss nuisances from independently reconstructed float64 features/folds, compares `p_raw`, clipped `p`, and `mhat` at absolute tolerance `1e-12`, and then recomputes every required downstream row formula (`src/verify_confirmatory.py:1932-2125`).
- PASS — The remaining 99 seeds and direct/PXR nuisance fits are not independently refitted. This is explicit in verifier output; their exact grids, vector hashes, sharing invariants, finite values, all frozen masks/oracle quantities, summaries, gates, and seed-41001 representative nuisance path are checked. No claim of exhaustive duplicate model execution is permitted.

### Authorization, write-once behavior, and manifests

- PASS — Decision 004 must bind this exact static-audit path/hash/verdict, both final code hashes, and the fixed smoke output path. Smoke cannot write inside the formal directory.
- PASS — Decision 005 must bind Decision 004, both final code hashes, the fixed smoke directory, hashes of `smoke.json` and its run manifest, and a hashed independent smoke-verification JSON whose internal identities and evidence hashes are checked (`src/run_confirmatory.py:202-281`; `src/verify_confirmatory.py:396-485`).
- PASS — Full mode requires the sole frozen formal path to be absent and creates it once with `exist_ok=False`; no deletion, overwrite, resume, or alternate formal path exists (`src/run_confirmatory.py:2188-2218`). The full manifest remains `COMPLETED_PENDING_INDEPENDENT_FULL_VERIFICATION` until the separate verifier passes.
- PASS — Exact output schemas, Cartesian keys, row counts, no-extra-file inventory, source/code/authorization hashes, environment identity, and pre-manifest output hashes are independently checked.

## Non-blocking advisories

1. P2 reporting — Admission cells correctly fail closed through `denominator_valid=false` and `pass=false`, but they do not persist the literal `NONADJUDICATIVE_FAIL` label named in Protocol Section 11.1 when the complete-case denominator or ratio is invalid. This cannot create a favorable gate result.
2. P2 reporting — The Biogen SMD helper excludes zero-variance features internally but the semisynthetic/direct tables do not explicitly persist the excluded-feature list; PXR does. The verifier requires the complete 11-feature Biogen SMD key set, so an unexpected zero-variance feature stops verification rather than passing silently.

Bootstrap ratio rows with an invalid denominator store finite zero placeholders only alongside `status=NONADJUDICATIVE_DENOMINATOR`, `denominator_valid=false`, and `point_pass=false`; those placeholders must never be interpreted as numerical estimates.

## Runtime and operational boundary

The dominant formal loop contains about 9,600 fold-level nuisance fits, roughly 4,800 of them high-dimensional, plus direct-mask, PXR, and prediction fits. Caching bounds memory and avoids redundant prediction/scaling work, but wall time may be multiple hours. The formal directory is created before computation; a crash can therefore leave an empty or partial path that deliberately blocks retry. Do not delete or reuse such a path under the existing authorization—record the failure and obtain a new additive recovery decision.

## Clearance boundary

Static code review is complete with no unresolved blocker. The next permissible action is creation of Decision 004 bound to this report and the two final code hashes, followed by the one fixed-path, non-adjudicative synthetic smoke and its independent verification. Full execution remains prohibited until a hash-bound Decision 005 exists.
