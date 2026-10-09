# `confirm_001` formal execution failure record 001

Record date: 2026-09-11

Status: FAILED_BEFORE_SCIENTIFIC_DECISION
Exit code: 1
Failed output path: artifacts/confirmatory/confirm_001/full
Failed output created UTC: 2026-09-10T23:13:07.6964668Z
Last failed-output write UTC: 2026-09-10T23:21:37.9403533Z
Decision 004 SHA-256: 0e9ab1a55bbe4fa99409c17d1a2c8d9087b9f84f4eedc81ba819edc3b8af39d8
Decision 005 SHA-256: 753068d66e5ca6dfd78dbcd1e2a56a81be87ea7908da214eb1a7ddd68675c025
Failed runner SHA-256: 9ce5f27e25e043db622b550446683ca830ec5536cc3e3f54897ed13fd4b22499
Failed verifier SHA-256: fc26475e8b99e0c64da5ff6b1073bec82f9a59e1247bb8021223f7f9c04fbe37

## Command and terminal evidence

The authorized command was:

```text
python src/run_confirmatory.py --mode full --decision-004-sha 0e9ab1a55bbe4fa99409c17d1a2c8d9087b9f84f4eedc81ba819edc3b8af39d8 --decision-005-sha 753068d66e5ca6dfd78dbcd1e2a56a81be87ea7908da214eb1a7ddd68675c025
```

The process returned exit code 1. No raw stdout/stderr log path had been designated before launch, so a byte-exact full terminal log is unavailable. The captured terminal tail showed this call chain and final exception:

```text
main
  -> run_full (line 1968)
  -> durable_integrity (line 1860)
  -> compare_frames (line 1727)
  -> require (line 168)
RuntimeError: factorial contrasts numeric mismatch: ratio
```

Repeated scikit-learn `FutureWarning` messages about the future deprecation of the explicit logistic `penalty` argument preceded the exception. No `ConvergenceWarning` or model-fit exception was observed.

## Failure boundary

The runner stopped during its reload-and-compare durability check. It did not write `summary.json`, `verification.md`, or `run_manifest.json`; consequently there is no completed formal result, no integrity PASS, no H1/H2 adjudication, and no frozen scientific decision. The partial directory is retained byte-for-byte and must not be deleted, overwritten, resumed, or treated as a completed result.

## Read-only diagnosis

Recomputing `factorial_contrasts.csv` from the stored `semisynthetic_scenarios.csv` reproduced exactly two absolute-tolerance failures in the exploratory augmentation-ratio rows under the frozen `1e-12` absolute comparison:

| endpoint | mechanism | loss representation | absolute ratio difference | relative ratio difference |
|---|---|---|---:|---:|
| `LOG_RLM` | `mar` | `full_2266` | `1.4551915228366852e-11` | approximately `1.43e-16` |
| `LOG_RLM` | `mcar` | `full_2266` | `7.275957614183426e-12` | approximately `1.28e-16` |

The corresponding numerator or denominator round-trip differences were at most about `5.1e-17`; the ratios were approximately `101480.84` and `56982.39`, so division by a small exploratory denominator amplified a one-unit-in-last-place serialization difference beyond the absolute threshold. No threshold, estimator, mask, seed, model, hypothesis, or scientific gate was changed or examined to obtain this diagnosis.

## Partial output inventory

Exactly 16 files were present after failure:

| file | bytes | SHA-256 |
|---|---:|---|
| `actual_mask_scenarios.csv` | 84954 | `15b89765bf9d3247645a8b97c5bf21bf2b3447622aaca2ab02827896847a02d4` |
| `admission_gate.json` | 14580 | `3bbb6886a6713ef43ead4463acdf81c11d339de775ecce079c0be38fa7fb0e51` |
| `availability_map.csv` | 997 | `8291139c6048263012959b90d2ab616f02d4ccf5a524090e4d649034c7d26b2f` |
| `bootstrap_intervals.csv` | 3118 | `632f2c2f098f43422f31070fe17504bb52401caf20e9cec26a211c949ccee6c2` |
| `cohort_audit.csv` | 1270 | `9f271e285b0edfd63d8456a6d7423d8289be8345aae0b5c0a8d5cf616aee0453` |
| `factorial_contrasts.csv` | 6018 | `13d465916d644586a7bfa607114575c29031fa4169e0c811b7183c54d49c2986` |
| `generator_vectors.parquet` | 99662 | `a6b817b9a3db207dd84f1a816a7374a5a21b186c853b73c2f2c422a1bfc830f7` |
| `nuisance_fold_assignments.parquet` | 100191 | `03dc22c344471ed09e51451fc9a79cd38b70f183565bf9034965f9b8b0e268a3` |
| `predictions.parquet` | 112752 | `7b996d697054c0144afe2f284dd3a5b83818c0b299cd804b41314931c7678bfb` |
| `pxr_boundary_decision.json` | 905 | `aa4f17bcea4f42958968b0a8b8f53223bc60bca0d06d3d4834f04bcf55ada2b0` |
| `pxr_linkage_audit.csv` | 4473735 | `d44c0b510d236e28c2b1d0afb120b4af73ed3d6a22e673ea836306650b5cbb34` |
| `pxr_overlap_summary.csv` | 2680 | `879191d9f7a9adcf08f24afb4cc3e04dc393007ceda68d58c17ade65139506fd` |
| `pxr_selection_strata.csv` | 157 | `0896b31c96259383d2fd08536c00833c6a1f08ca08eafc65bd37f4478eb48793` |
| `row_replay_seed41001.parquet` | 3666061 | `97674bcee17dcdb28f91d8ab23b87e6130dc79840a50985fad8cc6ba4507673e` |
| `semisynthetic_scenarios.csv` | 13017555 | `6283a264ed7aa519a36f01e031f7c54dcf04390b526c3071d3d098ae40ad3b12` |
| `structure_audit.json` | 4758 | `a5233c5e7d9bfe74d201e761a59c50556d323ea67ab2d1b36f0477b1e6d4b6ea` |

## Failed-code snapshot

The exact failed implementation was copied before any recovery edit:

| snapshot | SHA-256 |
|---|---|
| `artifacts/confirmatory/confirm_001/failure_001_code_snapshot/run_confirmatory_failed.py` | `9ce5f27e25e043db622b550446683ca830ec5536cc3e3f54897ed13fd4b22499` |
| `artifacts/confirmatory/confirm_001/failure_001_code_snapshot/verify_confirmatory_failed.py` | `fc26475e8b99e0c64da5ff6b1073bec82f9a59e1247bb8021223f7f9c04fbe37` |

Any recovery must use an additive amendment, new code hashes, an independently audited clearance, and a new absent output directory. The partial output above is failure evidence only and may not be used to alter scientific specifications or rescue a result.
