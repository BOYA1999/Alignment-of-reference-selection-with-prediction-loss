# Independent audit of `confirm_001` Amendment 001

Audit time: 2026-09-11T05:57:23+08:00  
Verdict: `CLEAR_TO_IMPLEMENT`

No model, prediction, mask realization, nuisance fit, bootstrap estimate, PXR propensity fit, or scientific outcome was computed in this audit. No frozen file was edited.

## Frozen identities

All byte-level SHA-256 values were independently recomputed and matched.

| Artifact | SHA-256 |
|---|---|
| Base protocol `artifacts/confirmatory/confirm_001/protocol.md` | `1819fa04713ee1a37df66c847825e5de3adcba736fba47c5c40b1c8db7796d7a` |
| Base config `artifacts/confirmatory/confirm_001/config.json` | `04fc27fc407a562853dd362d31a4c4cbc01bdb754199d5a18014cdfcbc7a5c23` |
| Base source manifest `data/external/source_freeze_20260911/source_manifest.json` | `be55555aed0b2a5bed095ae8d052500becd527c97c809d6db5b7a14e7548fd78` |
| Decision 003 `artifacts/decisions/003_go_confirm_001.md` | `0b784a409faa5858fecdd70bf77e58fbb71ab85076bf343b89bc8c7a138df9ec` |
| Protocol amendment `protocol_amendment_001.md` | `989eac8ebd6e26bdcb5cee454ec51f4c00e76ad93bf4b245b92e0a59e9ff4a49` |
| Source-manifest amendment `source_manifest_amendment_001.json` | `5b6b3337a9b7b1a4d359062ab07807527accc37b5cb63798db69e5f23837dfe5` |
| Config amendment `config_amendment_001.json` | `9008878d6a8c23aae059cf3b4aefcdc813bdd7387cbf1aa84fd6d38db620ec3a` |

The config amendment records and resolves both amendment paths and hashes. Its six-entry `formal_input_precedence` includes the base protocol, base config, base source manifest, source-manifest amendment, protocol amendment, and config amendment in that order.

## Six corrected vectors

For every row below, the amendment field path resolves in the immutable base manifest, the historical value occurs exactly once in both the base manifest and base config, the historical hash independently reproduces under its recorded terminal-LF encoding, and the corrected hash independently reproduces under Protocol Section 9.

| Vector (n) | Historical SHA-256 | Corrected SHA-256 |
|---|---|---|
| Biogen `Internal ID\tSMILES` text (3,521) | `5b31164bdda1b74e01aa7ba5c2a8da86e71471a0264db80df6fd3ff9f39b2d85` | `040e6d7526d1ebc4f0be2f955f85ec06d7f37614b0ee2c039f5d7981461939f9` |
| Biogen `cluster` text (3,521) | `c85847e68f2d726c325b0255722180201c0cca464fe0e2cbabbdb35734020a79` | `b30601402d774e5064f7024027f6cf39a66366986c2a703c67c6ae0a691bee0c` |
| Biogen `ds` text (3,521) | `8ea3ebb25b17919946dd12111f4aa9d1f8f6c619b7ff9e28b1752501bbb89a79` | `7f0fe717a49feaad22355269ed584979ef68b6a029304c74ceb09b69da4ee827` |
| PXR parent `OCNT_ID` text (10,747) | `df3e6a7832c051eeb8b98520c646541cb9207d2e28b07e5c19537f3206238402` | `7b9b757a49fbfb8d1848721aae251d34908670c0acda9c36246c71283293e168` |
| PXR availability raw `uint8` (10,747; sum 2,722) | `cd91cb263eadd8577ab4923b645eede5238c0afa906bf56160b4647c1d7b3d62` | `74efc08da49d87b9bd33f746df6ad683e2dfcdbcbc1a740c32172ce08e4a2faa` |
| PXR formal-hit raw `uint8` (10,747; sum 873) | `38958812c6cffbc216fc991caf23af9f9e8c640a96fe407ed84bb0150cc64b02` | `584db5b185e75848a5f5c9f99d5d1e12bc83b92bd7adc4f8bf8791398fe3515d` |

Exhaustiveness check: the base source manifest contains exactly these six precomputed source-derived ordered/vector hashes. The base config contains the same six plus only the unaffected bootstrap seed-index matrix hash. The amendment has six unique old hashes and six unique corrected hashes; no seventh source-derived vector requires this serialization correction.

## Descriptor, bootstrap, paths, and environment

- `sorted(Descriptors._descList, key=lambda item: item[0])` returned 217 names, all unique (`AvgIpc` through `qed`). LF joining without a terminal LF gave `84cea0bf91c3c74a704770fd5d5477dc689fc3c09e0997cd268d070acc61be73`, exactly matching the amendment.
- `default_rng(52001).integers(0,100,size=(10000,100),endpoint=False,dtype=np.int64)`, serialized as C-order little-endian `int64`, gave `0a0fa03a209b619393428cea5c45ac58e7c4dd2d2f70c573425b93dfd35ba481`, exactly matching the base config.
- All 29 frozen path-bound dependencies were present and matched their recorded hashes: base protocol and decisions (4/4), forensic inputs/results (7/7), collision/source-freeze records (4/4), Biogen/PXR staged inputs and derivative (7/7), Biogen construction-provenance files (4/4), and ExpansionRx lineage-only files (3/3). The separately listed base config and all three amendments also matched.
- The live analysis environment exactly matched the contract: Python 3.12.0; numpy 1.26.4; pandas 2.3.3; scipy 1.17.1; scikit-learn 1.8.0; RDKit 2026.03.3; XGBoost 3.3.0; pyarrow 17.0.0; Windows-11-10.0.26200-SP0.
- `artifacts/confirmatory/confirm_001/full/` was absent at the end of the audit.

## Scope and execution gate

The amendments change only hash serialization, descriptor enumeration, and per-row string-hash disambiguation. They do not override source bytes, row order, cohorts, splits, endpoints, features, models, folds, seeds, random-number scopes, estimands, formulas, gates, thresholds, output grids, decision precedence, or anti-claims. Therefore the frozen scientific design is unchanged.

`CLEAR_TO_IMPLEMENT` authorizes construction of the runner and independent verifier against the six-file precedence chain. `CLEAR_TO_RUN` remains conditional on an additive post-implementation record showing that both programs independently recompute the six corrected vectors, enforce the 217-name and bootstrap hashes, pass synthetic smoke and source-code audit, reproduce all frozen output-grid formulas, and re-confirm that the formal directory is absent immediately before launch.
