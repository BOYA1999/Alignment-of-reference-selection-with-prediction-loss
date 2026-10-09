# Alignment of reference selection with prediction loss

Analysis code and saved results for *Alignment of reference selection with prediction loss and descriptor extrapolation in clearance predictor evaluation*.

The study examines how selective reference availability and descriptor extrapolation affect estimates of clearance prediction error. It compares complete-case, Hájek and augmented inverse-probability weighted (AIPW) estimates, including fitted auxiliary models and different treatments of the RDKit Ipc descriptor (information content based on the characteristic polynomial).

The results depend on the fitted configuration and target cohort. AIPW did not consistently improve on fitted Hájek. Across ten structure-grouped partitions of the same Biogen collection, the raw-Ipc failure recurred in all 32 eligible MAR cells containing Mol2946 in the target cohort, and in none of the 48 corresponding cells where that compound stayed in prediction training. These partitions assess robustness within one dataset; they do not provide an external validation cohort.

## Run the saved-result checks

Use Python 3.12 with assertions enabled. From the repository root:

```bash
python -m pip install -r requirements-check.txt
python -B reproducibility/verify_minor_revision.py
```

This command checks the file inventory and SHA-256 hashes, then recalculates summaries and conditional Monte Carlo intervals from the saved estimates. It checks 448 fitted summaries, 64 case contributions, 3,200 fixed-rule choices, 64 paired ratio intervals, 80,000 multisplit records and the 32 Figure 3 inputs. It reads local files and prints two JSON reports. It does not download assay data or fit models.

The expected status strings are `PASS_LOCAL_CANDIDATE_SAVED_RESULT_QA` and `PASS_MINOR_REVISION_SAVED_RESULTS_AND_FIGURE_INPUTS`. These names come from the original checkers and describe the saved-result checks. The reported `public_access_created: false` means that the checker itself creates no public access.

## Find the code and results

| Path | Contents |
| --- | --- |
| `reproducibility/verify_minor_revision.py` | Entry point for the current saved-result checks. |
| `reproducibility/verify_candidate.py` | Inventory, hash and aggregate checks used by the entry point. |
| `reproducibility/src/` | Revision-003 analysis and row-level verification scripts. |
| `reproducibility/artifacts/` | Saved revision-002 estimates and revision-003 choices and diagnostics. |
| `reproducibility/addenda/multisplit/` | Ten-partition analysis code, protocols, saved estimates and checks. |
| `reproducibility/historical/revision002/src/` | Earlier estimators, refitting scripts, input retrieval and staged reconstruction. |
| `reproducibility/paper/` | Figure exports, electronic tables and plotting inputs under `figures/`, `tables/` and `source_data/`. |
| [docs/REPRODUCING.md](docs/REPRODUCING.md) | Figure generation, input reconstruction and the limits of each check. |
| [docs/RELEASE_SCOPE.md](docs/RELEASE_SCOPE.md) | Source archive, publication inventory and licensing scope. |

## Scope

The repository contains analysis code, protocols, selected historical records and saved results. Raw assays, molecular structures, full compound-level predictions, row-vector caches, manuscript text, author Word files and submission correspondence are excluded. Input URLs and checksums identify the sources needed for further reconstruction.

The current check does not repeat the historical alignment Parquet arithmetic or the original S18 AUROC bootstrap contrasts. Full predictor and auxiliary-model reconstruction requires the source inputs and separate scientific environments described in the run instructions; that complete sequence was not rerun for this code publication.

The existing [MIT license](LICENSE) is retained for project software. Third-party code keeps its own notices, and source datasets remain subject to their original terms. See [the release scope](docs/RELEASE_SCOPE.md) before redistributing data.
