# Reproducibility files

These files come from the CILS formatting revision of 9 October 2026. The repository publishes its analysis-code subset, with manuscript and submission material omitted as listed in `PUBLICATION_SCOPE.json`. The original local archive remains unchanged. The repository README and `docs/REPRODUCING.md`, one directory above this payload, give the current entry points.

The source archive derived from `CILS_language_polished_candidate_20261009.zip` (SHA-256 `4a3982408d8e649fc616dec5e5e023ebfb0a5e8b9b004a43a149e9816018171f`). That editorial revision synchronized equation numbering, Ipc definitions and end-matter order in the separately delivered manuscript. Word table formatting was handled outside the archive. This code subset retains the analysis scripts, estimates, electronic tables and figure files from that revision.

## Check saved results

With Python 3.12, NumPy and pandas installed, run from the repository root:

```text
python -B reproducibility/verify_minor_revision.py
```

The command first runs `verify_candidate.py`, then checks the current figure inputs, the five-record multisplit grid, the relocated PXR rows and the retained oracle figure. `MANIFEST.json` lists the size and SHA-256 of each payload file. Inventory and hashes are checked again after the calculations. The scripts read local files without accessing the network or fitting models.

The checks reaggregate 448 fitted summaries and 64 case contributions, verify the 80,000 multisplit records and 3,200 fixed choices, recalculate S18 utility intervals and all 64 S19 paired ratio intervals, and check the retained negative results. Repeated masks and bootstrap intervals describe conditional Monte Carlo precision. The default entry does not recalculate the historical alignment Parquet arithmetic for 192 summaries or the original S18 AUROC bootstrap contrasts; those files and earlier records are retained.

## Generate Figure 3

The current main Figure 3 shows 32 AIPW-to-fitted-Hajek ratios and paired conditional intervals from Table S19. Its generating script, plotting style, full-precision input CSV and four exports are included. From the repository root, with Matplotlib, NumPy and pandas installed:

```text
python -B reproducibility/paper/figures/scripts/Figure_3_fitted_hajek.py
```

The script reads `paper/source_data/Figure_3_fitted_hajek.csv` within this payload and writes PNG, PDF, SVG and TIFF files. Check a clean copy first, then generate figures in another working copy: dependency versions can change binary exports and therefore their hashes. The recorded plotting environment used Python 3.12, Matplotlib 3.10.8, pandas 2.3.3 and NumPy 1.26.4. The saved-result check uses its separately recorded environment.

The former main Figure 3 is retained as supplementary Figure S5 with unchanged image bytes. The current main Table 1 covers clearance; its two former PXR rows are also supplied in Table S5C. The original electronic Table 1 is retained as a source record.

## Reconstruct inputs and fits

Raw assays, structures and full compound-level predictions are excluded. `historical/revision002/DATA_ACCESS.md` records source access and pinned identities; its `RUNBOOK.md` describes the staged reconstruction. Source terms still apply when retrieving or redistributing inputs. Work in a separate reconstruction tree and retain the frozen records. Input reconstruction was not rerun for the formatting revision or this code publication.

Full refitting also requires the scientific environments, original project layout and omitted inputs. Revision-003 row-level verification needs row-vector NPZs and reconstructed predictions. Multisplit fitting needs its split and feature inputs. The relevant code and protocols remain under `historical/revision002/` and `addenda/multisplit/`. The complete predictor and auxiliary reconstruction sequence was not rerun for this publication.

## Distribution and historical records

Raw assay copies, SMILES/scaffolds, full compound-level predictions and vectors, PXR linkage records, feature caches, row-vector NPZs, user/reviewer attachments, author Word files and local environments remain excluded. The permitted Mol2946 locator, limited case arithmetic and split/source-row indices are retained. PXR supports a separate selection-record analysis, with no external clearance-risk validation.

Historical protocols, licenses, status files and graphical abstracts describe their dated analyses and deliveries. Some historical paths refer to manuscript or submission files omitted from this code subset. The regenerated inventories cover the files distributed here; original execution records continue to identify the executed source bytes. No DOI or submission is created by the saved-result checks or this publication. The legacy graphical abstract remains historical material rather than a synchronized current submission figure.
