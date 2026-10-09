# Reproducing the results

## Check the files and saved estimates

From the repository root, use Python 3.12 and the versions in `requirements-check.txt`:

```bash
python -m pip install -r requirements-check.txt
python -B reproducibility/verify_minor_revision.py
```

Keep assertions enabled; `python -O` disables checks. The `-B` option prevents bytecode files from changing the checked inventory. Run this command before any data retrieval or figure generation, since those steps write files.

The check covers:

- Exact paths, sizes and SHA-256 hashes of the reproducibility payload, including the two subdirectory manifests.
- Reaggregation of 448 fitted summaries and 64 saved case contributions.
- The 3,200 fixed-rule choices, 128 setting summaries and S18 utility changes and conditional intervals.
- The 64 S19 median-error ratios and paired conditional intervals, including the 3/16 and 1/16 point comparisons with fitted Hájek.
- Agreement between the ten per-partition files and all 80,000 multisplit estimates, including the conditional 32/32 and 0/48 recurrence counts.
- The 32 current Figure 3 inputs, the five-estimate grid per multisplit unit and the retained PXR summary rows.

The repeated masks and bootstrap intervals describe Monte Carlo precision conditional on the saved fits and cohorts. The checker does not recalculate the historical alignment Parquet's 192 summaries or the original S18 AUROC bootstrap contrasts. Their saved files remain available with the earlier verification records.

## Generate Figure 3

Use a separate working copy so the checked exports remain available for comparison. The plotting versions are recorded in `requirements-plot.txt`; install them in a separate environment from the saved-result check if needed.

```bash
python -m pip install -r requirements-plot.txt
python -B reproducibility/paper/figures/scripts/Figure_3_fitted_hajek.py
```

The script reads `reproducibility/paper/source_data/Figure_3_fitted_hajek.csv` and writes PNG, PDF, SVG and TIFF exports. Binary figure files can differ across dependency versions, so regeneration can make the frozen hash check fail even when the numerical inputs agree. The earlier oracle figure is retained as Figure S5 with its original image bytes.

## Reconstruct inputs or refit models

Start with the pinned sources in `reproducibility/historical/revision002/DATA_ACCESS.md` and the staged commands in `reproducibility/historical/revision002/RUNBOOK.md`. These are dated records of the original analysis. The software, split-reconstruction and document environments have separate requirements files; the small root requirements files cover only the checks and Figure 3 above.

| Task | Inputs and instructions |
| --- | --- |
| Recalculate the Mol2946 descriptor and saved contributions | Run `src/reproduce_ipc_case.py` from the historical revision directory. It retrieves the pinned official Biogen CSV in memory or accepts a checksum-matching local source. It replays 64 contributions without refitting predictors or auxiliaries. |
| Retrieve public inputs and reconstruct the grouped split | Use the historical fetch and split-reconstruction scripts with their pinned source revisions and separate split environment. |
| Reconstruct the original confirmation and revision-002 analyses | Follow the historical staged replay entries in a separate reconstruction tree. They write to `public_replay/` and require the documented inputs. The complete staged sequence was not rerun for this publication. |
| Repeat revision-003 row-level verification | Restore confirmation predictions and all 1,600 revision-002 row-vector NPZs at the script's expected paths. Those files are excluded here. |
| Refit all ten partitions | Restore the pinned Biogen derivative, historical confirmation inputs and original project layout expected by `addenda/multisplit/code/`. Its README and manifests describe the prerequisites. |

The historical standalone release checker expects the complete earlier candidate, including manuscript files omitted from this code repository. Use the root saved-result entry for the distributed subset. The historical manuscript and package builders also refer to files from the author workspace and are not current repository entry points. Restore the documented prerequisites in another tree before running generation scripts; do not remove saved evidence to satisfy an empty-output-directory check.

PXR files support a separate selection-record analysis. They do not add external clearance-risk validation. Source retrieval and full model fitting were not rerun as part of the language edits and code publication.
