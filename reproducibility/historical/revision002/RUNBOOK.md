# Execution levels and verification boundaries

Use Python 3.12 and requirements-lock.txt for scientific work. The split reconstruction has a separate requirements-split.txt environment; do not substitute it for the model-analysis environment. Run scripts without Python optimization and with -B to avoid changing a clean package inventory.

## Level 1 — integrity and saved aggregates

In a fresh extracted archive, run:

```text
python -B src/verify_revision2_release.py
```

This checks every distributed file, protected original code, permitted data fields, 448 fitted summaries, 192 alignment summaries and saved single-case arithmetic. It does not train any models or inspect excluded row-vector data.

## Level 2 — publication assets

```text
python -B src/rebuild_revision2_assets.py --output ../revision002_asset_replay
```

The output directory must not exist. The script compares every resulting file byte-for-byte against the distributed artifacts. It regenerates 16 figure files (four current main figures in three formats plus the graphical abstract in four), new tables and source data. Twelve retained supplementary figure files are copied and checked, not claimed as independently regenerated. Older build/release scripts are historical workflow records and are not the current aggregate-replay entry.

## Level 3 — minimal public-source Ipc case

```text
python -B src/reproduce_ipc_case.py
```

The pinned official CSV is retrieved in memory and its SHA-256 checked. Alternatively use --source with the checksum-matching local official CSV. Zero-based row 2945 is Mol2946. The script computes descriptors and replays 64 saved contributions; it does not refit an auxiliary or the fixed RF/XGB predictors. The source remains subject to its original terms.

## Level 4 — source retrieval and split provenance

```text
python -B src/fetch_public_inputs.py
python -B src/replay_confirmation.py --check-only
```

The retrieval entry downloads pinned source inputs and normalizes only documented line endings. For separate de novo Biogen derivative construction, follow src/rebuild_biogen_split.py --help in the requirements-split.txt environment; this reconstructs the source derivative and grouped split, not model performance. Source freeze manifests and DATA_ACCESS.md document versions and terms. Fetching inputs changes the package inventory, so Level 1 should be run before retrieval or in a separate fresh extraction.

## Level 5 — full staged scientific reconstruction

```text
python -B src/replay_confirmation.py --diagnostics
python -B src/replay_revision2.py --check-only
python -B src/replay_revision2.py
```

Use a fresh extracted package, never the author workspace. Stage 1 reconstructs the protected original scientific computations into public_replay/confirmation and creates the revision-001 forensic cache. Stage 2 requires that 19-file output and checksum-matching cache, then refits the original auxiliary grid and molecular-weight bridge, runs alignment controls and numerical checks into public_replay/revision_002. Original archived evidence is not overwritten. Stage 2 does not refit fixed predictors; those are Stage 1 inputs. The prior ZIP is not nested in this package, and Stage 2 records parent_archive_locally_verified=false instead of pretending to verify it.

**This complete Level 5 sequence has not been newly executed end to end in revision 002.** The author workspace completed the original-grid and new bridge refits against the existing fixed-predictor losses, plus the independent checks documented in verification.json. Stage 1 input preflight/static checks are historical evidence, and Stage 2 prerequisite/path checks are not a full public rerun. Do not report Level 1–3 success as Level 5 success. Differences must be reported, not repaired by overwriting the frozen evidence.
