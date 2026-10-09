# CILS structure-split robustness addendum

This addendum contains the prespecified ten-partition robustness analysis added to the manuscript *Alignment of reference selection with prediction loss and descriptor extrapolation in clearance predictor evaluation*.

## Scientific scope

The analysis repeats the original structure-grouped splitting rule with seeds 61,453–61,462, refits the two fixed predictor configurations and the selected nuisance-model grid within every partition, and evaluates 100 masks per endpoint–anchor–mechanism–predictor cell. It is a same-dataset robustness analysis of one public Biogen collection. It is not an external validation cohort.

The retained result is conditional: the raw-Ipc failure recurred in all 32 eligible MAR cells from the four partitions placing Mol2946 in the target cohort, but in none of 48 corresponding cells from the six partitions retaining that compound in prediction training. XGBoost had lower finite-cohort MAE than the random forest in all 20 partition–endpoint comparisons.

## Contents

- `protocol/`: frozen protocol and split-screen records.
- `code/`: split generation, complete analysis, asset generation and independent verification code.
- `results/`: the 80,000 setting-level risk estimates, ten per-split estimate files, aggregate summaries, decision record and verification records. Row-level prediction files are intentionally excluded.
- `paper_assets/`: figure source data, Table S17 source tables and publication figure exports.
- `manifest-sha256.txt`: SHA-256 for every payload file other than the manifest itself.

## Replay boundary

The analysis code expects the pinned public Biogen source derivative and the earlier frozen confirmatory assets described by their hashes in `results/manifest.json`. Those source-bearing inputs are not redistributed here. With those inputs restored to the paths stated in the code, run `revision4_split_vectors.py`, then `revision4_multisplit.py`, and finally `verify_revision4_multisplit.py`. The supplied results permit independent checking of the retained aggregate decisions without redistributing molecular structures.

## Verification

Run `python verify_package.py` from this directory to check the package manifest, row counts, the ten-file aggregation, grid completeness, decision statistics and manuscript-facing asset hashes. A successful run prints `PASS_CILS_MULTISPLIT_ADDENDUM`.
