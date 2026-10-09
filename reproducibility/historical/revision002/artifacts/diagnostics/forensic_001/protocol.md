# Frozen protocol — `forensic_001`

Frozen: 2026-09-11, after `pilot_001 = STOP_PILOT_NEGATIVE` and before any `forensic_001` result is generated.

Amendment 1: 2026-09-11, still before any `forensic_001` result. A pre-result audit removed ambiguities in oracle probabilities, shared folds/masks, row-level replay, signature aggregation, and the F1 source keys. No threshold or scientific direction was changed.

## Status and purpose

This is an exploratory failure-mechanism diagnostic, not a retry or amendment of `pilot_001`. Its only purpose is to distinguish:

1. an implementation/formula mismatch;
2. instability from the estimated propensity model;
3. instability from the loss-augmentation model;
4. instability associated with the full high-dimensional molecular observation representation.

The `pilot_001` gate decision remains `STOP_PILOT_NEGATIVE` regardless of this diagnostic.

## Frozen inputs

- Use `config/pilot.json` at SHA-256 `469f9b1a57303a8b0bf3c6212bb433f320da774f09750a59152463961e4a557b`.
- Pinned source hashes are train `c5214ad8c8a4d7d4d09082fcce24fe16c97a980a1e0873b91b8c8e474b79f6e4`, test `e21c0ef19795317216287e924e00ca3ef9388b7032a7fa4c7f27e2580372c85e`, and raw `f674ec74cca1146bc386f832a32d4b8d921d3c312f92cb436cc005901c724a3c`.
- The frozen RDKit version is `2026.03.3`.
- Use the frozen late-cohort predictions from `artifacts/experiment/pilot_001/full/predictions/predictions.parquet`.
- Truth endpoints: `LogD`, `KSOL`.
- Prediction models: `ecfp_rf`, `ecfp_desc_xgb`.
- Availability anchor: `HLM`, chosen because it is the first prespecified MAR anchor, not because of its result.
- Seeds: all original 40 mask seeds, `31001` through `31040`.
- Mechanisms:
  - `mcar_hlm`: along the fixed ascending `row_index` truth-row order, set `p_oracle = mean(row_audit.available_HLM[truth_idx])`, then generate exactly once per endpoint and seed as `default_rng(seed).random(n) < p_oracle`;
  - `mar_hlm`: along the same fixed truth-row order, set `p_oracle_i = row_audit.propensity_clipped_HLM[truth_idx]`, then generate exactly once per endpoint and seed as `default_rng(seed).random(n) < p_oracle_i`.

Neither oracle probability is re-estimated from `lowdim_12` or clipped again. Within an endpoint x mechanism x seed cell, the generated `R` vector is reused unchanged across both feature sets and both prediction models.

## Frozen feature sets

`full_2267`:

- 2,048-bit ECFP, radius 2, chirality on;
- all 217 sorted finite RDKit 2D descriptors;
- maximum Tanimoto similarity to official train;
- normalized numeric identifier order.

`lowdim_12`:

- `MolWt`, `MolLogP`, `TPSA`, `NumHDonors`, `NumHAcceptors`, `NumRotatableBonds`, `RingCount`, `FractionCSP3`, `HeavyAtomCount`, `MolMR`;
- maximum Tanimoto similarity to official train;
- normalized numeric identifier order.

The low-dimensional block is fixed before its results are known. It is a diagnostic contrast, not a newly optimized estimator.

## Frozen nuisance fitting

- Three-fold scaffold-grouped out-of-fold fitting within each truth cohort.
- Generate one endpoint-specific fold assignment before any scenario fitting. Reuse those exact row-level fold IDs across feature sets, mechanisms, seeds, and prediction models; save the fold ID with every replay row and save a SHA-256 of the ordered fold-ID vector for each endpoint.
- Fold-specific `StandardScaler`, fit on all nuisance-training rows.
- Propensity: L2 logistic regression, `C=1`, `liblinear`, no class weighting, `max_iter=1000`, `random_state=20260910`.
- Estimated probabilities clipped to `[0.05, 0.95]`.
- Loss nuisance: ridge `alpha=10`, solver `lsqr`, tolerance `1e-6`, fit only on observed training rows; predictions are not truncated.
- The shared endpoint-specific fold IDs are used for both propensity and loss nuisance fits.

## Frozen estimators and decomposition

For each endpoint x mechanism x seed x feature set x model, record:

- complete full-cohort MAE truth;
- complete-case mean loss;
- outcome-regression mean `mean(mhat)`;
- estimated-propensity Hájek IPW;
- estimated-propensity AIPW;
- oracle-propensity Hájek IPW;
- oracle-propensity AIPW using the same fitted `mhat`;
- fitted-propensity clipping fraction and ESS/observed;
- oracle-propensity ESS/observed;
- scaffold-cluster influence-function SE and interval for both AIPW versions.

Oracle probabilities are the uniquely defined generator probabilities above. They are used only because the masks are synthetic and their generator is known.

For seed `31001`, save a long-format table covering all 16 endpoint x mechanism x feature-set x prediction-model cells. Expected rows are `2,270 x 8 + 2,170 x 8 = 35,520`, with unique key `{endpoint, mechanism, seed, feature_set, model, row_index}`. Each row must contain `molecule_id`, `scaffold_group`, `scaffold_sha256`, `fold_id`, `R`, `loss`, `mhat`, `p_fitted`, `p_oracle`, and fitted/oracle versions of `augmentation`, `pseudo_outcome`, and centered `influence`.

For `v in {fitted, oracle}`, freeze the row formulas as

`augmentation_i^v = R_i * (loss_i - mhat_i) / p_i^v`,

`Q_i^v = mhat_i + augmentation_i^v`,

`phi_i^v = Q_i^v - mean(Q^v)`.

With `N` rows and `G` distinct scaffold groups, freeze the cluster standard error as

`SE^v = sqrt((G / (G - 1)) * sum_g (sum_{i in g} phi_i^v)^2) / N`.

For `v in {fitted, oracle}`, freeze the Hájek estimator as

`IPW^v = sum_i [R_i * loss_i / p_i^v] / sum_i [R_i / p_i^v]`.

## Fixed integrity checks

F1 — reproduction: every `full_2267` CC, fitted IPW, fitted AIPW, SE, ESS fraction, and clipping fraction must reproduce the corresponding 320 `pilot_001` HLM MCAR/MAR rows within absolute tolerance `1e-12`. The source keys are exactly `mechanism=mcar, mask=HLM, anchor=null` for 160 rows and `mechanism=mar, mask=null, anchor=HLM` for 160 rows. Both source and diagnostic tables must have the exact expected key grid and no duplicate keys. Also compare `n_total`, `n_observed`, `full_truth`, and `generator_mean_probability`.

F2 — row replay: the saved seed-31001 table must contain exactly the 16 cells and 35,520 unique rows specified above. All fitted/oracle IPW, AIPW, and cluster-SE values must recompute from the frozen row formulas within `1e-12`; fold and scaffold keys/hashes must also match their frozen ordered vectors.

Failure of F1 or F2 reclassifies the mechanism diagnostic as `IMPLEMENTATION_MISMATCH`; it does not create a positive scientific result.

## Prespecified descriptive mechanism signatures

All comparisons use pooled median absolute bias across both endpoints, both prediction models, and all 40 seeds, evaluated separately for MCAR and MAR.

- `PROPENSITY_INSTABILITY_SUPPORTED`: using only `full_2267`, median oracle-propensity AIPW absolute bias is at most 50% of median fitted-propensity AIPW absolute bias in both mechanisms.
- `AUGMENTATION_INSTABILITY_SUPPORTED`: using only `full_2267`, median oracle-propensity AIPW absolute bias is greater than twice median oracle-propensity IPW absolute bias in both mechanisms.
- `HIGH_DIMENSION_INSTABILITY_SUPPORTED`: in both mechanisms, median `lowdim_12` fitted-propensity AIPW absolute bias is at most 50% of the corresponding `full_2267` median, and `median(clip_fraction_full_2267) - median(clip_fraction_lowdim_12) >= 0.20`.

If F1/F2 pass and at least one signature is met, the diagnostic decision is `MECHANISM_SIGNAL_DETECTED`. Otherwise it is `MECHANISM_INDETERMINATE`.

These thresholds classify evidence for a new hypothesis only. They cannot change `pilot_001`, select a favorable final estimator, or authorize confirmatory expansion without a separate live collision audit and decision record.

## Required outputs

- `scenario_decomposition.csv`: exact 640-row diagnostic grid.
- `row_replay_seed31001.parquet`: row-level formula inputs and outputs for the fixed replay seed.
- `summary.json`: F1/F2 results, pooled mechanism summaries, signature flags, and diagnostic decision.
- `run_manifest.json`: source/input/output hashes, versions, command, and runtime.
- `verification.md`: independent replay scope and result.

Raw source CSV files must not be copied into a public artifact.
