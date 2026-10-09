# Pilot protocol — selective reference availability in ADMET model evaluation

Protocol ID: `pilot_001`

Frozen: 2026-09-10, before any project-specific predictor training or missingness-estimator result was generated.

This is a time-stamped local protocol, not a prospectively registered study. Public test labels became available before this project and published third-party model results have been read; that limitation will be disclosed.

## Study aim

Estimate how selectively available experimental references can change the estimated late-cohort error and ranking of frozen ADMET prediction models. Validate observable-data estimators in near-complete endpoints where the full analysis-cohort loss is known, then restrict genuinely sparse endpoints to assumption-indexed sensitivity analysis.

## Data contract

Source: `openadmet/openadmet-expansionrx-challenge-data`, DOI `10.57967/hf/9687`.

Pinned revision: `6b898ccc43d10d25b230fb09e22a6e30c30022b5`.

Expected files and SHA-256:

| File | Rows | SHA-256 |
|---|---:|---|
| `expansion_data_train.csv` | 5,326 | `c5214ad8c8a4d7d4d09082fcce24fe16c97a980a1e0873b91b8c8e474b79f6e4` |
| `expansion_data_test.csv` | 2,282 | `e21c0ef19795317216287e924e00ca3ef9388b7032a7fa4c7f27e2580372c85e` |
| `expansion_data_raw.csv` | 7,618 | `f674ec74cca1146bc386f832a32d4b8d921d3c312f92cb436cc005901c724a3` |

The organizer-defined train/test split is early versus late optimization. Numeric identifiers within a split are treated only as identifier order, not elapsed time.

The public package will not redistribute these raw CSV files. It will contain retrieval instructions, revision, hashes, and derived summaries.

## Cohorts and endpoints

- Training cohort: all official train molecules with a clean in-range label for the modeled endpoint.
- Target late cohort: all 2,282 official test molecules with valid structures.
- Known-truth pilot cohorts: test molecules with clean LogD or KSOL values, analyzed separately. Expected coverage is 99.5% and 95.1%, respectively.
- Empirical reference-availability masks: clean numeric availability for HLM, MLM, Caco-2 Papp A>B, MPPB, and MBPB on the same late cohort. Caco-2 efflux is omitted from the pilot because its availability pattern duplicates Papp; MGMB is omitted because only 209 test references are available.
- Raw modifiers are used to distinguish censored/out-of-range values from blanks when possible. Unresolved blanks are called unavailable references; no individual cause is assigned.

## Prediction models

Models are fit once on the official train split and frozen before missingness experiments.

1. `median`: endpoint training median.
2. `ecfp_ridge`: ridge regression, alpha 1.0, on 2,048-bit Morgan/ECFP radius 2 fingerprints with chirality.
3. `ecfp_rf`: random forest, 400 trees, `min_samples_leaf=2`, `max_features="sqrt"`, seed 20260910.
4. `ecfp_desc_xgb`: XGBoost on the same fingerprint plus finite RDKit 2D descriptors; 800 trees, depth 6, learning rate 0.03, subsample 0.8, column sample 0.8, L2 1.0, seed 20260910.

No model or hyperparameter is selected using official test outcomes. LogD is modeled on its released scale. KSOL is modeled as `log10(y + 1)`, matching the official challenge transformation apart from irrelevant constant unit shifts.

## Estimands and metrics

For model `m`, the primary target is full analysis-cohort MAE:

`R_m = mean(|Y_i - prediction_m,i|)`.

The primary pairwise target is:

`Delta_ab = mean(loss_a,i - loss_b,i)`.

Required pilot metrics:

- full-cohort MAE and normalized MAE relative to the median predictor;
- complete-case risk estimate and absolute/relative bias versus known truth;
- stabilized IPW and cross-fitted AIPW risk estimates;
- pairwise loss-difference sign and selected best model;
- empirical 95% interval coverage across repeated known-truth masks;
- propensity AUROC, clipped-weight fraction, and effective sample size.

R2 and Spearman correlation are descriptive only.

## Observation covariates and cross-fitting

Only variables known before evaluating the reference result are eligible: fixed molecular structure descriptors, ECFP fingerprint, test-to-train fingerprint similarity, and normalized numeric identifier order. Endpoint outcomes and other post-reference values are excluded.

- Actual-mask predictability is estimated by five-fold scaffold-grouped cross-fitting.
- Pilot estimators use three-fold scaffold-grouped cross-fitting.
- Propensity nuisance model: L2 logistic regression on fold-fitted standardized observation features.
- Loss nuisance model: ridge regression, alpha 10, fit only among observed rows in each training fold.
- Propensities are clipped to `[0.05, 0.95]` before weighting.
- AIPW uses fold-held-out nuisance predictions for every molecule.
- If weighted effective sample size is below 25% of observed n, the full-cohort estimate is not accepted; an overlap-population result is descriptive only.

## Known-truth stress tests

Predictors remain fixed. For LogD and KSOL separately:

1. Apply each empirical sparse-endpoint availability mask directly and compare complete-case estimates with known full-cohort truth.
2. Generate coverage-matched MCAR controls.
3. Fit cross-fitted observation models for HLM, MLM, and Caco-2 Papp. Their out-of-fold probabilities, clipped to `[0.05, 0.95]`, define three empirical MAR anchor surfaces.
4. For each anchor, generate 40 Bernoulli masks with seeds 31001–31040.
5. Generate two MNAR stress variants by shifting each anchor logit by `gamma * z(Y)`, with `gamma=-1` and `gamma=+1`, and recentering the intercept to preserve the anchor mean coverage.

MCAR quantifies sampling variance; MAR tests recoverability under observed-covariate selection; MNAR demonstrates non-identification when reference availability depends on the unobserved outcome.

The primary estimator comparison is complete case versus AIPW for `ecfp_rf` and `ecfp_desc_xgb`. Other model results are supporting diagnostics.

## Uncertainty

- Molecules are observational units; Bemis–Murcko scaffolds are resampling/variance clusters.
- Model losses are paired on the same molecules.
- AIPW normal intervals use cluster-robust influence-function standard errors.
- Mask repetitions are Monte Carlo repetitions and are not treated as independent chemical datasets.
- Real-mask point comparisons use paired scaffold bootstrap with 1,000 replicates in the full stage; the pilot may use 200 to verify wiring.

## Pilot Go/Stop rules

All gates are evaluated as written:

1. Parent cohort, structures, phase, and observation/censoring indicators reconstruct for at least 7,608 ML-ready molecules.
2. Late-cohort LogD and KSOL evaluable coverage is at least 90%.
3. At least two non-dummy models have normalized MAE below 1.0 on both truth endpoints.
4. At least two empirical sparse masks have grouped cross-fitted propensity AUROC at least 0.65 and weighted ESS at least 25% of observed n after clipping.
5. At least one truth-endpoint × empirical/empirical-MAR mask produces complete-case absolute bias at least 5% of full-cohort risk or changes the best-model identity.
6. Across MAR stress tests, AIPW reduces median absolute bias by at least 25% versus complete case and its empirical 95% interval coverage is at least 85%.
7. No train/test structural overlap, outcome leakage, non-finite required metric, or source/version mismatch remains.

If gates 1–3 or 7 fail: stop the route. If gates 4–6 fail: retain the negative pilot and do not expand. Ambiguous causes of real missingness do not invalidate the known-truth experiment but limit real sparse endpoints to estimand and sensitivity reporting.

## Prohibited interpretations

- Do not claim that IPW/AIPW identifies real MNAR endpoint performance.
- Do not label unavailable references as inactive, negative, low, failed, or unmeasured without source evidence.
- Do not infer a causal assay cascade.
- Do not call this independent external or prospective experimental validation.
- Do not generalize one public campaign to all ADMET benchmarks.

## Amendments

Any correction after this freeze must be appended here with time, reason, affected outputs, and whether it changes a gate. Silent edits are prohibited.

### A01 — 2026-09-10, raw-file hash transcription correction

The initially recorded raw-file hash had 63 hexadecimal characters and was therefore not a valid SHA-256 value. A fresh download from pinned revision `6b898ccc43d10d25b230fb09e22a6e30c30022b5` produced a 686,791-byte file with SHA-256 `f674ec74cca1146bc386f832a32d4b8d921d3c312f92cb436cc005901c724a3c`; the repository API independently reports the same file size and Git blob identity `112a951e8774614ff6a6fc7486b430648b973a35`. The final `c` was missing from the audit handoff. This amendment changes only the expected raw-file hash and does not alter any scientific gate, cohort, model, metric, or analysis rule.

### A02 — 2026-09-10, descriptor and gate-6 aggregation clarification

Before any project-specific predictor was trained, implementation preflight found that “finite RDKit 2D descriptors” and the gate-6 aggregation were not machine-explicit. Under RDKit `2026.03.3`, all 217 names in `rdkit.Chem.Descriptors._descList` produced finite values for every one of the 7,608 ML-ready structures. The predictor and observation-feature matrices therefore use all 217 registry descriptors, sorted by descriptor name; the run manifest must record the complete name list, count, and RDKit version. No descriptor was selected using an endpoint value or result.

Gate 6 is computed over the 480 equally weighted primary MAR rows (2 truth endpoints × 3 empirical MAR anchors × 40 seeds × 2 primary models). Bias reduction is `1 - median(abs(AIPW - truth)) / median(abs(CC - truth))`, with each median taken over the same 480 rows. Coverage is the proportion of those 480 AIPW 95% intervals containing the corresponding full-cohort truth. This clarification changes no threshold, cohort, model class, stress mechanism, or success rule; no result existed when it was made.

### A03 — 2026-09-10, propensity and estimator formula clarification

Before any project-specific outcome predictor or pilot estimator was run, the nuisance-model defaults and the phrase “stabilized IPW” were made machine-explicit. Every propensity fold uses `LogisticRegression(penalty="l2", C=1.0, solver="liblinear", class_weight=None, max_iter=1000)` after fitting the feature standardizer on that training fold only. No propensity hyperparameter is tuned.

For clipped out-of-fold propensity `p_i`, observation indicator `R_i`, loss `L_i`, and held-out loss-nuisance prediction `m_i`, the stabilized IPW estimate is the Hájek ratio `sum(R_i L_i / p_i) / sum(R_i / p_i)`. The AIPW estimate is `mean(m_i + R_i (L_i - m_i) / p_i)`. Weighted effective sample size is `(sum(R_i / p_i))^2 / sum((R_i / p_i)^2)`, and its gate-4 fraction uses the observed count `sum(R_i)` as denominator. For AIPW, define centered influence values `phi_i = m_i + R_i (L_i - m_i) / p_i - AIPW`; with `G` scaffold clusters and cluster sums `S_g = sum(phi_i in g)`, the standard error is `sqrt((G / (G - 1)) * sum(S_g^2)) / n`, followed by a two-sided normal 95% interval. Acyclic molecules remain singleton structure groups. These formulas only resolve implementation ambiguity and change no gate, threshold, cohort, feature, or scientific claim.

### A04 — 2026-09-10, gate-5 setting and repetition aggregation clarification

Before any project-specific outcome predictor or pilot estimator was run, gate 5 was made immune to a chance pass from one of the 40 Monte Carlo repetitions. Gate 5 uses only the two prespecified primary models, `ecfp_rf` and `ecfp_desc_xgb`. For each of the 10 actual-mask settings (2 truth endpoints × 5 empirical masks), a setting triggers if either primary model has `abs(CC - truth) / truth >= 0.05`, or if the complete-case winner differs from the full-truth winner. For each of the six empirical-MAR settings (2 truth endpoints × 3 anchors), the 40 seeds are repetitions of one setting: its bias branch triggers when `max_model median_seed(abs(CC - truth) / truth) >= 0.05`, and its ranking branch triggers only when at least 21 of 40 complete-case winners differ from the full-truth winner. Gate 5 passes if any actual or MAR setting triggers. MCAR and MNAR do not enter gate 5. If the two full-truth risks are exactly tied, the ranking branch for that truth endpoint is non-adjudicative; an exact complete-case tie is not counted as an identity change. This clarification preserves the frozen 5% threshold and “at least one setting” rule while preventing seed-wise multiple-opportunity passing.

### A05 — 2026-09-10, mask-generation and gate-integrity clarification

Before any project-specific outcome predictor or pilot estimator was run, mask generation was made fully deterministic apart from the frozen Bernoulli seeds. Gate-4 AUROC is computed from un-clipped five-fold scaffold-grouped out-of-fold probabilities on all 2,282 valid late-cohort structures; clipping to `[0.05, 0.95]` is applied only for weights and empirical MAR anchor probabilities. For a given truth cohort, each MCAR control uses the corresponding empirical mask's observed fraction within that truth cohort as its constant Bernoulli probability. Each empirical MAR mask is an independent Bernoulli draw from the clipped anchor probability restricted to that truth cohort; simulated counts are not forced to equal a target count.

For MNAR generation, `z(Y)` is calculated on the truth endpoint's analysis scale (released LogD; `log10(KSOL + 1)`) using that truth cohort's population mean and population standard deviation (`ddof=0`). For every truth-endpoint × anchor × gamma setting, a scalar intercept shift is solved so that the mean MNAR probability equals the mean clipped MAR anchor probability in that truth cohort. The resulting expit probability is not clipped again. `gamma=+1` makes larger analysis-scale outcomes more likely to have a reference; `gamma=-1` makes them less likely. Outcomes enter only this MNAR generator and never the observation features or fitted nuisance estimators. Every actual or simulated mask used for IPW/AIPW is re-estimated with the frozen three-fold scaffold-grouped nuisance procedure; generator probabilities are diagnostics, not estimator inputs.

Gate 6 retains all 480 prespecified MAR model rows. It is not filtered by endpoint, anchor, model, estimated ESS, interval width, or result direction. A missing or non-finite row, zero complete-case median bias denominator, or a row that cannot produce the prespecified nuisance estimates makes gate 6 non-passing and is surfaced to gate 7; it is never silently dropped. AIPW interval width is reported alongside coverage. These clarifications change no threshold, cohort, feature, or claim.

### A06 — 2026-09-10, loss-nuisance scaling and ESS acceptance clarification

Before any project-specific outcome predictor or pilot estimator was run, the loss-nuisance implementation and the existing ESS acceptance rule were made machine-explicit. Within each three-fold nuisance split, one `StandardScaler` is fit on all rows of the nuisance-training fold and applied to its training and held-out rows. The propensity model uses all standardized training rows; the loss-nuisance `Ridge(alpha=10, solver="lsqr", tol=1e-6)` uses only the observed rows from that same standardized training fold. Held-out loss predictions are not truncated, and AIPW estimates are not forced to be nonnegative.

The protocol's prespecified statement that a full-cohort estimate is not accepted when weighted ESS is below 25% of observed n applies to every primary-model MAR row. All 480 rows remain in the output and in the descriptive bias and coverage arithmetic, but gate 6 can pass only if all 480 rows also satisfy `ESS / observed_n >= 0.25`. An ESS failure is recorded separately from a numerical failure and is not relabeled as gate 7 unless it reflects a missing, non-finite, leakage, or version-integrity problem. This clarification changes neither the 25% threshold nor the 480-row aggregation.

### A07 — 2026-09-11, deterministic nuisance fit and gate-grid integrity clarification

Static code review, still before any project-specific model or estimator was run, found that scikit-learn's `liblinear` interface accepts a random state and that row-count checks alone cannot prove the gate inputs are the intended settings. Every propensity fit therefore uses `random_state=20260910` in addition to the A03 parameters. Gate 5 must verify the unique Cartesian products for the 10 actual settings and six MAR settings with all 40 frozen seeds and both primary models; gate 6 must verify the unique 2 × 3 × 40 × 2 MAR Cartesian product, not merely a total of 480 rows. Any missing, duplicated, unexpected, or non-finite key/value makes the relevant gate non-adjudicative and non-passing. This is implementation hardening only; it changes no analysis setting or threshold.
