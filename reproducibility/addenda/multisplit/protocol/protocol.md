# Revision 004: prespecified structure-split perturbation

**Frozen before outcome analysis:** 2026-09-18  
**Role:** targeted robustness follow-up to the CILS manuscript  
**Status:** READY_FOR_OUTCOME_BLIND_SCREEN

## Question

Does the paper's central descriptor-to-estimator failure recur under alternative, fully structure-grouped Biogen partitions, or is its magnitude mainly a consequence of the single derived split used in the original study?

This is the highest-priority feasible gap because the current manuscript explicitly rests on one derived train/test split. It is not an external validation and will not be described as one.

## Frozen inputs and invariants

- The pinned 3,521-row Biogen source derivative, its 1,905 BitBIRCH-Lean cluster labels, endpoint values, molecular structures, feature definitions, RF400/XGB800 configurations, nuisance models, clipping rule, three-fold scaffold cross-fitting, anchors and mask seeds remain unchanged.
- Each perturbation holds out complete existing clusters until the nominal 30% target is reached, using the exact upstream `split_by_cluster` algorithm.
- Ten split seeds are fixed before screening: **61453 through 61462 inclusive**. Seed 61453 is the original split. Every seed will be retained; no seed may be substituted or selected by its result.
- Mask seeds remain 41001 through 41100 and are paired across loss specifications within a split.
- The two sentinel endpoints (`LOG_HLM`, `LOG_RLM`), two regular anchors (`LOG_SOL`, `LOG_MDR1_ER`), two mechanisms (MCAR, MAR) and two fixed predictors are all retained.

## Staged execution

1. **Outcome-blind screen.** Recreate the ten cluster partitions and record split size, endpoint/anchor availability counts, whether source row 2945 is assigned to the target split, and descriptor ranges. No predictor or risk outcome is fitted at this stage. All ten partitions proceed unless a mechanical requirement fails.
2. **Smoke run.** Run split 61453 with mask seed 41001 and verify it against the frozen original predictions and revision-002 risk estimates at scalar tolerance 1e-10 where deterministic equality is expected.
3. **Full run.** Refit RF400 and XGB800 for every split, recreate endpoint-specific target cohorts, refit low-dimensional propensity and three conditional-loss models for every mask: low-dimensional, raw full-dimensional, and full-dimensional with raw Ipc removed. Retain negative and null results.

## Primary summaries

For each split x endpoint x anchor x mechanism x predictor cell, summarize 100 masks by median absolute finite-cohort risk error for complete case, Hajek, and AIPW. The fitted-propensity AIPW comparisons of interest are:

1. raw full-dimensional loss versus dropped-Ipc full-dimensional loss;
2. raw full-dimensional loss versus low-dimensional loss;
3. the identity and absolute share of the largest pseudo-outcome contributor;
4. RF400 versus XGB800 finite-cohort MAE and the stability of their ordering across splits.

MAR is the primary selection mechanism; MCAR is a prespecified sensitivity analysis. No significance test across ten deterministic partitions is planned. Results are descriptive partition-robustness evidence.

## Interpretation rules

- **Conditional mechanism recurrence:** among partitions where source row 2945 belongs to the finite endpoint target cohort, raw-loss median error must exceed both processed-loss alternatives by at least 10-fold in at least 80% of fixed MAR cells, and row 2945 must be the largest absolute contributor in at least 95% of raw-loss estimates.
- **Partition sensitivity:** if raw-loss median error differs by at least 10-fold between partitions with versus without row 2945 in the target cohort, the manuscript must explicitly describe the catastrophic scale as partition-dependent.
- **Broader raw-Ipc instability:** this may be claimed only if raw-loss error exceeds dropped-Ipc error by at least 10-fold in at least half of MAR cells from partitions where row 2945 is absent from the target cohort.
- Failure of any rule is retained as a negative result and narrows, rather than rescues, the claim.

## Stop condition

Stop after the ten frozen partitions, verification, aggregate reporting and manuscript mapping. Do not add alternative seeds, models, descriptors, thresholds or datasets in response to the observed result. Independent external validation remains a separate future-data requirement.
