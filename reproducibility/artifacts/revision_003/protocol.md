# Incremental information after ordinary processing

Locked before the new diagnostic calculations on 2026-09-13. This is a post-result, retrospective protocol, not preregistration or independent confirmation. Parent: revision_002 B1/B2 and JPA-01–04 review; user approved substantive strengthening and emphasized information beyond ordinary preprocessing/simple methods.

## Question and scope

Does observable individual-contribution concentration provide information or a useful choice beyond ordinary preprocessing and simple comparators? Do not infer utility merely from repairing raw Ipc. Do not make an external dataset an administrative requirement or add a neural architecture for appearance.

C1 is executable now: read existing B2 estimates, with fitted low-dimensional propensity and low/drop-Ipc loss auxiliaries only. Sixteen settings, 100 existing paired masks each. No new seeds, predictors, features, targets, fits or thresholds. Original B1 low/drop under low P is supporting retrospective transport across masking mechanisms; raw/full-P cells are contextual stress evidence, not primary evidence of incremental practical value.

The current host has Ryzen 9 9950X (16 cores/32 threads), approximately 31 GiB RAM, approximately 18 GiB available at intake, and approximately 515 GiB free on D. Existing Python 3.12 has NumPy 1.26.4, pandas 2.3.3, SciPy 1.17.1, scikit-learn 1.8.0, RDKit 2026.03.3. C1 uses CSVs only, single process, no paid services/GPU. No manuscript edits until the evidence decision. All prior scientific outputs and ZIPs stay read-only.

## Prespecified diagnostics and evaluation

1. Primary observable diagnostic: top-one absolute pseudo-outcome share, already saved as top1_abs_share. It uses model auxiliaries, observation states and observed losses; no missing reference is needed to calculate it. Negative/cancelling pseudo-outcomes prevent interpreting it as a probability or uncertainty interval.
2. Simple diagnostic comparators: maximum normalized observed weight, inverse ESS fraction (rank-equivalent to negative ESS), maximum absolute SMD. Clipping is reported, with undefined constant-score associations left undefined.
3. Within each setting and auxiliary variant, report Spearman correlation with absolute risk-estimation error and AUROC for error exceeding the paired complete-case error. These evaluation outcomes use hidden full-cohort truth only for this retrospective test. Report all settings, sign counts, undefined cases and denominators. Do not pool 1,600 masks as independent chemical datasets.
4. Concrete diagnostic-use test: for each paired mask, choose the low or drop estimate with smaller top-one share; exact ties choose low. This is a deliberately simple fixed heuristic testing one form of practical utility, not a newly validated general estimator. Compare with always-low, always-drop, and complete cases. The lower-error oracle choice is a labeled unattainable reference only.
5. Primary utility contrast: equal-setting average of mean absolute risk-estimation error of the heuristic minus always-low. Negative favors the diagnostic heuristic. Also report median-error ratios, signed bias, RMSE, choice frequency and fraction of paired masks improved/worsened. Never divide by zero or silently drop a setting.
6. Conditional uncertainty: paired seed-block bootstrap, 2,000 resamples, seed 93001; the same resampled seed indices apply jointly across settings and variants. Report 2.5/97.5 percentiles for the primary utility contrast and mean AUROC difference versus each simple diagnostic. This quantifies mask Monte Carlo variability conditional on these cohorts/settings, not population or independent-dataset uncertainty. No multiple-testing-adjusted general superiority claim.

## Decision boundaries

- Evidence for useful choice requires the primary contrast to favor the heuristic with its conditional interval below zero, with full heterogeneity and safety costs disclosed. Otherwise this heuristic's utility is unsupported, not proof that every possible diagnostic is useless.
- Ranking information requires consistent direction and an advantage over the simple diagnostics within ordinary processed variants; inspect and report numerical results rather than inventing a post-result success threshold.
- If concentration only flags the already known raw-Ipc catastrophe, no incremental practical benefit is established.
- If both ranking and the fixed choice test fail to establish incremental benefit, do not tune thresholds, add heuristic variants or extend seeds to obtain a favorable result. Record the failed route and ask for a scientific route decision before constructing a materially different method.
- A positive retrospective finding is a candidate for a separately frozen independent-data test, not independent validation. AstraZeneca HLM CHEMBL3301370 is an identified possible source because it measures a matching endpoint; metadata/licence/overlap/measurement eligibility must pass before any such evaluation. Do not select among datasets by the observed diagnostic performance. No full external run is preclaimed here.

## Outputs and manuscript mapping

C1: src/revision3_incremental.py; artifacts/revision_003/incremental/{diagnostic_by_setting.csv,choice_by_setting.csv,choice_rows.csv,summary.json,manifest.json}. Supplementary quantitative evidence for the Discussion/practical-value boundary, with main-text promotion only after an explicit evidence decision. Update paper/paper_experiment_matrix.md with the actual outcome. Current publication claims and original conclusions are not silently replaced.

An independent checker must recompute the choice contrasts from saved choice rows, test no missing-truth use in the heuristic, confirm complete settings/seeds, and verify input hashes are unchanged. Failure/negative evidence remains in the record.
