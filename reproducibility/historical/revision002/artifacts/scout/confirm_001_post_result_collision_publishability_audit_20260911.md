# `confirm_001` post-result collision and publishability audit

**Search cutoff:** 2026-09-11  
**Audit state:** post-result, before manuscript drafting  
**Verdict:** `GO_WRITE_NARROW_FALSIFICATION_BOUNDARY`  
**Article type:** Original Research Article, not Review  

## 1. Frozen result audited

This audit evaluates the realized result, not the pre-result proposal. The accepted evidence is bound to Decision 010 (`cd6c20a085ccb233cdb0b1c83ce61d1d1bd9e1f9c0fe3e7884f6ef8653e439ff`), the immutable full-run manifest (`8da9b9c1d0425fff118d983a8491aceb57ae8b6085300eacf18cd8e0d7e6eac4`), the full summary (`784132d1df5ab6a99a934bf2c53d17c7d11f650beb4b733be222011f060986ae`), the external verification record (`b250b08c731b01488e4337c12915cd3b85d72f11df9c188166ad6662edeef5ea`), and the post-verification audit (`86909436314e036c4b424e18711b2dab857959a11972542484639078a488a182`).

The result shape is:

1. the prespecified propensity-instability and high-dimensional-propensity signatures reproduced in all four Biogen sentinel-by-mechanism cells for H1 and H2;
2. the sole low-dimensional candidate passed all eight MCAR admission cells but failed all eight MAR admission cells, so the 16-cell admission gate failed;
3. all eight regular direct-mask cells passed overlap eligibility, but direct translation was not adjudicated after admission failure and was only descriptively favorable in four of eight cells;
4. sparse hPPB/rPPB masks failed overlap and could not rescue the regular-anchor result; and
5. the public PXR screen-to-dose-response cascade triggered the prespecified `ABSTAIN_NONIDENTIFIABLE_FULL_LIBRARY_DOSE_RESPONSE` action.

The publishable object is therefore a falsification and decision-boundary study. It is not a positive correction study.

## 2. Result-targeted search and closest collisions

The live search combined the realized result terms with fixed-predictor risk, selectively observed outcomes, ADME endpoint availability, ExpansionRx, Biogen, known-truth MCAR/MAR sentinels, nuisance misspecification, assay cascades, dose-response promotion, overlap, negative controls, and abstention. Primary studies, preprints, and official source records were prioritized. The earlier pre-result audits were rechecked but were not treated as sufficient for this post-result gate.

| Closest work | Occupied territory | Remaining difference after seeing the result |
|---|---|---|
| Seal et al., [*Model Validation Protocols for Machine Learning in Small Molecule Drug Discovery*](https://doi.org/10.64898/2026.08.19.745868), bioRxiv, 2026 | Strongest direct framing and data collision: ExpansionRx and Biogen, unequal ADME assay coverage, realistic splitting, representation and evaluation failure modes | Does not estimate the finite-population loss of a fixed predictor under selectively available references, run the frozen propensity/loss nuisance factorial, require known-truth admission before direct-mask interpretation, or use PXR as an identification-based abstention case |
| Lim, [*Predicting Quantifiability from Primary Screens to Prioritize Dose-Response Profiling*](https://arxiv.org/abs/2608.26538), 2026 | Treats experimental measurability after a primary screen as a predictable quantity and a triage target | Predicts whether follow-up will be quantifiable; it does not correct or falsify fixed-predictor risk estimates, use an ExpansionRx-to-Biogen confirmation chain, or impose the present admission rule |
| Ulichney and Coston, [*Beyond the Training Distribution: Evaluating Predictions Under Distribution Shift and Selection Bias*](https://arxiv.org/abs/2606.14506), 2026 | Arbitrary fixed black-box risk under selective labels and covariate shift using influence-function/double-machine-learning estimation | Closes estimator novelty. It does not provide the ADMET-specific reference-availability map, independent sentinel admission experiment, or staged-assay abstention boundary |
| Xue et al., [*Probabilistic Contrastive Pretraining for Multi-task ADME Property Prediction*](https://arxiv.org/abs/2606.11508), 2026 | Joint use of ExpansionRx and Biogen, multi-task ADME prediction, and incomplete labels | Closes any claim of first combining these datasets or first handling their missing labels. Its target is predictor training, not fixed-model population-risk evaluation under selective reference availability |
| Boyer et al., [*Estimating and Evaluating Counterfactual Prediction Models*](https://doi.org/10.1002/sim.70287), 2025 | Fixed candidate prediction models, loss-based evaluation, inverse weighting, doubly robust estimation, and nuisance misspecification | General method territory is occupied; the ADMET measurement-process experiment, prospective admission failure, and PXR boundary remain different |
| Tian et al., [*Semi-supervised Regression Analysis with Model Misspecification and High-dimensional Data*](https://arxiv.org/abs/2406.13906), 2024 | AIPW, high-dimensional propensity/outcome nuisance functions, and model misspecification | Closes high-/low-dimensional nuisance-factorial novelty. The target parameter and ADMET-specific use are different |
| Liu and Kashima, [*Chemical Property Prediction Under Experimental Biases*](https://www.nature.com/articles/s41598-022-12116-5), 2022 | Molecular propensity scores, inverse weighting, and counterfactual regression under experimental sampling bias | Closes first-use-of-propensity-correction claims in molecular property prediction. It corrects training selection rather than estimating the risk of an already fixed predictor from selectively observed evaluation outcomes |
| Fang et al., [*Prospective Validation of Machine Learning Algorithms for ADME Prediction*](https://doi.org/10.1021/acs.jcim.3c00160), 2023 | Source study and industrial context for the six-endpoint Biogen collection and genuine prospective prediction assessment | The released 3,521-compound table does not encode that prospective sequence. The local 2,463/1,058 structure-grouped split is derived and cannot be called prospective, temporal, or official |
| OpenADMET, [PXR challenge description](https://openadmet.ghost.io/announcing-the-next-openadmet-blind-challenge-predicting-pxr-induction/) and [pinned dataset](https://doi.org/10.57967/hf/9731), 2026 | Documents the single-concentration-to-dose-response cascade, formal-hit rule, and additional analogue/diversity promotion routes | Supports the abstention case but cannot identify full-library dose-response risk because truth is absent for most non-advanced compounds and not every promotion variable is public |
| OpenADMET, [multi-fidelity, cost-aware active-learning record](https://doi.org/10.5281/zenodo.21789539), 2026 | Uses the same PXR single-point/dose-response hierarchy for multi-fidelity learning and experiment selection | Closes first-use-of-the-PXR-cascade claims; it does not estimate historical-selection-adjusted fixed-model risk or implement the present admission/abstention chain |

No exact end-to-end duplicate was located. This is a negative search result rather than proof of global novelty.

## 3. Claims closed by the collision audit

The manuscript must not claim:

- a new AIPW, IPW, double-machine-learning, or missing-outcome estimator;
- a new general ADMET model-validation framework;
- the first treatment of missing labels, assay coverage, or experimental selection in molecular property prediction;
- the first joint use of ExpansionRx and Biogen or the first computational use of the PXR assay cascade;
- MNAR correction, proof that PXR is MNAR, or recovery of full-library PXR pEC50/Emax risk;
- prospective, temporal, official, or transportability validation from the derived Biogen split;
- general superiority of low-dimensional propensity models or successful correction of real sparse endpoints; or
- clinical, wet-laboratory, or prospective drug-development utility.

## 4. Defensible contribution after narrowing

The remaining contribution is a linked, ADMET-specific empirical test of when an available statistical correction must not be interpreted:

> A development-derived selective-reference instability signature was tested under a frozen protocol in an independent public ADME collection. Although the mechanism signatures reproduced, the only prespecified candidate failed every MAR admission cell. Adequate overlap under regular empirical masks therefore did not authorize a correction claim, and a documented PXR assay cascade required abstention when full-library dose-response risk was not identified.

This contribution is narrower than a method or framework paper, but it is more than a data note or review because it contains a completed confirmatory computation, prespecified hypotheses and thresholds, known-truth controls, an independent confirmation dataset, retained negative results, and a real assay-process boundary analysis.

## 5. Reviewer-facing publishability test

| Likely objection | Required response in the manuscript |
|---|---|
| Seal et al. already use ExpansionRx and Biogen to study validation failure | Name this closest neighbor and distinguish the estimand and decision rule. Do not use broad “validation framework” language |
| AIPW and selective-label risk estimation are established | State explicitly that the estimators are probes, not contributions; cite the missing-outcome and double-robust literature |
| The candidate failed under the benign MAR generator | Make this the main result: an overlap- and coverage-looking analysis can still fail prospective known-truth admission |
| Transferred masks may not reproduce another endpoint's biology | Call them coverage- and shape-matched controls only; never use them to identify the real endpoint mechanism |
| The Biogen split is not prospective | Call it a frozen structure-grouped derived split and an independent public confirmation dataset, not prospective validation |
| PXR lacks truth for non-advanced compounds | Report selection strata and overlap diagnostics only, then the prespecified abstention action; do not impute a full-library effect |
| Positive H1/H2 could be mistaken for method success | Always report H1/H2 jointly with the failed admission, non-adjudicated direct translation, and PXR abstention |
| The result depends on high-dimensional numerical instability | Retain the complete factorial, exact row replay, frozen negative/full-loss results, and independent bit-level verification in the supplement and reproducibility record |

The result remains relevant to the JPA special issue because AI-assisted ADME development depends on experimental references that are generated by selective assay processes. The paper's practical lesson is a measurement and evaluation rule: diagnostic replication is insufficient; correction requires a separate, prospective admission test, and an unidentified target requires an explicit abstention.

## 6. Decision and writing contract

`STOP_DUPLICATED` is not supported because no located study combined the fixed-risk estimand, development-confirmation separation, known-truth admission failure, direct-mask non-translation, and PXR identification-based abstention.

`GO_WRITE_NARROW_FALSIFICATION_BOUNDARY` is supported subject to all of the following:

1. use an Original Research Article structure and keep literature synthesis subordinate to the empirical study;
2. title and abstract must foreground falsification/admission rather than correction performance;
3. report the complete joint outcome: H1/H2 pass, admission fail, direct overlap pass but translation not adjudicated, and PXR abstention;
4. present Seal et al. as the closest data/framing neighbor and Ulichney and Coston as the closest general fixed-risk method neighbor;
5. label 2026 preprints as preprints and refresh their status before submission;
6. include the complete negative factorial and stress-anchor evidence in the supplement; and
7. make code/data release and formal submission separate author-controlled gates.

This audit authorizes manuscript and reproducibility-package construction. It authorizes no new confirmatory model run, no result-dependent estimator revision, no raw-data redistribution beyond documented licenses, no public release, and no journal submission.
