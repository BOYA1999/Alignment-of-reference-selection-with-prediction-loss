# `confirm_001` exact-claim collision addendum

**Search cutoff:** 2026-09-11  
**Search type:** live, claim-targeted primary-paper/preprint and official-source check  
**Collision-only verdict:** `CLEAR_NARROW_CLAIM`  
**Run authorization:** none; protocol/config/source-freeze review remains a separate gate

## Exact claim audited

The candidate contribution is not a new estimator. It is the following linked, ADMET-specific validation architecture:

1. evaluate the population MAE of already fixed molecular predictors when reference values are selectively available;
2. use complete ADMET sentinel endpoints to create known-truth MCAR and covariate-conditional MAR controls;
3. cross a low- versus high-dimensional propensity representation with a low- versus high-dimensional conditional-loss representation in a complete `2 x 2` factorial, while keeping masks, folds, losses and prediction models shared;
4. treat the factorial as a specification diagnostic, with only the prospectively named propensity contrasts eligible for confirmation;
5. apply a prespecified sentinel admission rule before interpreting correction under direct empirical reference-availability masks, with passage explicitly necessary but not sufficient for MAR, positivity or cross-endpoint transport; and
6. use the recorded OpenADMET PXR screen-to-dose-response process as an abstention/non-identifiability boundary, never as proof of MNAR or recovery of unobserved full-library dose-response risk.

The live search combined terms for fixed-model risk, selective labels, missing outcomes, AIPW/DML, propensity and outcome nuisance specifications, high-dimensional variable selection, synthetic or negative controls, ADMET reference availability, assay cascades, admission rules, overlap and abstention. The searches also rechecked the closest works identified in the broader collision audit. No confirmation result was inspected or generated for this addendum.

The frozen primary statistic `median(abs(estimate - finite-population truth))` is a median absolute **estimation error**, not statistical bias in the usual repeated-sampling sense. Paper-facing claims should therefore use “absolute error,” “truth-recovery error,” or “instability,” unless signed repeated-sampling bias is separately defined and frozen.

## Collision map

| Primary source | What is already occupied | What remains distinct in the exact claim |
|---|---|---|
| Kang and Schafer, [Demystifying Double Robustness](https://arxiv.org/abs/0804.2958) | Population-mean estimation with incomplete outcomes; separate propensity and outcome models; finite-sample failure when both working models are wrong; sensitivity to small fitted propensities | A `2 x 2` nuisance comparison and finite-sample AIPW instability are not novel. The remaining contribution must be the ADMET measurement-process use of prospectively gated, known-truth sentinel controls, not the factorial itself |
| Ulichney and Coston, [Beyond the Training Distribution](https://arxiv.org/abs/2606.14506) | Closest fixed-black-box setting: general-loss target risk, selective labels, covariate shift, influence-function/DML correction, cross-fitting and stress tests | Closes method novelty. It does not supply an ADMET reference-availability map, an independent ADMET sentinel admission rule, or a staged-assay abstention case |
| Gronsbell et al., [Efficient Evaluation of Prediction Rules in Semi-Supervised Settings under Stratified Sampling](https://pmc.ncbi.nlm.nih.gov/articles/PMC9586151/) | Prediction-rule performance with nonuniformly sampled labels, weighting, augmentation, cross-validation and full-population evaluation | Closes generic missing-label performance-estimation and augmentation novelty. Sampling probabilities arise from a known design, and there is no molecular feature factorial, empirical availability-mask admission gate or assay-cascade boundary |
| Corbin, Baiocchi and Chen, [Avoiding Biased Clinical Machine Learning Model Performance Estimates in the Presence of Label Selection](https://arxiv.org/abs/2209.09188) | Fixed-model evaluation under selection completely at random, selection at random and selection not at random; known-truth simulations; recovery using properly specified weights | Closes generic selective-label and known-probability recovery claims. The outcome type, molecular availability masks, nuisance-factorial diagnosis and admission-before-translation workflow differ |
| Li et al., [Evaluation of predictive model performance of an existing model in the presence of missing data](https://pmc.ncbi.nlm.nih.gov/articles/PMC8985431/) | Evaluation of an existing prediction model with IPW, AIPW and imputation under MCAR/MAR/MNAR and nuisance-model specification experiments | A close evaluation-design precedent, but its missingness concerns prediction variables rather than selectively unavailable reference outcomes, and it has no molecular representation/admission/PXR chain |
| Kossen et al., [Active Testing](https://proceedings.mlr.press/v139/kossen21a.html), and Ashouritaklimi et al., [Prediction-Powered Active Testing](https://arxiv.org/abs/2607.08347) | Risk estimation when test labels are deliberately acquired with known, positive sampling probabilities; unbiased weighting and variance reduction | Historical assay availability is not investigator-randomized active testing. The candidate must not claim a new label-acquisition or prediction-powered estimator |
| Zou, Witten and Williamson, [Generalized Prediction-Powered Inference](https://arxiv.org/abs/2602.10332) | Partially observed outcomes, semiparametric efficiency and prediction-powered evaluation of classifier metrics | Further closes generic partially observed-outcome inference novelty; no ADMET sentinel gate or assay-cascade identification audit |
| Liu and Kashima, [Chemical Property Prediction Under Experimental Biases](https://www.nature.com/articles/s41598-022-12116-5) | Molecular propensity estimation and inverse-propensity correction under simulated experimental sampling bias | Corrects biased model-training data rather than estimating the risk of a fixed predictor from selectively observed evaluation outcomes; no sentinel admission rule |
| Tian, Wu and Tan, [Semi-supervised Regression Analysis with Model Misspecification and High-dimensional Data](https://arxiv.org/abs/2406.13906) | High-dimensional propensity and outcome nuisance estimation for AIPW under semi-supervised/covariate-shift settings | Closes any claim that combining high-dimensional nuisance models with AIPW is new; target parameters and ADMET validation architecture differ |
| Cho and Yang, [Variable Selection for Doubly Robust Causal Inference](https://arxiv.org/abs/2301.11094), and Tang et al., [Ultra-high Dimensional Variable Selection for Doubly Robust Causal Inference](https://arxiv.org/abs/2007.14190) | Separate roles of variables and feature selection in propensity and outcome models; high-dimensional DR inference | Reinforces that low/full feature contrasts are specification choices, not a novel causal mechanism or generally optimal representation |
| Wyss et al., [Synthetic Negative Controls](https://pubmed.ncbi.nlm.nih.gov/35439779/) and [Tuning LASSO Models for Propensity Score Weighting and Using Synthetic Negative Control Exposures for Residual Bias Detection](https://pmc.ncbi.nlm.nih.gov/articles/PMC13055431/) | Synthetic controls used to screen propensity analyses; balance diagnostics; high-dimensional propensity regularization; bias and coverage operating characteristics | Closest precedent for the admission logic. It concerns causal treatment-effect analyses in healthcare databases, not fixed-predictor risk under molecular reference availability. Therefore the gate itself cannot be presented as a general methodological invention |
| Coston, [Falsifying Discriminant Validity of Predictive Algorithms](https://arxiv.org/abs/2601.17146) | Falsification before downstream interpretation using proxy outcomes and calibrated predictive losses | Closes broad claims of inventing predictive-model falsification. The candidate controls the validity of a missing-reference risk estimator, not discriminant validity of a predictor |
| Dervovic and Cashmore, [Model Evaluation in the Dark](https://proceedings.mlr.press/v258/dervovic25a.html) | Classifier performance evaluation with missing labels using multiple imputation and uncertainty distributions | Establishes another direct missing-label evaluation neighborhood, but not propensity/AIPW fixed-risk estimation or a drug-analysis admission architecture |
| Martin and Bowen, [The Missing Label Problem](https://doi.org/10.26434/chemrxiv.9758423.v1) | Drug-screening missing-label bias, especially treating untested binding pairs as negatives, with mitigation for ligand-based virtual screening | Closes broad “missing labels in drug discovery” novelty. It addresses training/ranking assumptions rather than selectively observed reference outcomes for a fixed predictor's population risk |
| Seal et al., [Model Validation Protocols for Machine Learning in Small Molecule Drug Discovery](https://www.biorxiv.org/content/10.64898/2026.08.19.745868v1) | Broad small-molecule model-validation workflow, ADME assay coverage and split/representation/evaluation diagnostics | Closes a general ADMET validation-framework claim. It does not estimate fixed-model risk under selectively missing references or impose the proposed sentinel admission boundary |
| Fang et al., [Prospective Validation of Machine Learning Algorithms for ADME Prediction](https://pubs.acs.org/doi/10.1021/acs.jcim.3c00160) | The Biogen public six-endpoint industrial ADME collection and prospective predictive-model assessment | Supplies the external chemical data context, not the proposed selective-reference risk estimand, known-truth controls or admission rule. A locally derived split must not be called the paper's prospective split |
| OpenADMET, [PXR challenge description](https://openadmet.ghost.io/announcing-the-next-openadmet-blind-challenge-predicting-pxr-induction/) and [pinned public dataset](https://huggingface.co/datasets/openadmet/pxr-challenge-train-test) | A documented primary-screen to dose-response cascade: formal low-dose hits are followed up, with additional analogue/diversity routes; public screen and dose-response files | Supports an empirical abstention example. It cannot supply full-library dose-response truth or validate a correction estimator because most non-advanced compounds have no dose-response reference and not all promotion variables are public |

## Direct-duplicate finding

No primary paper or preprint located in this live check implemented the full combination of:

- fixed-predictor ADMET risk under selectively available references;
- independent real ADME sentinel endpoints with known full-cohort truth;
- a shared-mask `2 x 2` propensity-feature by conditional-loss-feature factorial;
- a prospective, necessary-not-sufficient admission gate before direct-mask interpretation; and
- a documented drug-discovery assay cascade whose prespecified action is abstention when support or identification fails.

This is a negative search finding, not proof of global novelty. Every statistical component has precedent. The defensible unit of novelty is the **linked ADMET measurement-process validation and abstention architecture**, evaluated prospectively on an independent public collection.

## PXR boundary required by the exact claim

The official challenge description states that formal hits satisfy `log2 fold change > 1` and `FDR < 0.05`, that every formal hit at the lowest nominal screening dose was advanced to dose response, and that analogue and diversity branches were also advanced. The pinned public-file audit therefore supports a public-file-specific availability indicator joined by `OCNT_ID`, not reconstruction of every official promotion decision.

The PXR component remains non-identifying even if a structure or screen-feature propensity model appears well calibrated. Dose-response truth is absent for most non-advanced compounds, and unrecorded analogue, diversity, QC or discretionary routes may affect advancement. The admissible paper-facing conclusion is a selection/overlap map followed by `ABSTAIN_NONIDENTIFIABLE_FULL_LIBRARY_DOSE_RESPONSE` or an explicitly overlap-restricted description. It is not a corrected full-library pEC50/Emax estimate.

Before protocol freeze, any model-based PXR overlap output must uniquely name its cohort, predictors, preprocessing, folds, clipping rule, ESS denominator, balance features and prespecified strata. The most direct empirical strata are the frozen formal-hit and non-hit groups; additional response/FDR bins must be fixed before their results are calculated. If those items remain unspecified, only the raw advancement-by-formal-hit table can be treated as protocol-determined.

## Claim language allowed after a successful confirmation

Subject to the frozen gates actually passing, the narrow claims may be:

1. the two development-derived propensity-instability signatures replicated on the prespecified Biogen confirmation design;
2. under the study's known-truth MCAR and benign low-dimensional MAR generators, the prespecified feature representation affected finite-sample fitted-propensity risk correction while the loss representation was held fixed;
3. the prespecified sentinel gate admitted or rejected the candidate specification in this dataset; and
4. the PXR public files illustrated why the workflow abstains from full-library dose-response correction under incomplete support and identification.

Even after a pass, the study may not claim that low-dimensional propensity models are generally superior, that the transferred empirical mask is MAR, that gate passage proves transportability, that PXR is physically MNAR, or that the architecture is universally validated.

## Result-after narrative firewall

- ExpansionRx `pilot_001` and `forensic_001` are development evidence only; their observations, seeds and thresholds cannot be counted as external confirmation.
- The `2 x 2` factorial must not be called a newly discovered method. Only the propensity contrasts named before Biogen results may adjudicate confirmation; loss-representation and interaction results remain exploratory regardless of direction.
- The two propensity signatures share the same high-dimensional fitted-propensity denominator and are not independent confirmations. Their conjunction is a stringent decision rule, not two statistically independent discoveries.
- A favorable oracle/fitted or low-/high-dimensional error ratio does not by itself establish practically important harm: the numerator, denominator, and their absolute difference must be reported, and a ratio driven by uniformly negligible errors may not be narrated as consequential instability.
- A direct-mask correction may not be interpreted merely because it is numerically favorable. It must first pass the known-truth admission gate, and passage remains necessary rather than sufficient.
- A failed endpoint, anchor, model or mechanism may not be rescued by pooled reporting, stress-anchor substitution, clipping/penalty changes or post-result model selection.
- PXR abstention is a design boundary already implied by documented selection and missing full-library truth; it cannot be advertised after the fact as an unexpected positive empirical discovery.
- Any manuscript-level novelty statement must be refreshed against papers appearing after this 2026-09-11 cutoff.

## Collision decision

`STOP_DUPLICATED` is not supported: no directly isomorphic ADMET study was found. `CLEAR_NARROW_CLAIM` is supported only for the linked domain-specific architecture above. This collision result does not by itself make the current draft protocol freeze-ready; exact computation and source-contract ambiguities must be resolved independently before `CLEAR_TO_FREEZE`.
