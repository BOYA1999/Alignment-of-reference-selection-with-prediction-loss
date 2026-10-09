# Post-result C4 loss-nuisance forensic analysis

Parent: reviewer item C4, frozen `confirm_001/full_recovery_001`; manuscript analysis A3, claim C2. Role: supporting failure-mechanism evidence, main-text brief result and supplementary decomposition. All new results are post-result diagnostics and cannot replace the frozen confirmation or change its gates.

Question: which stored loss-nuisance predictions drive the extreme high-dimensional AIPW estimates, and do an independent numerical path and a bounded single-factor sensitivity reproduce or remove the extreme behavior?

Stage 1 inspects all 64 stored seed-41001 cells, decomposes mean predicted loss and mean residual correction, reports predicted-loss distributions and the largest absolute centered pseudo-outcomes, and reconstructs the high-dimensional features from the pinned source. Fold-level feature scaling and extrapolation are inspected without changing masks, folds, losses, or predictors.

Stage 2, after the stage-1 findings fix the diagnostic target, independently refits only the eight endpoint × anchor × mechanism cells for seed 41001, with two fixed loss columns. The original LSQR ridge is replayed, and a centered dual-system ridge solution using a dense positive-definite solver is compared. A single identified feature can be removed in a separately labeled sensitivity. No parameter search, model selection, or gate optimization is performed. Exact stage-2 definition is recorded before execution.

Resources: Windows CPU, 32 logical processors, 33.4 GB physical RAM, more than 500 GB free storage; dense problems contain at most 929 rows × 2266 features. BLAS threads capped at 2. No fixed-predictor retraining; one replay seed only. Full 100-seed numerical refits are outside this bounded diagnostic.

Inputs: frozen row replay, predictions, scenario table, fold assignments, feature contract, pinned Biogen derivative. Outputs remain under this directory; the 19-file frozen directory and its runner/verifiers remain read-only. Reconstructed features are local diagnostic intermediates containing structure-derived data and are not automatically public-package eligible.

Completion: decomposition, independent numerical comparison, sensitivity result, input/output hashes, unchanged frozen-file hashes, and a reviewer-facing interpretation with explicit limitations. Stop expanding when the failure mechanism is sufficiently localized or an unresolved discrepancy requires a distinct investigation.

## Stage-2 definition, recorded after inspection and before refits

The dominant extrapolation feature is `Ipc`: the dominant held-out compound has about 1.95 × 10^17, versus a maximum 8.81 × 10^9 among the corresponding nuisance-training rows; standardized held-out values exceed 5 × 10^8. Other inspected features have maximum absolute standardized held-out values below 84. The independent original refit uses the frozen alpha=10, LSQR tolerance=1e-6, shared two-output losses, and all original features. Numerical verification uses explicit observed-training centering and solves `(Xc Xc.T + 10 I) dual = yc` in double precision, followed by `coef = Xc.T dual`; its ridge normal-equation residual is recorded.

The ONE stability sensitivity clips only held-out `Ipc` to the minimum/maximum of that fold's complete nuisance-training feature rows before applying the unchanged fitted scaler and ridge. It does not change training inputs, labels, coefficients, any other feature, propensity, masks, folds, fixed predictors, or frozen results. This inspects the contribution of one identified extrapolation factor. It is a post-result sensitivity, not an optimized estimator, proof of general recovery, or new confirmation. Eight cells × three folds × two loss outputs are refit. The current stage uses two BLAS threads and only a few hundred MB; available system memory reported by root is approximately 16 GB.
