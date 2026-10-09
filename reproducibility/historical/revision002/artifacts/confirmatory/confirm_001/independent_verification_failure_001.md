# `confirm_001` independent full-verification failure 001

Recorded at: 2026-09-11T08:39:33+08:00

State: INDEPENDENT_VERIFICATION_FAILED_BEFORE_RESULT_ACCEPTANCE
Study ID: confirm_001
Runner output state: COMPLETED_PENDING_INDEPENDENT_FULL_VERIFICATION
Scientific interpretation state: STOP_IMPLEMENTATION_OR_DATA
Independent verification accepted: false

## Bound execution and evidence

Full output path: artifacts/confirmatory/confirm_001/full_recovery_001
Full output file count: 19
Full run manifest SHA-256: 8da9b9c1d0425fff118d983a8491aceb57ae8b6085300eacf18cd8e0d7e6eac4
Runner summary SHA-256: 784132d1df5ab6a99a934bf2c53d17c7d11f650beb4b733be222011f060986ae
Runner verification SHA-256: 7ef5b2cdb2d1885b074250402efd313aa9758b337b4f486147ddc10e0042af99
Generator vectors SHA-256: a6b817b9a3db207dd84f1a816a7374a5a21b186c853b73c2f2c422a1bfc830f7
Runner SHA-256: 650a8b8f6b9bf65db69efc7ba60340aafa99c3c4b8b5bba957a8e6c6fa800c5c
Verifier SHA-256: 740e31590ec40f40f7783e59889d265c3d5d20296901c3968d4b593986238715
Decision 006 SHA-256: 96cd052f18c9442ebe9b3b3524da7b88539c0896f6541dbc675cf9cb094cee0b
Decision 007 SHA-256: faea63b749016aaaf570888ca5cbc99fda856088c53ac02f489fca3df7f45fe7

The runner exited successfully; its manifest records `duration_seconds = 475.48431849479675`, while the console `DONE` line reported 475.648 seconds of outer elapsed time. The manifest binds 51 inputs and the 18 pre-manifest outputs. The full directory contains exactly the 19 expected files, all 18 stored output hashes match the current files, and no file has been added, removed, or changed after completion. The manifest status remains `COMPLETED_PENDING_INDEPENDENT_FULL_VERIFICATION`. Any runner-internal scientific decision remains quarantined and is not accepted evidence.

## Failed independent verification

Command: python src/verify_confirmatory.py --mode recovery-full --decision-006-sha 96cd052f18c9442ebe9b3b3524da7b88539c0896f6541dbc675cf9cb094cee0b --decision-007-sha faea63b749016aaaf570888ca5cbc99fda856088c53ac02f489fca3df7f45fe7
Exit code: 1
First failing check: generator probability hash mismatch
Failure location: src/verify_confirmatory.py:2073
Independent full-verification output file created: false

The first rejected cell was `LOG_HLM / LOG_SOL / mar`. All four MCAR generator probability vectors matched bit-for-bit. All four MAR vectors failed exact hash comparison:

| Endpoint / anchor | Stored hash | Reconstructed hash | Maximum absolute probability difference |
|---|---|---|---:|
| LOG_HLM / LOG_SOL | `13d646c9c3cee7bbc6bc9eb2f7fa5ac02f06b6ed713f4323cfb7a78bc5df95c1` | `955e9109c8f6c5add2bfcaf32e87bb25031db4935bd255f1dfe16dd3d757fead` | 4.3520742565306136e-14 |
| LOG_HLM / LOG_MDR1_ER | `fd0dc595d0441ea3efec295d229a6b78cc02a3ad0833bb957392c88ff7836c86` | `1917ed828fb52c7b8762b53e717c2412bbdc8a8030ebf15d49a36c9253a1d5b1` | 4.808616838047897e-10 |
| LOG_RLM / LOG_SOL | `c4d0a289a44125989614e20d96ed8ec6dd3d3afa44843ceca6a94ca38d3a8029` | `257aac98db5a68d88d588e02fd8f1c91bf2d71a6fada5a0e0399984c9f691dd5` | 4.440892098500626e-15 |
| LOG_RLM / LOG_MDR1_ER | `2a1ce55271a766e2d0fb99d4055a5fe3d0542eaa001e9111fa310ef9be5532f4` | `98637258c7be3870f258afe926f99b7ae28bbc66b875a604eac03730550ebdd9` | 4.8717696543576494e-12 |

## Read-only root-cause isolation

The runner selects target rows and then the ten descriptor columns, producing an F-contiguous `float64` low-dimensional matrix with strides `[8, 7432]` for HLM and `[8, 7376]` for RLM. The verifier selects the ten columns before target rows, producing a C-contiguous matrix with strides `[88, 8]`. The two matrices are elementwise bit-identical. Their different reduction order in `StandardScaler`, followed by amplification in liblinear, creates the MAR differences above.

Forcing the independently reconstructed matrix to F-contiguous order reproduces all four stored MAR score vectors, intercept shifts, probability vectors, and hashes bit-for-bit. Keeping it C-contiguous reproduces the rejected verifier hashes. Availability dtype, source values, row ordering, prediction-model execution order, OpenMP state, Parquet serialization, and hash encoding were excluded. Each stored probability hash is self-consistent with its Parquet float64 bytes, and stored `score + shift -> expit -> clip` regenerates the stored probabilities bit-for-bit.

Across the four MAR cells, 100 frozen seeds, and 370,200 Bernoulli threshold comparisons, the C- versus F-layout probability differences caused zero mask flips. This supports an implementation-only diagnosis but does not waive the frozen exact-hash and `1e-12` requirements.

## Required recovery boundary

`full_recovery_001` is immutable and may not be rerun, resumed, overwritten, or used as a destination. Decision 007 is consumed. Recovery may only proceed through an additive verifier-only amendment, independent amendment and code audits, and a new hash-bound decision. The original runner and verifier remain unchanged. A supplemental verifier may read the frozen sources and the immutable full directory, make the source-derived low-dimensional matrix explicitly F-contiguous without changing any element, preserve every model, seed, formula, hash check, and tolerance, and create at most one external verification record with exclusive creation. A second verification failure must remain a failure and requires another additive record; it may not trigger result repair or threshold relaxation.
