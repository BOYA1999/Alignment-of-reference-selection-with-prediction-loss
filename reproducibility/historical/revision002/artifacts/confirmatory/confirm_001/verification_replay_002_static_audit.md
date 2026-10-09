# `confirm_001` verification replay 002 static audit

Verdict: CLEAR_TO_VERIFICATION_REPLAY_002
Supplemental verifier 002 SHA-256: 4818856cd480688b4a85178c519023971ac4f48ad4ab6011a81b88827def1d9a
Supplemental verifier 001 SHA-256: 4e95a20d5a3c74d5c96122b80fc37cf5d44892c4438e846b781899fbf2212a57
Original verifier SHA-256: 740e31590ec40f40f7783e59889d265c3d5d20296901c3968d4b593986238715
Original runner SHA-256: 650a8b8f6b9bf65db69efc7ba60340aafa99c3c4b8b5bba957a8e6c6fa800c5c
Failure record 002 SHA-256: 455c8e3928099fb31452801b5cec1c12d3ab716eaba1a03dd4ca8ee59a94e634
Protocol Amendment 004 SHA-256: 00482a584e320e7b8188fe29a6db9f1f5bf1349874062c473e1d407273e8f3e3
Config Amendment 004 SHA-256: 5be83a824a04f26729e1d51fe1aa547ebf750430ad1469a9dee7270b364a65df
Amendment 004 audit SHA-256: 104484171d0e567463bf8584c49fcf84ea2bd8d36e0f9bf29a4b413f0e040373
Full run manifest SHA-256: 8da9b9c1d0425fff118d983a8491aceb57ae8b6085300eacf18cd8e0d7e6eac4

Three independent read-only audits covered the code and authorization state machine, the numerical and array-layout contract, and the original-verifier flow plus full-directory immutability. They found exactly two loader monkeypatches and no replacement of `compare_frames`, `smd_map`, `verify_pxr_linkage`, or `verify_full`. The runner is hash-bound but is never imported or invoked, and the completed full directory has no write path.

The inherited Amendment-003 Biogen low-dimensional conversion remains restricted to `lowdim` and `representations.lowdim_11`. The sole added Amendment-004 conversion remains restricted to `features.pxr_structure_10` and `features.pxr_structure_screen_12`. Both corrections require source-derived `float64` values, frozen C- and F-order strides, unchanged shape, and elementwise bit identity; the PXR arrays additionally require bit identity with their independently source-derived `expected_linkage` columns. Stored-generator probability rehashing and MAR score-plus-shift bit replay remain bit exact. The unchanged original `verify_full` is called in full with the original exact SMD JSON, hash, grid, formula, threshold, boundary, and decision checks.

The failed verification-output path and the new verification-output path were both absent during static clearance and must remain absent before execution. Only a separately hashed Decision 009 may authorize one read-only replay, and only complete success may exclusively create the new external PASS record. This static clearance authorizes no runner execution, no mutation of `full_recovery_001`, and no scientific interpretation; the quarantined result remains uninterpretable until replay 002 passes and its output and immutable state pass an independent post-verification audit.
