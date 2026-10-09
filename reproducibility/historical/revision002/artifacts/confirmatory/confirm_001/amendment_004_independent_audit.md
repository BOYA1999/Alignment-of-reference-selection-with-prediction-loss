# `confirm_001` Amendment 004 independent audit

Verdict: CLEAR_TO_IMPLEMENT_VERIFICATION_REPLAY_002
Failure record 002 SHA-256: 455c8e3928099fb31452801b5cec1c12d3ab716eaba1a03dd4ca8ee59a94e634
Protocol Amendment 004 SHA-256: 00482a584e320e7b8188fe29a6db9f1f5bf1349874062c473e1d407273e8f3e3
Config Amendment 004 SHA-256: 5be83a824a04f26729e1d51fe1aa547ebf750430ad1469a9dee7270b364a65df
Consumed Decision 008 SHA-256: 3cc6426bdebeb739d31d8d0cee1092d7de7ac8c4625b6b2004011390fec43216
Supplemental verifier 001 SHA-256: 4e95a20d5a3c74d5c96122b80fc37cf5d44892c4438e846b781899fbf2212a57
Original verifier SHA-256: 740e31590ec40f40f7783e59889d265c3d5d20296901c3968d4b593986238715
Original runner SHA-256: 650a8b8f6b9bf65db69efc7ba60340aafa99c3c4b8b5bba957a8e6c6fa800c5c
Full run manifest SHA-256: 8da9b9c1d0425fff118d983a8491aceb57ae8b6085300eacf18cd8e0d7e6eac4

Three independent read-only audits covered the authorization state machine and write boundaries, the numerical and memory-layout diagnosis, and preservation of the complete verification flow. They agreed that the sole new correction is a bit-preserving F-order conversion of the two independently source-reconstructed PXR feature matrices at their fixed context keys. The Amendment-003 Biogen low-dimensional F-order correction and its stored-generator durability checks must be inherited unchanged.

The new supplemental verifier must execute the unchanged original `verify_full` flow completely. It may not patch or replace `compare_frames`, `smd_map`, `verify_pxr_linkage`, any exact comparator, any tolerance, any threshold, or any scientific decision function. The original runner, original verifier, first supplemental verifier, consumed decisions, prior amendments and audits, and the immutable 19-file full output remain hash-bound and unchanged.

The failed first external verification path remains permanently absent. Any successful second verification record must use the separate Amendment-004-fixed write-once path and may be created exclusively only after every supplemental and original check passes. This audit clears implementation and subsequent static code review only. It does not authorize Decision 009, verifier execution, runner execution, any write inside the frozen full directory, creation of a PASS record, scientific interpretation, manuscript drafting, release, or submission.
