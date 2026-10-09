# `confirm_001` verification replay static audit

Verdict: CLEAR_TO_VERIFICATION_REPLAY
Supplemental verifier SHA-256: 4e95a20d5a3c74d5c96122b80fc37cf5d44892c4438e846b781899fbf2212a57
Original verifier SHA-256: 740e31590ec40f40f7783e59889d265c3d5d20296901c3968d4b593986238715
Original runner SHA-256: 650a8b8f6b9bf65db69efc7ba60340aafa99c3c4b8b5bba957a8e6c6fa800c5c
Failure record SHA-256: 3a6ac0a67988d3037831f9f21a084b814f78e701ec464c5e224b09018e74c8a1
Protocol Amendment 003 SHA-256: 458d0217da500b89a89c28018e07f2df18b3d741f1fb13e068451cef0a447adb
Config Amendment 003 SHA-256: 690875bca0afcd95a44946c57aa552c3287c1c80efb690cb412318565999ab60
Amendment 003 audit SHA-256: 61a67e578d57ce790182646c61ba8da18ce8f49d7fb9f6853f4033b21156fe75
Full run manifest SHA-256: 8da9b9c1d0425fff118d983a8491aceb57ae8b6085300eacf18cd8e0d7e6eac4

Three independent read-only audits covered complementary surfaces. The code, security, and state-machine audit checked hash-bound authorization, fail-closed exception handling, pre-existing-output rejection, and exclusive creation of the external verification record only after every check succeeds. The numerical and layout audit checked that the sole correction is the F-contiguous reconstruction of the two Biogen low-dimensional matrices, with frozen pre- and post-conversion strides, unchanged shape and `float64` dtype, and elementwise bit identity; it also checked bit-exact rehashing of all eight stored generator probability vectors and bit-exact stored-score-plus-shift reconstruction for all four MAR cells. The verifier-flow and immutability audit checked that the unchanged original `verify_full` path is called in full without skipped or replaced downstream checks, that the runner is hash-bound but never imported or invoked, and that the immutable full directory is read only behind the exact 19-file inventory and pinned-manifest gates.

This clearance authorizes one verifier-only replay after a separately hashed Decision 008. It authorizes no runner execution and no write, append, overwrite, resume, or other mutation inside `full_recovery_001`. The quarantined computational result remains scientifically uninterpretable until the supplemental replay completes with PASS and the resulting external verification record and frozen-directory immutability both pass an independent post-verification audit.
