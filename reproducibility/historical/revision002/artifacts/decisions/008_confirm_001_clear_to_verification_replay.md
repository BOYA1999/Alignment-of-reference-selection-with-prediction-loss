# `confirm_001` Decision 008

Decision: CLEAR_TO_VERIFICATION_REPLAY
Mode: verification-replay
Supplemental verifier SHA-256: 4e95a20d5a3c74d5c96122b80fc37cf5d44892c4438e846b781899fbf2212a57
Original runner SHA-256: 650a8b8f6b9bf65db69efc7ba60340aafa99c3c4b8b5bba957a8e6c6fa800c5c
Original verifier SHA-256: 740e31590ec40f40f7783e59889d265c3d5d20296901c3968d4b593986238715
Decision 006 SHA-256: 96cd052f18c9442ebe9b3b3524da7b88539c0896f6541dbc675cf9cb094cee0b
Consumed Decision 007 SHA-256: faea63b749016aaaf570888ca5cbc99fda856088c53ac02f489fca3df7f45fe7
Failure record SHA-256: 3a6ac0a67988d3037831f9f21a084b814f78e701ec464c5e224b09018e74c8a1
Protocol Amendment 003 SHA-256: 458d0217da500b89a89c28018e07f2df18b3d741f1fb13e068451cef0a447adb
Config Amendment 003 SHA-256: 690875bca0afcd95a44946c57aa552c3287c1c80efb690cb412318565999ab60
Amendment 003 audit SHA-256: 61a67e578d57ce790182646c61ba8da18ce8f49d7fb9f6853f4033b21156fe75
Verification-replay static audit SHA-256: df4c666f97390258e164d7c65cda527d388d7e82f7fefbfc9ba44cd6010ffe35
Verification-replay static audit verdict: CLEAR_TO_VERIFICATION_REPLAY
Full output path: artifacts/confirmatory/confirm_001/full_recovery_001
Full output file count: 19
Full run manifest SHA-256: 8da9b9c1d0425fff118d983a8491aceb57ae8b6085300eacf18cd8e0d7e6eac4
Verification output path: artifacts/confirmatory/confirm_001/full_recovery_001_independent_verification.json

## Basis

The Amendment 003 independent audit cleared implementation of the narrowly scoped supplemental verifier, and three independent read-only static audits subsequently cleared its authorization state machine, numerical and memory-layout correction, complete reuse of the unchanged original `verify_full` flow, failure behavior, and immutable-output boundary. The supplemental verifier, all upstream authorizations and amendments, the original code pair, and the frozen full-run manifest are hash-bound above.

## Authorized action

This decision authorizes exactly one verifier-only replay against the existing immutable `full_recovery_001` directory. The replay may read and independently reconstruct the frozen evidence and, only after every original and supplemental check passes, exclusively create the external verification record named above.

This decision authorizes no runner execution and no write, append, overwrite, resume, rename, deletion, or other mutation inside `full_recovery_001`. It does not accept or interpret the quarantined scientific result. Any scientific interpretation remains prohibited until the replay returns PASS and an independent post-verification audit confirms both the external record and continued frozen-directory immutability.
