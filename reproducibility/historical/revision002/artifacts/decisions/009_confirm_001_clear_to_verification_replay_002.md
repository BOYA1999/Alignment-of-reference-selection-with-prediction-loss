# `confirm_001` Decision 009

Decision: CLEAR_TO_VERIFICATION_REPLAY_002
Mode: verification-replay-002
Supplemental verifier 002 SHA-256: 4818856cd480688b4a85178c519023971ac4f48ad4ab6011a81b88827def1d9a
Supplemental verifier 001 SHA-256: 4e95a20d5a3c74d5c96122b80fc37cf5d44892c4438e846b781899fbf2212a57
Original runner SHA-256: 650a8b8f6b9bf65db69efc7ba60340aafa99c3c4b8b5bba957a8e6c6fa800c5c
Original verifier SHA-256: 740e31590ec40f40f7783e59889d265c3d5d20296901c3968d4b593986238715
Decision 006 SHA-256: 96cd052f18c9442ebe9b3b3524da7b88539c0896f6541dbc675cf9cb094cee0b
Consumed Decision 007 SHA-256: faea63b749016aaaf570888ca5cbc99fda856088c53ac02f489fca3df7f45fe7
Consumed Decision 008 SHA-256: 3cc6426bdebeb739d31d8d0cee1092d7de7ac8c4625b6b2004011390fec43216
Failure record 001 SHA-256: 3a6ac0a67988d3037831f9f21a084b814f78e701ec464c5e224b09018e74c8a1
Failure record 002 SHA-256: 455c8e3928099fb31452801b5cec1c12d3ab716eaba1a03dd4ca8ee59a94e634
Protocol Amendment 003 SHA-256: 458d0217da500b89a89c28018e07f2df18b3d741f1fb13e068451cef0a447adb
Config Amendment 003 SHA-256: 690875bca0afcd95a44946c57aa552c3287c1c80efb690cb412318565999ab60
Amendment 003 audit SHA-256: 61a67e578d57ce790182646c61ba8da18ce8f49d7fb9f6853f4033b21156fe75
Protocol Amendment 004 SHA-256: 00482a584e320e7b8188fe29a6db9f1f5bf1349874062c473e1d407273e8f3e3
Config Amendment 004 SHA-256: 5be83a824a04f26729e1d51fe1aa547ebf750430ad1469a9dee7270b364a65df
Amendment 004 audit SHA-256: 104484171d0e567463bf8584c49fcf84ea2bd8d36e0f9bf29a4b413f0e040373
Verification-replay 001 static audit SHA-256: df4c666f97390258e164d7c65cda527d388d7e82f7fefbfc9ba44cd6010ffe35
Verification-replay 002 static audit SHA-256: dc52aaa99a5d3d06f377046878a10feb756da8bb807754c410861bd889ec26ae
Verification-replay 002 static audit verdict: CLEAR_TO_VERIFICATION_REPLAY_002
Full output path: artifacts/confirmatory/confirm_001/full_recovery_001
Full output file count: 19
Full run manifest SHA-256: 8da9b9c1d0425fff118d983a8491aceb57ae8b6085300eacf18cd8e0d7e6eac4
Failed verification output path: artifacts/confirmatory/confirm_001/full_recovery_001_independent_verification.json
Failed verification output status: ABSENT
Verification output path: artifacts/confirmatory/confirm_001/full_recovery_001_independent_verification_002.json

## Basis

Failure record 002 freezes the Decision-008-authorized replay failure before result acceptance and proves continued full-output immutability. Amendment 004 and its independent audit isolate the failure to a layout-triggered reduction-order defect in two bit-identical PXR source matrices. Three independent static audits cleared the new verifier's strictly bounded, bit-preserving PXR correction, inherited low-dimensional correction, complete reuse of the unchanged original verifier, authorization chain, failure behavior, and write boundary.

## Authorized action

This decision authorizes exactly one execution of `src/verify_confirmatory_replay_002.py` against the existing immutable `full_recovery_001` directory. It authorizes no runner execution. It authorizes no write, append, overwrite, resume, rename, deletion, or other mutation inside the full directory. The failed first verification path must remain absent.

Only after every supplemental check and every unchanged original `verify_full` check succeeds may the verifier exclusively create the separate external verification record named above. A PASS record remains pending independent post-verification audit and does not by itself authorize scientific interpretation, manuscript drafting, public release, or submission. If the execution fails, Decision 009 is consumed and may not be reused.
