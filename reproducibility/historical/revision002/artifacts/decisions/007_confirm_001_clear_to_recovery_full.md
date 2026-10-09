# Decision 007 — `confirm_001` recovery-full clearance

Decision date: 2026-09-11

Decision: CLEAR_TO_RECOVERY_FULL
Mode: recovery-full
Runner SHA-256: 650a8b8f6b9bf65db69efc7ba60340aafa99c3c4b8b5bba957a8e6c6fa800c5c
Verifier SHA-256: 740e31590ec40f40f7783e59889d265c3d5d20296901c3968d4b593986238715
Decision 006 SHA-256: 96cd052f18c9442ebe9b3b3524da7b88539c0896f6541dbc675cf9cb094cee0b
Recovery smoke status: PASS
Independent recovery smoke verifier status: PASS
Recovery smoke output path: artifacts/confirmatory/confirm_001/recovery_smoke_001
Recovery full output path: artifacts/confirmatory/confirm_001/full_recovery_001
Failed output path: artifacts/confirmatory/confirm_001/full
Recovery smoke JSON SHA-256: e94107b10faa947a883f909a0bf60605799be72ad18dcf41076883fd157b1839
Recovery smoke run manifest SHA-256: 160d8a21d6a97bf2ed424a25bd69c901673ea1746a0ff51f8437e2bda4e296ff
Recovery roundtrip fixture SHA-256: 3961cd5b2acc70c90a10f94fd5ea17e6a382291c873401ef4beeb1765b5b6dfb
Independent recovery smoke verification SHA-256: 1e24af7b2659a8e4cf67a58ecab656c7e306567a0a80acf0351d9116ecb36055
Independent recovery smoke verification path: artifacts/confirmatory/confirm_001/recovery_smoke_001_independent_verification.json

## Basis

The Decision-006-authorized recovery smoke completed at the sole fixed path and remained synthetic and non-adjudicative. Its ten checks passed, the large-ratio fixture restored all component floats bit-exactly with zero derived-ratio difference, and the independent verifier reconstructed the contract. Three subsequent read-only audits independently closed the exact three-file inventory, 46 frozen input hashes, current code and environment identities, Decision 006, smoke payload and manifest, and the exclusive fixed-path verification record. Both programs rejected recovery-full while this decision was absent.

## Scope

This decision authorizes exactly one full recomputation from the frozen sources at `artifacts/confirmatory/confirm_001/full_recovery_001`, using the frozen scientific contract and the code identities above. It does not predetermine a favorable result. Negative, null, abstention, falsification, or stop outcomes must be retained exactly as computed. The old failed `full/` directory remains immutable evidence and is not a scientific input.

The recovery-full path is write-once. If execution fails after its creation, it may not be deleted, overwritten, resumed, or reused; a new additive failure record and recovery decision would be required. Independent full verification is mandatory before interpretation, writing, release, or submission.
