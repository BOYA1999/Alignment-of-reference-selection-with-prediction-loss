# `confirm_001` verification replay 002 post-verification audit

Verdict: ACCEPT_VERIFICATION_REPLAY_002_PASS
Verification record SHA-256: b250b08c731b01488e4337c12915cd3b85d72f11df9c188166ad6662edeef5ea
Verification status: PASS
Recorded acceptance state: PASS_PENDING_POST_VERIFICATION_AUDIT
Verified scientific decision: GO_WRITE_FALSIFICATION_BOUNDARY
Decision 009 SHA-256: 6120e6d5332a295a1396613b64758d29302595f662bcddbd9f1ea5b1205e92f3
Supplemental verifier 002 SHA-256: 4818856cd480688b4a85178c519023971ac4f48ad4ab6011a81b88827def1d9a
Full run manifest SHA-256: 8da9b9c1d0425fff118d983a8491aceb57ae8b6085300eacf18cd8e0d7e6eac4
Full output inventory: 19 files, 0 directories
Frozen input hashes: 51/51 PASS
Frozen output hashes: 18/18 PASS
Failed verification output status: ABSENT

Three independent read-only post-verification audits covered the external PASS record and its authorization chain, the frozen full-directory inventory and complete input/output hash maps, and the numerical plus code-identity evidence. They independently confirmed all eight stored generator probability-vector hashes, all four stored MAR score-plus-shift formula replays, the bit-preserving Biogen low-dimensional and PXR feature-layout conversions, and zero maximum absolute differences for the seed-41001 `p_raw`, `p_fitted`, and `mhat` nuisance refits. They also confirmed that the unchanged original full verification flow completed without skipped checks, the four frozen code files retained their authorized hashes, the immutable full directory remained unchanged, and the earlier failed verification-output path remained absent.

This audit accepts the bounded computational evidence and permits progression only to the post-result collision assessment and the subsequent writing gate. It does not authorize claims beyond the verified design, clinical or experimental validation, unrestricted generalization, public release, repository publication, journal submission, or any other external action.
