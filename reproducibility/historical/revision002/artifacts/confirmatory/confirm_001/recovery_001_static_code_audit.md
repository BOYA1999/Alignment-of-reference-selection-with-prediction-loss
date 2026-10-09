# `confirm_001` Recovery 001 static code audit

Audit date: 2026-09-11

Verdict: CLEAR_TO_RECOVERY_SMOKE
Runner SHA-256: 650a8b8f6b9bf65db69efc7ba60340aafa99c3c4b8b5bba957a8e6c6fa800c5c
Verifier SHA-256: 740e31590ec40f40f7783e59889d265c3d5d20296901c3968d4b593986238715

Three independent review tracks passed: implementation/static fidelity, authorization/verifier closure, and numerical/scientific-contract fidelity. The repaired code closes the recovery-smoke and D007 evidence schemas, fixed paths, code/input/environment identities, synthetic checks, independent-verification scope, and exclusive write-once creation. Frozen scientific computations and source readers remain unchanged; generated result CSVs use `%.17g` writes and `round_trip` reads without relaxing formula tolerance. The failed full directory remains evidence-only and byte-identical.

Non-blocking P2 boundary: the failed-output runtime inventory checks exact files, byte counts, and hashes but would not reject an added empty directory or empty reparse directory. The audited tree currently contains neither, and failed outputs are not scientific inputs.

This verdict is a static-code clearance only. Recovery smoke execution still requires a separately hashed Decision 006 and remains synthetic and non-adjudicative; it does not authorize the recovery full run, interpretation, writing, release, or submission.
