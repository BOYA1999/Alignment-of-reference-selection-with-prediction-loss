# Independent pre-result audit of `confirm_001`

Date: 2026-09-11  
Verdict: `CLEAR_TO_FREEZE`

## Audited stable snapshot

- Decision 003 SHA-256: `9e4afd7786878eb9764cd9f91385ad97c990cfcc39f4b3155f164edfec1dd3f7`
- Protocol SHA-256: `c1bc35298b0c15cec6a479032d37622f424e1d80a9f2b73ea5ef61939310d351`
- Config SHA-256: `95492f7b57a5c5d77f5c54de7906ab9a1508f46b6d0b1a959948efc9953de9b8`
- Exact-claim collision addendum SHA-256: `d141c86b381c4febdc050cf47713afe847bad84a01f352f50bfd30ea24149350`

The config parsed without duplicate keys and its recorded protocol hash matched the audited protocol. The formal result directory did not exist and no Biogen prediction, mask realization, nuisance estimate, risk estimate, bootstrap result, or PXR propensity fit existed at audit time.

## Scope and conclusion

The audit independently checked source and license identities, development-confirmation separation, the 2 x 2 propensity-feature by loss-feature factorial, H1/H2 aggregation without duplicate oracle rows, random-number scope, MAR intercept solver, estimator and interval formulas, unscaled-feature SMD, clipping and weight definitions, direct-mask output versus interpretation gates, decision priority, PXR parent/join/feature/stratum contracts, row identities, vector serialization, expected output grids, and post-result claim restrictions.

All identified ambiguities and arithmetic errors were corrected before this verdict. The final remaining vector-hash ambiguity was resolved by specifying little-endian C-order raw bytes for numeric vectors and UTF-8 LF joining without a terminal LF for ordered text vectors.

`CLEAR_TO_FREEZE` authorized only the mechanical transition from the audited `READY_FOR_PROTOCOL_AUDIT` snapshot to `FROZEN_PRE_RESULT` with `formal_run_authorized=true`. It did not authorize a run before synthetic smoke, independent source-code audit, replay-verifier readiness, and an absent or empty formal directory. It did not broaden the narrow novelty claim or alter any scientific threshold.
