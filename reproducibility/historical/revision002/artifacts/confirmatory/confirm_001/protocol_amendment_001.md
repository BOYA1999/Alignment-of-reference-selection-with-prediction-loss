# `confirm_001` Pre-result Protocol Amendment 001

Timestamp: 2026-09-11T05:46:04+08:00  
State: `PENDING_INDEPENDENT_AMENDMENT_AUDIT`

## Trigger

A read-only readiness audit performed after the protocol/config freeze and before any `confirm_001` smoke or formal result detected an internally impossible hash contract. Six source-derived vectors in the earlier source manifest and config were hashed as UTF-8 line records with a terminal LF. Protocol Section 9 instead requires ordered text without a terminal LF and requires availability/formal-hit vectors to be serialized as raw `uint8` bytes. A runner cannot satisfy both representations at once.

The formal result directory did not exist at detection. No Biogen prediction, confirmatory mask realization, nuisance fit, factorial estimate, bootstrap result, PXR propensity fit, or gate result had been computed.

## Resolution

The immutable frozen protocol, config, source manifest, protocol audit, and decision record are retained byte-for-byte. The machine-readable overlay is `config_amendment_001.json`; the corrected source-vector expectations and complete old-to-new mapping are in `data/external/source_freeze_20260911/source_manifest_amendment_001.json`.

For `confirm_001`, the six corrected values supersede only the six named historical vector hashes. Both the old and corrected values remain recorded. Every source-file byte hash, source revision, row count, cohort definition, row order, model, feature, fold, seed, estimator, threshold, output grid, decision rule, and anti-claim remains unchanged.

The six corrections are mechanical:

| Vector | Historical SHA-256 | Corrected protocol-encoding SHA-256 |
|---|---|---|
| Biogen internal-ID plus SMILES text | `5b31164b...2d85` | `040e6d75...39f9` |
| Biogen cluster text | `c85847e6...a79` | `b3060140...ee0c` |
| Biogen split text | `8ea3ebb2...a79` | `7f0fe717...e827` |
| PXR parent OCNT-ID text | `df3e6a78...8402` | `7b9b757a...e168` |
| PXR follow-up availability bytes | `cd91cb26...3d62` | `74efc08d...2faa` |
| PXR formal-hit bytes | `38958812...4b02` | `584db5b1...515d` |

The historical text hashes were independently shown to equal the same ordered values with one terminal LF. The corrected hashes were independently recomputed from the frozen source files under the Section 9 encoding. The only other precomputed analysis vector, the bootstrap seed-index matrix, continues to match its frozen little-endian `int64` hash exactly.

For mechanical reproducibility under the already pinned RDKit `2026.03.3`, the 217-descriptor block is enumerated exactly as `sorted(Descriptors._descList, key=lambda item: item[0])`; the ordered descriptor-name vector, joined by LF without a terminal LF, has SHA-256 `84cea0bf91c3c74a704770fd5d5477dc689fc3c09e0997cd268d070acc61be73`. Each per-row `canonical_smiles_sha256` and `murcko_group_sha256` is the SHA-256 of that one UTF-8 string with no separator, newline, terminal LF, or metadata prefix.

## Execution gate

No smoke or formal execution is allowed until an independent audit verifies the overlay hashes, exhaustiveness of the six-item correction set, unchanged base artifacts, exact analysis environment, absent formal directory, and runner/verifier use of the amended expectations. Clearance must be recorded additively; this amendment does not itself authorize execution.
