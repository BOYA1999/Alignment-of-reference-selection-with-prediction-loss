# Decision 001 — topic audit and pivot

Date: 2026-09-10

## Decision

`GO_PILOT`

The original AI-DMV framework is `STOP_DUPLICATED`. The retained Original Article is a measurement-aware study of selectively available reference labels in temporal ADMET model evaluation.

Working title:

**Selective reference measurements in temporal validation of ADMET prediction models: an estimand-aware evaluation and sensitivity study**

## Candidate comparison

| Candidate | Novelty | Identifiability/data | JPA fit | Decision |
|---|---|---|---|---|
| Generic AI-DMV plus AIRI and literature audit | Directly crowded by WelQrate, MOOD, ToxBench, TDC audits, and validation-framework work | Mixed datasets lack common time/source semantics; score weights arbitrary | Topic fit, article-type mismatch | STOP |
| Validation-design multiverse and model-selection regret | MOOD and current ExpansionRx comparisons already cover most of the broad question | Feasible, but internal ID time semantics are only pseudotime | Good | HOLD / secondary only |
| Local support × model-family superiority | Direct and adjacent prior work already demonstrates similarity-dependent model heterogeneity | Feasible | Good | STOP as main claim |
| Selectively available reference labels in later-stage ADMET evaluation | No directly isomorphic molecular/ADMET study found; methods exist in clinical/statistical domains | Near-complete endpoints allow known-truth validation; real sparse endpoints require bounded sensitivity analysis | Strong pharmaceutical-analysis and ADMET fit | GO_PILOT |

## Why the retained question is an Original Article

It has prespecified estimands, frozen predictions, controlled missing-label mechanisms with known truth, quantitative estimator validation, real campaign application, falsifiable gates, and negative-result rules. Literature is used only to define the gap and comparators.

## Conditions carried forward

- Use “reference availability” until blank-value causes are documented.
- Separate unavailable, censored/out-of-range, and invalid/QC states wherever the raw file permits.
- Never present MAR-adjusted estimates as identified under MNAR.
- Do not redistribute raw structures until the absent license file/attribution wording is clarified; publish retrieval code, pinned revision, hashes, and derived summaries instead.
- Update the collision search before manuscript freeze because the field is moving rapidly.
