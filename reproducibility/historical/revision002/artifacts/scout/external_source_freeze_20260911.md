# External source freeze for `confirm_001`

Audit date: 2026-09-11 (Asia/Shanghai)  
Scope: source identity, licensing, schema, parent-cohort construction, and preprocessing replay only. No model, risk estimator, or outcome analysis was run.

## Decision

**GO_DATA_READY**, conditional on the five frozen handling rules below. Both candidate sources are official, public, revision-pinned, locally staged without overwriting prior files, and sufficient for protocol construction. This is a data-readiness gate, not authorization to infer results from either dataset.

Machine-readable authority: `data/external/source_freeze_20260911/source_manifest.json`.

## 1. Biogen public ADME set

### Official source and license

- Official repository: <https://github.com/molecularinformatics/Computational-ADME>
- Current `main` revision observed in this audit: [`b00df003de117ce9e5b381afd886095c5f2af2d5`](https://github.com/molecularinformatics/Computational-ADME/commit/b00df003de117ce9e5b381afd886095c5f2af2d5), committed 2024-02-13 UTC.
- Repository license at that revision: MIT. The pinned `LICENSE` file has SHA-256 `06557ca6eeee99c162df4e5e819167c6103f35faffb51f7b0a4833905280a740`.
- Official data file: `ADME_public_set_3521.csv`, Git blob SHA-1 `a1fc97b07b02b22c4255d0d8c5a8efec878e8015`, 388,279 bytes, SHA-256 `2cfabc2667740c224487876c33b23124159ef43294e0f9e4d926cb6276c95a3b`.

The official README says this table contains 3,521 commercially sourced compounds and six experimental log-property endpoints: HLM, RLM, solubility, MDR1-MDCK efflux ratio, human plasma protein binding, and rat plasma protein binding.

### Exact source schema and completeness

The source is 3,521 rows × 10 columns, in this order:

`Internal ID`, `Vendor ID`, `SMILES`, `CollectionName`, `LOG HLM_CLint (mL/min/kg)`, `LOG MDR1-MDCK ER (B-A/A-B)`, `LOG SOLUBILITY PH 6.8 (ug/mL)`, `LOG PLASMA PROTEIN BINDING (HUMAN) (% unbound)`, `LOG PLASMA PROTEIN BINDING (RAT) (% unbound)`, `LOG RLM_CLint (mL/min/kg)`.

Non-missing endpoint counts are HLM 3,087; RLM 3,054; solubility 2,173; MDR1-MDCK 2,642; hPPB 194; and rPPB 168. All identifiers and SMILES are present. There are no duplicate `Internal ID` values and no duplicate exact SMILES strings.

### Source-to-derived mapping and exact split replay

The existing local derived file is `evidence/upstream/expansion_ml_comparison/biogen_adme_3521.csv`, SHA-256 `cd310dd3911a9f9763d0581e6202bc20aaf8a810523d01e9ab5e403681a386b4`. It retains all 3,521 source rows in the same order. The following eight fields were compared cell-by-cell, including missingness, with zero mismatches:

| Official field | Derived field |
|---|---|
| `Internal ID` | `Name` |
| `SMILES` | `SMILES` |
| solubility | `LOG_SOL` |
| HLM | `LOG_HLM` |
| RLM | `LOG_RLM` |
| MDR1-MDCK | `LOG_MDR1_ER` |
| human PPB | `LOG_HPPB` |
| rat PPB | `LOG_RPPB` |

`Vendor ID` and `CollectionName` are intentionally not carried into the derived table. No molecule was removed for invalid SMILES or all-endpoint missingness.

The local preparation chain belongs to the separate `PatWalters/expansion-ml-comparison` repository at revision `06582771d48af8096f6b05a5a247b177189a8fa4`; it is not part of the official Biogen repository. The relevant frozen file hashes are recorded in the JSON manifest. Its construction is:

1. Generate 2,048-bit ECFP4 fingerprints and BitBIRCH-Lean clusters using the frozen `clustering.py` settings.
2. Compute cluster sizes, permute the `value_counts()` index with `numpy.random.default_rng(61453)`, and hold out entire clusters until at least 30% of rows are in test.
3. The resulting table has 1,905 clusters, 1,123 singletons, 555 held-out clusters, 2,463 train rows, and 1,058 test rows; no cluster straddles the split.

An independent clean replay under the repository-declared versions—Python 3.12, numpy 2.5.2, pandas 3.0.5, scipy 1.17.1, scikit-learn 1.9.0, RDKit 2026.03.5, and bblean 0.10.3—matched both the full cluster vector and the full `ds` vector exactly (0 mismatched rows).

Important reproducibility finding: replaying only the split under the workstation's older pandas 2.3.3 / numpy 1.26.4 did **not** select the same tied-size clusters. This is an environment-ordering effect, not a source-data mismatch. `confirm_001` must either consume the frozen derived table/vector or run the exact pinned preprocessing environment; silently regenerating the split in another environment is prohibited.

## 2. OpenADMET PXR assay cascade

### Official source, revision, license, and inventory

- Official dataset: <https://huggingface.co/datasets/openadmet/pxr-challenge-train-test>
- DOI: <https://doi.org/10.57967/hf/9731>
- Current revision observed in this audit: `8041adfc295cac7f0bf85b73a2d87f21c5302ab1`, last modified 2026-09-02 10:40:11 UTC.
- License displayed in the official dataset metadata and README front matter: CC BY 4.0. The repository has no standalone license file.
- Inventory: 196 files—12 top-level files plus 184 PDB files under `structure_ground_truth/`. The 12 top-level names are frozen in the JSON manifest.

Only three files were retrieved because they are the minimum needed to define the proposed selection-process confirmation: the pinned README, primary dose-response training table, and single-concentration parent-screen table. No blinded test tables, counter-assay table, HT-chem data, or structural files were downloaded.

| File | Shape | Bytes | SHA-256 |
|---|---:|---:|---|
| `pxr-challenge_TRAIN.csv` | 4,139 × 18 | 858,687 | `fdc2252bb5763ab5742e78cb1ba692bebe131a972cf83aa1bd45e43d8c5bc223` |
| `pxr-challenge_single_concentration_TRAIN.csv` | 21,003 × 19 | 6,290,720 | `e781f8405e919a248b9b832a7cb24e8a615d400ee9dc50beea0d70d8aeaa8c2f` |
| `README.md` | — | 6,311 | `8c695adef4e113d8a59ae2c20fcfc3a5f30a82fb6bd6cd89a141a531d08995d0` |

The DRC table has complete values in all 18 fields and 4,139 unique `OCNT_ID` and exact SMILES values. The single-concentration table has 10,870 unique `OCNT_ID`; only `Molecule Name` is incomplete (15,237 blank rows), so it must not be used as a join key.

### Documented screen-to-dose-response selection

The primary OpenADMET assay description documents a two-tier design: broad single-concentration screening; formal advancement to eight-point dose-response when `log2 fold change > 1` and Benjamini-Hochberg `FDR < 0.05` at the lowest nominal 10 µM screen; additional analog and diversity-library follow-up; and downstream QC. Primary source: <https://openadmet.ghost.io/announcing-the-next-openadmet-blind-challenge-predicting-pxr-induction/>.

The file stores the lowest library concentration as exactly `8.251e-06 M` (8.251 µM), although the source narrative calls it nominally 10 µM. The frozen parent cohort must therefore use the data value, not a rounded equality test:

`compound_class == "Library" AND concentration_M == 8.251e-06`

This produces 10,747 rows, each with a unique `OCNT_ID`, exact SMILES, and batch. Join dose-response availability by `OCNT_ID`:

- 2,722 / 10,747 parent compounds have a DRC record.
- The formal hit rule identifies 873 rows; 869 have a DRC record and 4 do not.
- A further 1,853 parent rows have a DRC record without satisfying the formal rule at this concentration, consistent with the source's additional selection routes and QC, but the files do not provide a row-specific reason.

`OCNT_ID` is the only defensible cross-table join key. Joining exact SMILES finds only 2,375 parent DRC matches and joining exact `OCNT Batch` finds 2,700. Across the two tables, 351 `OCNT_ID` values have more than one exact SMILES representation and 23 have more than one batch value. Those differences must not be mistaken for different compounds.

Scientific boundary: DRC availability is a documented **mixed-route, outcome-informed advancement indicator**. It is not a binary PXR activity label. Lack of a DRC record means unavailable follow-up or non-advancement in these files; it must never be relabeled as inactive, negative, assay failure, or a random missing value.

## Frozen handling rules

1. Pin every future run to the recorded repository revisions and file SHA-256 values.
2. For Biogen, use the frozen derived table or the exact pinned preprocessing environment. Do not regenerate `ds` with different dependency versions.
3. For PXR, use the exact `8.251e-06 M` parent filter and join stages by `OCNT_ID`.
4. Treat PXR availability as mixed-route and outcome-informed; any MNAR statement is a design boundary or sensitivity claim, not identified truth.
5. Keep all staged source files local. A public package may contain retrieval instructions, revisions, hashes, and required attribution, but not these raw CSVs or upstream README/LICENSE copies.

## Local staging result

New directory created without overwriting prior files: `data/external/source_freeze_20260911/`.

Contents are restricted to the two Biogen provenance files plus source CSV and the three minimal PXR files listed above, together with the machine-readable manifest. No modeling output was written.
