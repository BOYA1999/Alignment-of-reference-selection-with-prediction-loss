# Verified reference audit

Date checked: 2026-09-12

## Scope and selection

`references.bib` and `references_verified.csv` contain 33 matched records: 23 journal/proceedings papers, 4 explicitly labelled preprints, 2 official datasets, 1 official challenge announcement, 1 versioned software release, 1 primary software issue and 1 official descriptor documentation page. The original 28 core and three recommended records are retained; two primary software sources were added for revision 002. The set is limited to sources needed by this Original Article.

The four preprints are Seal et al. (bioRxiv 2026), Lim (arXiv 2026), Xue et al. (arXiv 2026), and Ulichney and Coston (arXiv 2026). Corbin et al. is now cited as the published AMIA 2023 paper (PMID 37350883, pages 81–90), replacing its older preprint record. The September 2026 review-linked primary-source audit also verified Li, Chen/Qin, Wyss, Liu/Kashima, Kretschmer and Ulichney/Coston; see `rebuttal/pxr_audit.md`. Li's missing-covariate setting is distinguished from this study's masked-outcome setting.

## Verification basis

- Journal and proceedings metadata were checked against DOI/Crossref records; publisher or PubMed records were used where issue-level publication details required clarification.
- Dataset, Zenodo, software-release, and arXiv metadata were checked against DataCite and the official repository or preprint record.
- The official JMLR article page and BibTeX record were used for Pedregosa et al.; no DOI is assigned there.
- Full author sequences are retained, including all 29 authors for Seal et al., all five authors for Xue et al., and all 30 creators on the RDKit 2026.03.3 release.
- Petersen et al. is cited as the 2012 volume/issue publication; its online-first date was 2010. Kretschmer et al. is cited as 2025 even though its DOI contains `2024`.

## Data-version boundary

The dataset DOI records identify the evolving repositories, not the exact analytical snapshots. The study pins ExpansionRx revision `6b898ccc43d10d25b230fb09e22a6e30c30022b5` and PXR revision `8041adfc295cac7f0bf85b73a2d87f21c5302ab1`. These full revisions and the already frozen file hashes should appear in Methods/Data Availability; they must not be replaced by a short version currently displayed on a mutable DOI landing page.

## Recommended manuscript roles

- Data and ADMET context: Fang; official ExpansionRx and PXR datasets; the PXR announcement.
- Direct novelty boundary: Seal; Lim; Xue; Ulichney and Coston; Corbin; Coston's published FAccT paper.
- Molecular bias and validation context: Liu; Kretschmer; Ash.
- Selective-label and fixed-model evaluation: Lakkaraju; Corbin; Li. Li concerns missing predictors/covariates and is only an adjacent precedent, not an exact missing-outcome antecedent.
- AIPW, double robustness, cross-fitting, and diagnostics: Robins; Bang and Robins; Kang and Schafer; Chernozhukov et al.; Cole and Hernán; Petersen et al.; Austin and Stuart; Chen and Qin.
- Falsification/admission logic: Chen and Qin; Coston; Wyss et al. Wyss is an analogy from causal pharmacoepidemiology, not direct ADMET evidence.
- Molecular implementation: Bemis and Murcko; Morgan; Rogers and Hahn; the exact RDKit release; Breiman; Chen and Guestrin; Pedregosa et al. Model/software citations should remain only if the corresponding components occur in the executed analysis.

## Claim constraints

Seal et al. closes a broad claim of a new general ADME validation framework. Ulichney and Coston closes generic novelty for fixed-predictor risk estimation under selective labels using influence-function/DML machinery. Xue et al. closes claims of first combining Biogen and ExpansionRx or first handling their incomplete labels. Lim establishes follow-up quantifiability as a prediction/triage target. The defensible contribution is therefore the preregistered ADMET-specific falsification/admission experiment, its negative correction result, and the resulting PXR abstention boundary—not invention of IPW/AIPW, DML, selective-label evaluation, or general ADMET validation.

## Exclusions

The ExpansionRx announcement, active-testing papers, MoleculeNet, time-split and memorization papers were left out to keep the bibliography focused. Martin and Bowen's ChemRxiv item was excluded because the reported DOI could not be independently resolved during this audit and its unknown-interaction-label problem is not the same as selectively unavailable ADMET reference outcomes. OpenADMET blog posts and the PXR tutorial repository are documentation rather than primary scientific evidence and should be cited only in Data/Code Availability if needed.

## Artifact validation

A strict read-only BibTeX structure parser accepted all 31 entries and required fields. The BibTeX and CSV key sets agree exactly, with 31 unique keys, 30 unique DOI values, and no duplicate fields or DOI values; the sole no-DOI record is the official JMLR scikit-learn article. All 30 DOI values were resolved against Crossref or DataCite on the audit date. No local BibTeX, Biber, Pybtex, or bibtexparser executable/library was available, so bibliography-style compilation remains a later manuscript-build check rather than part of this artifact audit.
