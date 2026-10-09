# Release scope and provenance

This repository publishes the analysis-code subset of `CILS_format_revision_candidate_20261009.zip`, whose SHA-256 is `fe480776593996333bf961488fee5ceb6e0c10ee255cb1266936ee238b3b8d96`.

The original local archive is retained unchanged. The repository places the selected files under `reproducibility/` so Git metadata and repository documentation stay outside the exact payload inventory. The analysis scripts, parameters, saved estimates, electronic tables and figure bytes retain their source-archive contents. Language edits concern the current README and run instructions. Upstream code and notices are preserved.

## Included and omitted material

The payload contains estimator and refitting code, saved estimates and summaries, figure and table inputs, protocols, input locators and selected historical verification records. It omits manuscript text, cover letters, author metadata, submission responses, a draft citation template, the historical Word builder and the author-workspace package builder. `reproducibility/PUBLICATION_SCOPE.json` records every omitted path and its source hash, along with changed and retained paths.

The root payload manifest and the historical revision-002 inventory were regenerated for this subset. Original execution manifests, code hashes and dated status records retain their historical meaning. A reference to an omitted input records a reconstruction requirement; it does not imply that the file is distributed. The historical standalone release checker needs the complete earlier candidate, whereas the current root saved-result checker verifies the subset published here.

Raw assays, SMILES and scaffolds, full compound-level predictions, PXR linkage records, local feature caches and row-vector NPZs remain excluded. The permitted Mol2946 locator, limited case arithmetic and split/source-row indices are retained. The repository also excludes author Word files, user attachments, local environments and earlier ZIP archives.

## Licenses and study status

The repository's existing MIT license is unchanged. It covers project software and associated documentation. Third-party software retains its original notices, including the upstream expansion-ml-comparison notice. Software licensing does not grant rights to redistribute third-party assay data; source terms and access instructions remain in `historical/revision002/DATA_ACCESS.md` under the payload. This publication does not assign a new license to the aggregate results.

Historical README and status files describe their dated local candidates. Their statements about access, manuscripts and submission approval apply to those earlier deliveries. This repository provides access to the code subset; it does not create a DOI, complete author declarations, establish journal acceptance or document formal submission.

The publication checks repeat saved-result arithmetic and file integrity. They do not establish complete independent model retraining or external validation.
