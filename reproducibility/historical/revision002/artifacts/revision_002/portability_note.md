# Post-execution packaging amendment

The numerical runs and independent checks completed before this packaging-only edit. Their executed source files are retained in executed_source/, and the recorded hashes in original/manifest.json, bridge/manifest.json and verification.json identify those executed bytes.

The distributed revision2_refit.py adds an explicit optional parent-archive verification switch. It remains enabled by default in the author workspace. The staged public wrapper disables only that archive check because the previous archive is not a scientific input and is not nested in the new archive; generated records disclose parent_archive_locally_verified=false. Frozen scientific inputs are still hashed before and after every phase. No estimator, seed, feature, fit, or metric expression changed.

verify_revision2.py exposes the existing cache path as a module variable so isolated reconstructed inputs can be checked. The staged entry first requires the 19-file confirmation and checksum-matching forensic cache. It never overwrites the distributed evidence or author workspace. Its prerequisite-failure and path redirection were checked locally; a new full public Stage 1 plus Stage 2 execution has not been completed in this revision.
