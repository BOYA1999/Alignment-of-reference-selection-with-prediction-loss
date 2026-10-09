import argparse
import hashlib
import json
import shutil
from pathlib import Path
import revision2_refit as fit
import revision2_alignment as alignment
import verify_revision2 as verify


def main():
    parser = argparse.ArgumentParser(description='Stage 2: revision refits after the separately reconstructed confirmation and forensic cache')
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    frozen = root/'public_replay/confirmation'
    cache = root/'public_replay/revision_diagnostics/revision_forensic'
    required = [frozen/x for x in ['predictions.parquet','row_replay_seed41001.parquet','generator_vectors.parquet',
                                  'semisynthetic_scenarios.csv','structure_audit.json']]
    required += [cache/'local_feature_cache.npz', cache/'refit_manifest.json']
    missing = [str(p.relative_to(root)) for p in required if not p.is_file()]
    assert not missing, f'First run fetch_public_inputs.py, then replay_confirmation.py --diagnostics. Missing: {missing}'
    assert len(list(frozen.iterdir())) == 19
    assert not (root/'dingzhen.md').exists(), 'Use a fresh extracted package, not the author workspace'
    manifest = json.loads((cache/'refit_manifest.json').read_text())
    assert fit.sha(cache/'local_feature_cache.npz') == manifest['output_sha256']['local_feature_cache.npz']
    report = {'stage':'revision_002','prerequisites_present':True,'original_frozen_file_count':19,
              'cache_hash_verified':True,'full_model_refit_executed':False}
    if args.check_only:
        print(json.dumps(report,indent=2))
        return
    out = root/'public_replay/revision_002'
    out.mkdir(exist_ok=False)
    shutil.copy2(root/'artifacts/revision_002/protocol.md',out/'protocol.md')
    fit.FROZEN, fit.BASE, fit.OUT = frozen, cache, out
    fit.run('original',False,verify_parent_archive=False)
    fit.run('bridge',False,verify_parent_archive=False)
    alignment.SOURCE, alignment.OUT = frozen, out/'alignment'
    with alignment.threadpool_limits(limits=1):
        alignment.main()
    verify.OLD, verify.BASE, verify.CACHE = frozen, out, cache/'local_feature_cache.npz'
    with verify.threadpool_limits(limits=1):
        verify.main()
    report.update(status='COMPLETED_REVISION_REFITS_AND_SPECIFIED_CHECKS',full_model_refit_executed=True,
                  scope='Stage 2 refits auxiliaries only; fixed RF/XGB predictions come from prior Stage 1',
                  parent_archive_locally_verified=False)
    (out/'replay_stage_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
