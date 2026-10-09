import argparse
import hashlib
import importlib.util
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.dont_write_bytecode = True


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    if not __debug__:
        raise RuntimeError('Do not run scientific reconstruction with Python assertion checks disabled')
    parser = argparse.ArgumentParser(description='Post-result public scientific replay, separate from historical authorization.')
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--diagnostics', action='store_true')
    args = parser.parse_args()
    root = args.root.resolve()
    historical = root / 'verification/accepted_aggregate_outputs/run_manifest.json'
    if not historical.exists() and args.check_only:
        historical = root / 'artifacts/confirmatory/confirm_001/full_recovery_001/run_manifest.json'
    original = json.loads(historical.read_text(encoding='utf-8'))
    expected_aggregates = {name: value for name, value in original['output_hashes_before_manifest'].items() if name.endswith('.csv')}
    assert expected_aggregates, 'Historical aggregate inventory is empty'
    assert all((historical.parent / name).is_file() and digest(historical.parent / name) == value
               for name, value in expected_aggregates.items()), 'Historical aggregate comparator missing or changed'
    names = ['src/run_confirmatory.py', 'src/verify_confirmatory.py',
             'artifacts/confirmatory/confirm_001/protocol.md', 'artifacts/confirmatory/confirm_001/config.json',
             'data/external/source_freeze_20260911/source_manifest.json',
             'data/external/source_freeze_20260911/source_manifest_amendment_001.json']
    fetch = load_module(root / 'src/fetch_public_inputs.py', 'public_input_inventory')
    expected = {name: original['input_hashes'][name] for name in names}
    expected.update({name: value for name, _, value in fetch.FILES})
    cfg = json.loads((root / names[3]).read_text(encoding='utf-8'))
    construction = cfg['source_freeze']['biogen']['split_provenance']['construction_files_sha256']
    expected.update({f'evidence/upstream/expansion_ml_comparison/{name}': value for name, value in construction.items()})
    missing = [name for name in expected if not (root / name).is_file()]
    assert not missing, f'Retrieve pinned public inputs first; missing {missing}'
    actual = {name: digest(root / name) for name in expected}
    assert actual == expected, 'Scientific input or protected original code hash mismatch'
    runner = load_module(root / 'src/run_confirmatory.py', 'protected_confirmatory_scientific_functions')
    versions = {'python': platform.python_version(), 'numpy': runner.np.__version__, 'pandas': runner.pd.__version__,
                'scipy': runner.scipy.__version__, 'scikit_learn': runner.sklearn.__version__,
                'rdkit': runner.rdkit.__version__, 'xgboost': runner.xgboost.__version__,
                'pyarrow': runner.pyarrow.__version__, 'platform': platform.platform()}
    required = {k: v for k, v in cfg['analysis_environment'].items() if k not in ['platform', 'split_provenance_environment_is_separate']}
    assert {k: versions[k] for k in required} == required, f'Use the analysis dependency pins: {versions}'
    report = {'status': 'PASS_SCIENTIFIC_REPLAY_PREFLIGHT', 'scientific_inputs_checked': len(actual),
              'historical_aggregate_csv_files_checked': len(expected_aggregates),
              'wrapper_sha256': digest(Path(__file__)),
              'input_sha256': actual, 'environment': versions,
              'historical_platform_matches': versions['platform'] == cfg['analysis_environment']['platform'],
              'scope': 'Checks required scientific inputs and code. Does not replay historical author authorization or execute any model fit.'}
    if args.check_only:
        print(json.dumps(report, indent=2))
        return
    assert (root / 'RELEASE_STATUS.md').is_file(), 'Full replay must run in a fresh extracted release package'
    assert not (root / 'dingzhen.md').exists(), 'Never execute full replay in the author workspace'
    assert not (root / 'artifacts/confirmatory/confirm_001/full_recovery_001').exists(), 'Preserve any historical result directory'
    output = root / 'public_replay/confirmation'
    assert not output.parent.exists(), 'Use a fresh package extraction; preserve prior replay output'
    output.mkdir(parents=True)
    amendment = json.loads((root / names[5]).read_text(encoding='utf-8'))
    started = time.time()
    runner.run_full(root, output, cfg, amendment, actual, versions, started,
                    datetime.now(timezone.utc).isoformat(), 'PUBLIC_POST_RESULT_SCIENTIFIC_REPLAY',
                    ['frozen scientific protocol and inputs; historical author workflow is not re-executed'],
                    [], {'base_protocol': names[2], 'base_config': names[3], 'source_manifest': names[4],
                         'source_manifest_amendment': names[5]})
    comparison = {}
    assert {path.name for path in output.glob('*.csv')} == set(expected_aggregates), 'Replayed aggregate inventory mismatch'
    for name in sorted(expected_aggregates):
        left = runner.pd.read_csv(historical.parent / name, float_precision='round_trip')
        right = runner.pd.read_csv(output / name, float_precision='round_trip')
        try:
            runner.pd.testing.assert_frame_equal(left, right, check_dtype=False, atol=1e-12, rtol=0)
            comparison[name] = 'MATCH_WITHIN_1e-12'
        except AssertionError as error:
            comparison[name] = 'MISMATCH: ' + str(error)[:250]
    report.update(status='CONFIRMATION_REPLAY_COMPLETE_DIAGNOSTICS_PENDING' if args.diagnostics else 'COMPLETED_POST_RESULT_REPLAY_NOT_INDEPENDENT_VALIDATION',
                  scope='Executed the original scientific predictor and nuisance fits into isolated post-result output; historical author authorization was not replayed.',
                  output_directory='public_replay/confirmation', aggregate_comparison=comparison,
                  every_compared_aggregate_matches=all(v == 'MATCH_WITHIN_1e-12' for v in comparison.values()),
                  elapsed_seconds=time.time() - started,
                  historical_authorization_chain_replayed=False, diagnostics_requested=args.diagnostics,
                  diagnostics={})

    def write_report():
        (output.parent / 'replay_report.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')

    write_report()
    assert report['every_compared_aggregate_matches'], 'Replayed results differ; inspect preserved comparison report before diagnostics'
    if args.diagnostics:
        destination = output.parent / 'revision_diagnostics'
        for name in ['revision_risk_diagnostics', 'revision_pxr_audit', 'revision_forensic']:
            report['diagnostics'][name] = 'RUNNING'
            write_report()
            try:
                module = load_module(root / f'src/{name}.py', f'public_{name}')
                module.ROOT = root
                module.OUT = destination / name
                module.OUT.mkdir(parents=True)
                if name == 'revision_risk_diagnostics':
                    module.SOURCE = output
                    module.main()
                else:
                    module.FROZEN = output
                    if name == 'revision_forensic':
                        with module.threadpool_limits(limits=2):
                            module.inspect()
                            module.refit()
                    else:
                        module.main()
            except Exception as error:
                report['status'] = 'DIAGNOSTIC_FAILED_CONFIRMATION_OUTPUT_RETAINED'
                report['diagnostics'][name] = 'FAILED: ' + str(error)[:250]
                write_report()
                raise
            report['diagnostics'][name] = 'COMPLETED'
            write_report()
        report['status'] = 'COMPLETED_REPLAY_AND_REQUESTED_DIAGNOSTICS_NOT_INDEPENDENT_VALIDATION'
        report['elapsed_seconds'] = time.time() - started
        write_report()
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
