import csv
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_NAMES = {'predictions.parquet', 'nuisance_fold_assignments.parquet', 'generator_vectors.parquet',
                   'row_replay_seed41001.parquet', 'row_replay_seed31001.parquet', 'pxr_linkage_audit.csv',
                   'ADME_public_set_3521.csv', 'biogen_adme_3521.csv', 'biogen_adme_reconstructed_local_only.csv',
                   'expansion_data_train.csv', 'expansion_data_test.csv', 'expansion_data_raw.csv',
                   'pxr-challenge_TRAIN.csv', 'pxr-challenge_single_concentration_TRAIN.csv',
                   'top_influence_rows.csv', 'numerical_refit_predictions.csv', 'ipc_coefficient_contributions.csv',
                   'local_feature_cache.npz', 'mechanism_calibration_replicates.parquet', 'pasted-text.txt'}
ROW_KEYS = {'smiles', 'canonical_smiles', 'molecule_id', 'ocnt_id', 'source_row_index', 'cohort_row_index',
            'scaffold_group', 'max_z_source_row', 'dominant_source_row_index', 'dominant_molecule_id'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_strings(value):
    if isinstance(value, dict):
        for key, item in value.items():
            assert key.lower() not in {'dominant_source_row_index', 'dominant_molecule_id'}
            yield from json_strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from json_strings(item)
    elif isinstance(value, str):
        yield value


def main():
    if not __debug__:
        raise RuntimeError('Run package verification without Python assertion checks disabled')
    manifest = {}
    for line in (ROOT / 'MANIFEST.sha256').read_text(encoding='utf-8').splitlines():
        value, name = line.split('  ', 1)
        assert name not in manifest and len(value) == 64 and '..' not in Path(name).parts
        manifest[name] = value
    paths = [p for p in ROOT.rglob('*') if p.is_file() and p.name != 'MANIFEST.sha256']
    assert {p.relative_to(ROOT).as_posix() for p in paths} == set(manifest)
    for path in paths:
        relative = path.relative_to(ROOT).as_posix()
        assert sha(path) == manifest[relative], relative
        assert path.name not in FORBIDDEN_NAMES and not any(part in {'.venv', '__pycache__', 'input', 'backups'} for part in path.parts)
        assert not path.name.endswith(('.pyc', '.npz', '.parquet'))
        if path.suffix == '.csv':
            with path.open(encoding='utf-8-sig', newline='') as stream:
                columns = next(csv.reader(stream))
            assert not ROW_KEYS.intersection(c.lower() for c in columns), relative
        if path.suffix == '.json':
            for text in json_strings(json.loads(path.read_text(encoding='utf-8'))):
                assert not re.search(r'(?<![A-Za-z0-9])[A-Za-z]:[\\/]', text), relative
                assert 'PRIVATE KEY-----' not in text, relative
        elif path != Path(__file__).resolve() and path.suffix in {'.md', '.txt', '.py', '.csv', '.svg'}:
            text = path.read_text(encoding='utf-8', errors='strict')
            assert not re.search(r'(?i)[A-Z]:[\\/]+Users[\\/]', text), relative
            assert 'PRIVATE KEY-----' not in text and 'api_key=' not in text.lower(), relative
            if path.suffix in {'.md', '.csv', '.svg'}:
                assert not re.search(r'\bMol\d+\b|\bsource row \d+\b', text), relative
    contract = json.loads((ROOT / 'release_contract.json').read_text())
    redactions = json.loads((ROOT / 'verification/distribution_redactions.json').read_text())
    assert all(sha(ROOT / row['path']) == row['distributed_sha256'] for row in redactions)
    for name, expected in contract['protected_original_code_sha256'].items():
        assert sha(ROOT / name) == expected
    figures = {p.name: sha(p) for p in (ROOT / 'paper/figures').iterdir() if p.is_file()}
    assert figures == contract['figure_files'] and len(figures) == 24
    assert {Path(name).suffix for name in figures} == {'.svg', '.pdf', '.png', '.tiff'}
    assert {p.name: sha(p) for p in (ROOT / 'paper/tables').glob('*.csv')} == contract['table_files']
    base = ROOT / 'verification/accepted_aggregate_outputs'
    original = pd.read_csv(base / 'factorial_contrasts.csv')
    assert len(original) == 56 and int(original.contrast.eq('P_x_M').sum()) == 16
    risk = ROOT / 'artifacts/revision_001/risk'
    cells = pd.read_csv(risk / 'original_candidate_16_cells.csv')
    assert len(cells) == 16 and int(cells.loc[cells.mechanism.eq('mar'), 'error_ratio_pass'].sum()) == 0
    assert int(cells.loc[cells.mechanism.eq('mar'), 'all_support_diagnostics_pass'].sum()) == 2
    calibration = pd.read_csv(risk / 'mechanism_calibration_summary.csv')
    assert len(calibration) == 192 and calibration.replicates.eq(10000).all()
    table_sources = {
        'Table_S8_original_estimator_summary.csv': 'risk/original_estimator_summary.csv',
        'Table_S8_candidate_gate_chain.csv': 'risk/original_candidate_16_cells.csv',
        'Table_S9_mechanism_calibration.csv': 'risk/mechanism_calibration_summary.csv',
        'Table_S10_forensic_estimates.csv': 'forensic/stable_sensitivity_estimates.csv',
        'Table_S11A_propensity_tails.csv': 'pxr/propensity_tails.csv',
        'Table_S11B_weighting_sensitivity.csv': 'pxr/weighting_sensitivity.csv',
        'Table_S11C_parameter_completeness.csv': 'pxr/parameter_completeness.csv',
        'Table_S12_nonminimum_sensitivity.csv': 'risk/nonminimum_sensitivity.csv',
    }
    for name, source in table_sources.items():
        left = pd.read_csv(ROOT / 'paper/tables' / name, float_precision='round_trip')
        right = pd.read_csv(ROOT / 'artifacts/revision_001' / source, float_precision='round_trip')
        if name == 'Table_S8_candidate_gate_chain.csv':
            right = right.sort_values(['endpoint', 'anchor', 'model', 'mechanism']).reset_index(drop=True)
        pd.testing.assert_frame_equal(left, right, check_dtype=False, atol=1e-12, rtol=0)
    figure4_source = pd.read_csv(ROOT / 'paper/source_data/figure4_fixed_nuisance_calibration.csv', float_precision='round_trip')
    full_calibration = pd.read_csv(risk / 'mechanism_calibration_summary.csv', float_precision='round_trip')
    pd.testing.assert_frame_equal(figure4_source, full_calibration[full_calibration.method.eq('true_p_fixed_replay_m_aipw')].reset_index(drop=True),
                                  check_dtype=False, atol=1e-12, rtol=0)
    stable = calibration[calibration.method.eq('true_p_fixed_replay_m_aipw')]
    assert stable.groupby('mechanism').original_error_component_pass.sum().to_dict() == {
        'loss_rank_beta_0.5': 8, 'loss_rank_beta_2.0': 8, 'original_mar': 0, 'original_mcar': 8}
    p1 = pd.read_csv(risk / 'complete_observation_variance_check.csv')
    assert len(p1) == 4 and p1.true_mask_se.eq(0).all() and p1.legacy_centered_pseudo_cluster_se.gt(0).all()
    reconstruction = json.loads((ROOT / 'artifacts/revision_001/reconstruction/reconstruction_report.json').read_text())
    assert reconstruction['csv_bytes_equal'] and reconstruction['accepted_vector_hashes_match']
    assert all(v == 0 for v in reconstruction['column_mismatch_rows'].values())
    print(f'PASS revised local release: {len(paths)} manifest entries, 24 figures, protected original code unchanged')


if __name__ == '__main__':
    main()
