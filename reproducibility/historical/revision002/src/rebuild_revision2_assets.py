import argparse
import hashlib
import json
import shutil
from pathlib import Path
import build_revision2_assets as build


def main():
    parser = argparse.ArgumentParser(description='Rebuild revision-002 figures/tables in an empty output directory')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    for name in ['figures','tables','source_data']:
        (out/'paper'/name).mkdir(parents=True)
    inputs = ['original/estimates.csv','bridge/summary.csv','bridge/design.csv','alignment/summary.csv',
              'alignment/design.csv','independent_solver.csv']
    for name in inputs:
        destination = out/'inputs'/name
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(root/'artifacts/revision_002'/name,destination)
    for source in (root/'paper/figures').glob('Figure_S*'):
        if source.suffix in ['.png','.pdf','.svg']:
            shutil.copy2(source,out/'paper/figures'/source.name)
    build.REV, build.PAPER = out/'inputs', out/'paper'
    build.main()
    checks = {}
    for path in (out/'paper').rglob('*'):
        if path.is_file():
            rel = path.relative_to(out)
            checks[rel.as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest() == hashlib.sha256((root/rel).read_bytes()).hexdigest()
    assert all(checks.values()), [p for p,v in checks.items() if not v]
    report = {'status':'PASS','files_compared':len(checks),'byte_matches':sum(checks.values()),
              'regenerated_figure_files':16,'retained_supplementary_figure_copies':12,
              'scope':'Regenerates current main figures, graphical abstract and revision-002 tables/source data; historical supplementary figures are copied, not regenerated',
              'checks':checks}
    (out/'asset_replay_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))


if __name__ == '__main__':
    main()
