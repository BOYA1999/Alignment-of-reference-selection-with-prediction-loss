import argparse
import hashlib
import urllib.request
from pathlib import Path


FILES = [
    (
        "data/raw/expansion_data_train.csv",
        "https://huggingface.co/datasets/openadmet/openadmet-expansionrx-challenge-data/resolve/6b898ccc43d10d25b230fb09e22a6e30c30022b5/expansion_data_train.csv",
        "c5214ad8c8a4d7d4d09082fcce24fe16c97a980a1e0873b91b8c8e474b79f6e4",
    ),
    (
        "data/raw/expansion_data_test.csv",
        "https://huggingface.co/datasets/openadmet/openadmet-expansionrx-challenge-data/resolve/6b898ccc43d10d25b230fb09e22a6e30c30022b5/expansion_data_test.csv",
        "e21c0ef19795317216287e924e00ca3ef9388b7032a7fa4c7f27e2580372c85e",
    ),
    (
        "data/raw/expansion_data_raw.csv",
        "https://huggingface.co/datasets/openadmet/openadmet-expansionrx-challenge-data/resolve/6b898ccc43d10d25b230fb09e22a6e30c30022b5/expansion_data_raw.csv",
        "f674ec74cca1146bc386f832a32d4b8d921d3c312f92cb436cc005901c724a3c",
    ),
    (
        "data/external/source_freeze_20260911/biogen/ADME_public_set_3521.csv",
        "https://raw.githubusercontent.com/molecularinformatics/Computational-ADME/b00df003de117ce9e5b381afd886095c5f2af2d5/ADME_public_set_3521.csv",
        "2cfabc2667740c224487876c33b23124159ef43294e0f9e4d926cb6276c95a3b",
    ),
    (
        "data/external/source_freeze_20260911/biogen/LICENSE",
        "https://raw.githubusercontent.com/molecularinformatics/Computational-ADME/b00df003de117ce9e5b381afd886095c5f2af2d5/LICENSE",
        "06557ca6eeee99c162df4e5e819167c6103f35faffb51f7b0a4833905280a740",
    ),
    (
        "data/external/source_freeze_20260911/biogen/README.md",
        "https://raw.githubusercontent.com/molecularinformatics/Computational-ADME/b00df003de117ce9e5b381afd886095c5f2af2d5/README.md",
        "056e4811aa4c0e6ca096b8d1ae5d33a5d91f3ba02d1aad1ebce0d779c67be425",
    ),
    (
        "evidence/upstream/expansion_ml_comparison/biogen_adme_3521.csv",
        "https://raw.githubusercontent.com/PatWalters/expansion-ml-comparison/06582771d48af8096f6b05a5a247b177189a8fa4/biogen_adme_3521.csv",
        "cd310dd3911a9f9763d0581e6202bc20aaf8a810523d01e9ab5e403681a386b4",
    ),
    (
        "data/external/source_freeze_20260911/pxr/README.md",
        "https://huggingface.co/datasets/openadmet/pxr-challenge-train-test/resolve/8041adfc295cac7f0bf85b73a2d87f21c5302ab1/README.md",
        "8c695adef4e113d8a59ae2c20fcfc3a5f30a82fb6bd6cd89a141a531d08995d0",
    ),
    (
        "data/external/source_freeze_20260911/pxr/pxr-challenge_TRAIN.csv",
        "https://huggingface.co/datasets/openadmet/pxr-challenge-train-test/resolve/8041adfc295cac7f0bf85b73a2d87f21c5302ab1/pxr-challenge_TRAIN.csv",
        "fdc2252bb5763ab5742e78cb1ba692bebe131a972cf83aa1bd45e43d8c5bc223",
    ),
    (
        "data/external/source_freeze_20260911/pxr/pxr-challenge_single_concentration_TRAIN.csv",
        "https://huggingface.co/datasets/openadmet/pxr-challenge-train-test/resolve/8041adfc295cac7f0bf85b73a2d87f21c5302ab1/pxr-challenge_single_concentration_TRAIN.csv",
        "e781f8405e919a248b9b832a7cb24e8a615d400ee9dc50beea0d70d8aeaa8c2f",
    ),
]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    for relative, url, expected in FILES:
        path = root / relative
        if not path.is_file() or digest(path) != expected:
            if args.check_only:
                raise RuntimeError(f"missing or mismatched input: {relative}")
            path.parent.mkdir(parents=True, exist_ok=True)
            urllib.request.urlretrieve(url, path)
        if relative == "evidence/upstream/expansion_ml_comparison/biogen_adme_3521.csv":
            data = path.read_bytes()
            if b"\r\n" not in data:
                path.write_bytes(data.replace(b"\n", b"\r\n"))
        observed = digest(path)
        if observed != expected:
            raise RuntimeError(f"hash mismatch for {relative}: {observed}")
        print(f"PASS {relative} {observed}")


if __name__ == "__main__":
    main()
