import argparse
import hashlib
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.special import expit, logit
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler


ENDPOINT_COLUMNS = {
    "LogD": "LogD",
    "KSOL": "KSOL",
    "HLM": "HLM CLint",
    "MLM": "MLM CLint",
    "Caco_Papp": "Caco-2 Permeability Papp A>B",
    "Caco_Efflux": "Caco-2 Permeability Efflux",
    "MPPB": "MPPB",
    "MBPB": "MBPB",
    "MGMB": "MGMB",
}
ALL_MODELS = ["median", "ecfp_ridge", "ecfp_rf", "ecfp_desc_xgb"]
TRUSTED_RUNNER_SHA256 = "C54C839B4363E3134AB9C8ACF05DB53822A66BF087CD5EC9810C84AB27EFC6CF"
TRUSTED_ESTIMATORS_SHA256 = "79DD68DA10DC85EF85DE3FFF18430D53392D900558DC3FC9A54522E74E2649E2"

warnings.filterwarnings(
    "ignore",
    message="'penalty' was deprecated",
    category=FutureWarning,
    module="sklearn.linear_model._logistic",
)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def transformed(values, endpoint):
    values = np.asarray(values, dtype=float)
    return np.log10(values + 1.0) if endpoint == "KSOL" else values


def raw_state(value):
    text = str(value).strip()
    if not text or text.lower() in {"nan", "na", "n/a", "none"}:
        return "blank"
    try:
        number = float(text)
        return "numeric" if np.isfinite(number) else "blank"
    except ValueError:
        return "censored" if re.match(r"^[<>≤≥]\s*[-+]?\d", text) else "other"


def mnar_probability(base_p, y, gamma):
    z = (y - y.mean()) / y.std(ddof=0)
    target = base_p.mean()
    shifted = logit(base_p) + gamma * z
    low, high = -50.0, 50.0
    for _ in range(100):
        mid = (low + high) / 2
        if expit(shifted + mid).mean() < target:
            low = mid
        else:
            high = mid
    shift = (low + high) / 2
    return expit(shifted + shift), shift


def clean(value):
    return None if pd.isna(value) else value


def scenario_key(row, include_model=False):
    seed = clean(row["seed"])
    if seed is not None:
        if float(seed) != np.floor(float(seed)):
            raise AssertionError("noninteger seed")
        seed = int(seed)
    gamma = clean(row["gamma"])
    gamma = None if gamma is None else float(gamma)
    key = (
        row["mechanism"],
        row["truth_endpoint"],
        clean(row["mask"]),
        clean(row["anchor"]),
        seed,
        gamma,
    )
    return (*key, row["model"]) if include_model else key


def expected_keys(cfg, include_model=False):
    keys = []
    models = cfg["primary_models"] if include_model else [None]
    for endpoint in cfg["truth_endpoints"]:
        for model in models:
            suffix = (model,) if include_model else ()
            for mask in cfg["mask_endpoints"]:
                keys.append(("actual", endpoint, mask, None, None, None, *suffix))
                for seed in cfg["mask_seeds"]:
                    keys.append(("mcar", endpoint, mask, None, int(seed), None, *suffix))
            for anchor in cfg["mar_anchors"]:
                for seed in cfg["mask_seeds"]:
                    keys.append(("mar", endpoint, None, anchor, int(seed), None, *suffix))
                    for gamma in cfg["mnar_gamma"]:
                        keys.append(("mnar", endpoint, None, anchor, int(seed), float(gamma), *suffix))
    return keys


def recompute_cohort(train, test, raw):
    id_col = "Molecule Name"
    if (
        train[id_col].duplicated().any()
        or test[id_col].duplicated().any()
        or raw[id_col].duplicated().any()
    ):
        raise AssertionError("duplicate molecule identifier")
    train_ids = set(train[id_col])
    test_ids = set(test[id_col])
    raw_ids = set(raw[id_col])
    if train_ids & test_ids or not (train_ids | test_ids).issubset(raw_ids):
        raise AssertionError("source identifier reconciliation mismatch")
    raw_by_id = raw.set_index(id_col)
    records = []
    numeric_match = True
    for phase, frame in [("train", train), ("test", test)]:
        raw_smiles = frame[id_col].map(raw_by_id["SMILES"]).astype(str).to_numpy()
        if not np.array_equal(frame["SMILES"].astype(str).to_numpy(), raw_smiles):
            raise AssertionError(f"raw/clean SMILES mismatch: {phase}")
        for endpoint, column in ENDPOINT_COLUMNS.items():
            clean_values = pd.to_numeric(frame[column], errors="coerce").to_numpy(float)
            raw_values = pd.to_numeric(
                frame[id_col].map(raw_by_id[column]), errors="coerce"
            ).to_numpy(float)
            clean_numeric = np.isfinite(clean_values)
            states = frame[id_col].map(raw_by_id[column]).map(raw_state)
            numeric_match &= bool(
                np.array_equal(clean_numeric, states.eq("numeric").to_numpy())
            )
            if not np.array_equal(clean_numeric, np.isfinite(raw_values)) or not np.array_equal(
                clean_values[clean_numeric], raw_values[clean_numeric]
            ):
                raise AssertionError(f"raw/clean numeric mismatch: {phase} {endpoint}")
            counts = states.value_counts().to_dict()
            records.append(
                {
                    "phase": phase,
                    "endpoint": endpoint,
                    "rows": len(frame),
                    "numeric": int(counts.get("numeric", 0)),
                    "censored": int(counts.get("censored", 0)),
                    "blank": int(counts.get("blank", 0)),
                    "other": int(counts.get("other", 0)),
                    "clean_numeric": int(clean_numeric.sum()),
                    "coverage": float(clean_numeric.mean()),
                }
            )
    summary = {
        "train_rows": len(train),
        "test_rows": len(test),
        "ml_ready_rows": len(train) + len(test),
        "raw_rows": len(raw),
        "raw_only_rows": len(raw_ids - train_ids - test_ids),
        "identifier_overlap": len(train_ids & test_ids),
        "raw_numeric_equals_clean_numeric": bool(numeric_match),
    }
    return pd.DataFrame(records), summary


def rebuild_observation_contract(train, test, cfg):
    smiles = pd.concat([train["SMILES"], test["SMILES"]], ignore_index=True)
    descriptor_registry = sorted(Descriptors._descList, key=lambda item: item[0])
    descriptor_names = [name for name, _ in descriptor_registry]
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=int(cfg["fingerprint"]["radius"]),
        fpSize=int(cfg["fingerprint"]["n_bits"]),
        includeChirality=bool(cfg["fingerprint"]["use_chirality"]),
    )
    fps = np.empty(
        (len(smiles), int(cfg["fingerprint"]["n_bits"])), dtype=np.float32
    )
    desc = np.empty((len(smiles), len(descriptor_registry)), dtype=np.float64)
    mols = []
    bitvectors = []
    canonical_iso = []
    canonical_noiso = []
    salts = []
    stereo = []
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            raise AssertionError(f"invalid source SMILES at row {i}")
        fp = generator.GetFingerprint(mol)
        arr = np.zeros(int(cfg["fingerprint"]["n_bits"]), dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fps[i] = arr
        desc[i] = [func(mol) for _, func in descriptor_registry]
        mols.append(mol)
        bitvectors.append(fp)
        canonical_iso.append(
            Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        )
        canonical_noiso.append(
            Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
        )
        salts.append("." in str(smi))
        stereo.append(bool(Chem.FindMolChiralCenters(mol, includeUnassigned=True)))
    if not np.isfinite(desc).all():
        raise AssertionError("nonfinite rebuilt descriptor")

    split = len(train)
    train_iso = canonical_iso[:split]
    test_iso = canonical_iso[split:]
    train_noiso = canonical_noiso[:split]
    test_noiso = canonical_noiso[split:]
    train_bitvectors = bitvectors[:split]
    structure = {
        "valid_train": len(train_iso),
        "valid_test": len(test_iso),
        "train_exact_duplicates": len(train_iso) - len(set(train_iso)),
        "test_exact_duplicates": len(test_iso) - len(set(test_iso)),
        "train_test_isomeric_overlap": len(set(train_iso) & set(test_iso)),
        "train_test_nonisomeric_overlap": len(set(train_noiso) & set(test_noiso)),
        "train_salts": int(np.sum(salts[:split])),
        "test_salts": int(np.sum(salts[split:])),
        "train_stereo": int(np.sum(stereo[:split])),
        "test_stereo": int(np.sum(stereo[split:])),
    }
    similarity = np.asarray(
        [
            max(DataStructs.BulkTanimotoSimilarity(fp, train_bitvectors))
            for fp in bitvectors[split:]
        ],
        dtype=np.float64,
    )
    identifier_order = (
        test["Molecule Name"].astype(str).str.extract(r"(\d+)$", expand=False)
    )
    if identifier_order.isna().any():
        raise AssertionError("unparseable molecule identifier order")
    identifier_order = identifier_order.astype(float).to_numpy()
    identifier_span = identifier_order.max() - identifier_order.min()
    if not np.isfinite(identifier_span) or identifier_span == 0:
        raise AssertionError("invalid molecule identifier order span")
    identifier_order = (identifier_order - identifier_order.min()) / identifier_span
    test_groups = []
    for i, mol in enumerate(mols[split:]):
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(
            mol=mol, includeChirality=True
        )
        test_groups.append(scaffold or f"ACYCLIC_SINGLETON_{i}")
    test_groups = np.asarray(test_groups, dtype=object)
    x_obs = np.hstack(
        [
            fps[split:],
            desc[split:],
            similarity[:, None],
            identifier_order[:, None],
        ]
    ).astype(np.float32)
    expected_shape = (
        len(test),
        int(cfg["fingerprint"]["n_bits"])
        + int(cfg["rdkit_descriptors"]["expected_count"])
        + 2,
    )
    if (
        len(descriptor_names) != int(cfg["rdkit_descriptors"]["expected_count"])
        or len(descriptor_names) != len(set(descriptor_names))
        or x_obs.shape != expected_shape
        or not np.isfinite(x_obs).all()
    ):
        raise AssertionError("rebuilt observation feature contract mismatch")
    return {
        "structure": structure,
        "descriptor_names": descriptor_names,
        "x_obs": x_obs,
        "groups": test_groups,
        "canonical_sha256": np.asarray(
            [hashlib.sha256(value.encode("utf-8")).hexdigest() for value in test_iso]
        ),
        "scaffold_sha256": np.asarray(
            [
                hashlib.sha256(str(value).encode("utf-8")).hexdigest()
                for value in test_groups
            ]
        ),
        "similarity": similarity,
        "identifier_order": identifier_order,
    }


def assert_close(left, right, label, tol=1e-12):
    if not np.isclose(float(left), float(right), rtol=0, atol=tol):
        raise AssertionError(f"{label}: {left} != {right}")


def regenerate_mask(key, endpoint_indices, y, row_audit, cfg):
    mechanism, _, mask, anchor, seed, gamma = key
    if mechanism == "actual":
        r = row_audit[f"available_{mask}"].to_numpy(int)[endpoint_indices]
        return r, None, None
    if mechanism == "mcar":
        observed = row_audit[f"available_{mask}"].to_numpy(int)[endpoint_indices]
        probability = float(observed.mean())
        r = (np.random.default_rng(seed).random(len(y)) < probability).astype(int)
        return r, probability, None
    base_p = row_audit[f"propensity_clipped_{anchor}"].to_numpy(float)[endpoint_indices]
    if mechanism == "mar":
        r = (np.random.default_rng(seed).random(len(y)) < base_p).astype(int)
        return r, float(base_p.mean()), None
    probability, shift = mnar_probability(base_p, y, gamma)
    r = (np.random.default_rng(seed).random(len(y)) < probability).astype(int)
    return r, float(probability.mean()), float(shift)


def verify(output):
    output = output.resolve()
    root = Path(__file__).resolve().parents[1]
    cfg = json.loads((root / "config" / "pilot.json").read_text(encoding="utf-8"))
    manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    decision = json.loads((output / "gate_decision.json").read_text(encoding="utf-8"))
    for name, expected in cfg["dataset"]["files"].items():
        observed = sha256(root / "data" / "raw" / name)
        if observed != expected:
            raise AssertionError(f"configured data hash mismatch: {name}")
    if (
        sha256(root / "src" / "run_pilot.py").upper() != TRUSTED_RUNNER_SHA256
        or sha256(root / "src" / "pilot_estimators.py").upper()
        != TRUSTED_ESTIMATORS_SHA256
    ):
        raise AssertionError("analysis source differs from audited trust anchor")
    required_outputs = {
        "cohort_audit.csv",
        "cohort_summary.json",
        "structure_audit.json",
        "model_metrics.csv",
        "actual_mask_diagnostics.csv",
        "scenario_results.csv",
        "pairwise_results.csv",
        "gate_decision.json",
        "summary.md",
        "predictions/row_audit.csv",
        "predictions/predictions.parquet",
    }
    required_sources = {
        "protocol.md",
        "config/pilot.json",
        "src/run_pilot.py",
        "src/pilot_estimators.py",
        "baseline/metrics.json",
        "data/raw/expansion_data_train.csv",
        "data/raw/expansion_data_test.csv",
        "data/raw/expansion_data_raw.csv",
    }
    if not required_outputs.issubset(manifest["output_hashes_before_manifest"]):
        raise AssertionError("manifest omits a required output")
    if not required_sources.issubset(manifest["hashes"]):
        raise AssertionError("manifest omits a required source")
    if (
        manifest["hashes"]["src/run_pilot.py"].upper() != TRUSTED_RUNNER_SHA256
        or manifest["hashes"]["src/pilot_estimators.py"].upper()
        != TRUSTED_ESTIMATORS_SHA256
    ):
        raise AssertionError("manifest source differs from audited trust anchor")
    for relative, expected in manifest["output_hashes_before_manifest"].items():
        if sha256(output / relative) != expected:
            raise AssertionError(f"output hash mismatch: {relative}")
    for relative, expected in manifest["hashes"].items():
        if sha256(root / relative) != expected:
            raise AssertionError(f"source hash mismatch: {relative}")
    if sha256(root / "protocol.md") != cfg["protocol_sha256"]:
        raise AssertionError("current protocol hash mismatch")
    if (
        manifest["protocol_id"] != cfg["protocol_id"]
        or manifest["mode"] != "full"
        or manifest["current_protocol_sha256"] != cfg["protocol_sha256"]
        or manifest["dataset_revision"] != cfg["dataset"]["revision"]
    ):
        raise AssertionError("manifest identity mismatch")
    if rdkit.__version__ != cfg["rdkit_descriptors"]["rdkit_version"]:
        raise AssertionError("RDKit version mismatch")

    risk = pd.read_csv(output / "scenario_results.csv")
    pair = pd.read_csv(output / "pairwise_results.csv")
    models = pd.read_csv(output / "model_metrics.csv")
    diagnostics = pd.read_csv(output / "actual_mask_diagnostics.csv")
    cohort = pd.read_csv(output / "cohort_audit.csv")
    row_audit = pd.read_csv(output / "predictions" / "row_audit.csv")
    predictions = pd.read_parquet(output / "predictions" / "predictions.parquet")
    cohort_summary = json.loads(
        (output / "cohort_summary.json").read_text(encoding="utf-8")
    )
    stored_structure = json.loads(
        (output / "structure_audit.json").read_text(encoding="utf-8")
    )
    train = pd.read_csv(root / "data" / "raw" / "expansion_data_train.csv")
    test = pd.read_csv(root / "data" / "raw" / "expansion_data_test.csv")
    raw = pd.read_csv(
        root / "data" / "raw" / "expansion_data_raw.csv",
        dtype=str,
        keep_default_na=False,
    )

    recomputed_cohort, recomputed_summary = recompute_cohort(train, test, raw)
    cohort_keys = ["phase", "endpoint"]
    expected_cohort_keys = set(
        recomputed_cohort[cohort_keys].itertuples(index=False, name=None)
    )
    observed_cohort_keys = set(cohort[cohort_keys].itertuples(index=False, name=None))
    if (
        len(cohort) != 18
        or len(recomputed_cohort) != 18
        or observed_cohort_keys != expected_cohort_keys
        or cohort.duplicated(cohort_keys).any()
    ):
        raise AssertionError("cohort grid mismatch")
    observed_cohort = cohort.set_index(cohort_keys).sort_index()
    expected_cohort = recomputed_cohort.set_index(cohort_keys).sort_index()
    for column in ["rows", "numeric", "censored", "blank", "other", "clean_numeric"]:
        if not np.array_equal(
            observed_cohort[column].to_numpy(int),
            expected_cohort[column].to_numpy(int),
        ):
            raise AssertionError(f"cohort reconstruction mismatch: {column}")
    if not np.allclose(
        observed_cohort["coverage"], expected_cohort["coverage"], rtol=0, atol=1e-12
    ):
        raise AssertionError("cohort reconstruction mismatch: coverage")
    if cohort_summary != recomputed_summary:
        raise AssertionError("cohort summary mismatch")

    observation_contract = rebuild_observation_contract(train, test, cfg)
    recomputed_structure = observation_contract["structure"]
    if stored_structure != recomputed_structure:
        raise AssertionError("structure audit mismatch")

    row_required = {
        "row_index",
        "molecule_id",
        "canonical_sha256",
        "scaffold_sha256",
        "max_train_similarity",
        "normalized_identifier_order",
    }
    for mask in cfg["mask_endpoints"]:
        row_required.update(
            {
                f"available_{mask}",
                f"propensity_raw_{mask}",
                f"propensity_clipped_{mask}",
            }
        )
    if not row_required.issubset(row_audit.columns) or len(row_audit) != len(test):
        raise AssertionError("row audit shape mismatch")
    row_indices = pd.to_numeric(row_audit["row_index"], errors="coerce").to_numpy(float)
    if (
        not np.isfinite(row_indices).all()
        or not np.equal(row_indices, np.floor(row_indices)).all()
        or not np.array_equal(row_indices.astype(int), np.arange(len(test)))
        or not np.array_equal(
            row_audit["molecule_id"].astype(str).to_numpy(),
            test["Molecule Name"].astype(str).to_numpy(),
        )
    ):
        raise AssertionError("row audit/test alignment mismatch")
    if not np.array_equal(
        row_audit["canonical_sha256"].astype(str).to_numpy(),
        observation_contract["canonical_sha256"],
    ) or not np.array_equal(
        row_audit["scaffold_sha256"].astype(str).to_numpy(),
        observation_contract["scaffold_sha256"],
    ):
        raise AssertionError("row audit structure hash mismatch")
    if not np.allclose(
        row_audit["max_train_similarity"].to_numpy(float),
        observation_contract["similarity"],
        rtol=0,
        atol=1e-12,
    ) or not np.allclose(
        row_audit["normalized_identifier_order"].to_numpy(float),
        observation_contract["identifier_order"],
        rtol=0,
        atol=1e-12,
    ):
        raise AssertionError("row audit observation covariate mismatch")
    lo, hi = map(float, cfg["propensity_clip"])
    x_obs = observation_contract["x_obs"]
    groups = observation_contract["groups"]
    folds = list(
        GroupKFold(n_splits=int(cfg["availability_auc_folds"])).split(
            x_obs, groups=groups
        )
    )
    prepared = []
    for train_index, valid_index in folds:
        if set(groups[train_index]) & set(groups[valid_index]):
            raise AssertionError("scaffold leakage in rebuilt propensity folds")
        scaler = StandardScaler()
        x_train = scaler.fit_transform(x_obs[train_index])
        x_valid = scaler.transform(x_obs[valid_index])
        prepared.append((train_index, valid_index, x_train, x_valid))
    propensity_refit_max_difference = 0.0
    for mask in cfg["mask_endpoints"]:
        expected_available = np.isfinite(
            pd.to_numeric(test[ENDPOINT_COLUMNS[mask]], errors="coerce").to_numpy(float)
        ).astype(int)
        observed_available_raw = pd.to_numeric(
            row_audit[f"available_{mask}"], errors="coerce"
        ).to_numpy(float)
        if (
            not np.isfinite(observed_available_raw).all()
            or not np.isin(observed_available_raw, [0.0, 1.0]).all()
        ):
            raise AssertionError(f"invalid availability indicator: {mask}")
        observed_available = observed_available_raw.astype(int)
        p_raw = row_audit[f"propensity_raw_{mask}"].to_numpy(float)
        p = row_audit[f"propensity_clipped_{mask}"].to_numpy(float)
        if not np.array_equal(observed_available, expected_available):
            raise AssertionError(f"availability mismatch: {mask}")
        if not np.isfinite(p_raw).all() or not np.isfinite(p).all():
            raise AssertionError(f"nonfinite propensity: {mask}")
        if not np.allclose(p, np.clip(p_raw, lo, hi), rtol=0, atol=1e-12):
            raise AssertionError(f"propensity clipping mismatch: {mask}")
        refit_raw = np.empty(len(test), dtype=float)
        for train_index, valid_index, x_train, x_valid in prepared:
            if len(np.unique(expected_available[train_index])) != 2:
                raise AssertionError(f"single-class propensity fold: {mask}")
            model = LogisticRegression(
                penalty=cfg["propensity_model"]["penalty"],
                C=float(cfg["propensity_model"]["C"]),
                solver=cfg["propensity_model"]["solver"],
                class_weight=cfg["propensity_model"]["class_weight"],
                max_iter=int(cfg["propensity_model"]["max_iter"]),
                random_state=int(cfg["propensity_model"]["random_state"]),
            ).fit(x_train, expected_available[train_index])
            refit_raw[valid_index] = model.predict_proba(x_valid)[:, 1]
        difference = float(np.max(np.abs(refit_raw - p_raw)))
        propensity_refit_max_difference = max(
            propensity_refit_max_difference, difference
        )
        if difference > 1e-12 or not np.allclose(
            p, np.clip(refit_raw, lo, hi), rtol=0, atol=1e-12
        ):
            raise AssertionError(f"propensity refit mismatch: {mask}")

    prediction_columns = {
        "row_index",
        "molecule_id",
        "truth_endpoint",
        "model",
        "prediction",
    }
    if not prediction_columns.issubset(predictions.columns):
        raise AssertionError("prediction column mismatch")
    prediction_indices = pd.to_numeric(
        predictions["row_index"], errors="coerce"
    ).to_numpy(float)
    if not np.isfinite(prediction_indices).all() or not np.equal(
        prediction_indices, np.floor(prediction_indices)
    ).all():
        raise AssertionError("invalid prediction row index")
    predictions["row_index"] = prediction_indices.astype(int)
    expected_prediction_keys = {
        (endpoint, model, i)
        for endpoint in cfg["truth_endpoints"]
        for model in ALL_MODELS
        for i in range(len(test))
    }
    observed_prediction_keys = set(
        predictions[["truth_endpoint", "model", "row_index"]].itertuples(
            index=False, name=None
        )
    )
    if (
        len(predictions) != len(expected_prediction_keys)
        or predictions.duplicated(["truth_endpoint", "model", "row_index"]).any()
        or observed_prediction_keys != expected_prediction_keys
        or not np.isfinite(predictions["prediction"].to_numpy(float)).all()
    ):
        raise AssertionError("prediction grid mismatch")
    expected_prediction_ids = test["Molecule Name"].astype(str).to_numpy()[
        predictions["row_index"].to_numpy(int)
    ]
    if not np.array_equal(
        predictions["molecule_id"].astype(str).to_numpy(), expected_prediction_ids
    ):
        raise AssertionError("prediction/test identifier mismatch")

    expected_model_keys = {
        (endpoint, model)
        for endpoint in cfg["truth_endpoints"]
        for model in ALL_MODELS
    }
    if (
        len(models) != len(expected_model_keys)
        or models.duplicated(["truth_endpoint", "model"]).any()
        or set(
            models[["truth_endpoint", "model"]].itertuples(index=False, name=None)
        )
        != expected_model_keys
    ):
        raise AssertionError("model metric grid mismatch")

    observed_risk_keys = [scenario_key(row, True) for _, row in risk.iterrows()]
    observed_pair_keys = [scenario_key(row, False) for _, row in pair.iterrows()]
    if len(observed_risk_keys) != 2260 or len(set(observed_risk_keys)) != 2260 or set(observed_risk_keys) != set(expected_keys(cfg, True)):
        raise AssertionError("risk grid mismatch")
    if len(observed_pair_keys) != 1130 or len(set(observed_pair_keys)) != 1130 or set(observed_pair_keys) != set(expected_keys(cfg, False)):
        raise AssertionError("pair grid mismatch")
    if not (
        (pair["model_a"] == cfg["primary_models"][0]).all()
        and (pair["model_b"] == cfg["primary_models"][1]).all()
    ):
        raise AssertionError("paired model identity mismatch")

    prediction_lookup = {
        (endpoint, model): part.sort_values("row_index")["prediction"].to_numpy(float)
        for (endpoint, model), part in predictions.groupby(["truth_endpoint", "model"])
    }
    model_lookup = models.set_index(["truth_endpoint", "model"])
    recomputed_metrics = {}
    endpoint_coverage = {}
    for endpoint in cfg["truth_endpoints"]:
        train_y = pd.to_numeric(
            train[ENDPOINT_COLUMNS[endpoint]], errors="coerce"
        ).to_numpy(float)
        test_y = pd.to_numeric(
            test[ENDPOINT_COLUMNS[endpoint]], errors="coerce"
        ).to_numpy(float)
        observed = np.isfinite(test_y)
        y = transformed(test_y[observed], endpoint)
        endpoint_coverage[endpoint] = float(observed.mean())
        endpoint_mae = {
            model: float(
                np.abs(y - prediction_lookup[(endpoint, model)][observed]).mean()
            )
            for model in ALL_MODELS
        }
        denominator = endpoint_mae["median"]
        for model in ALL_MODELS:
            expected = {
                "train_n": int(np.isfinite(train_y).sum()),
                "test_n": int(observed.sum()),
                "test_coverage": endpoint_coverage[endpoint],
                "full_truth": endpoint_mae[model],
                "normalized_mae": endpoint_mae[model] / denominator,
            }
            recomputed_metrics[(endpoint, model)] = expected
            stored = model_lookup.loc[(endpoint, model)]
            if int(stored["train_n"]) != expected["train_n"] or int(
                stored["test_n"]
            ) != expected["test_n"]:
                raise AssertionError(f"model sample count mismatch: {endpoint} {model}")
            for column in ["test_coverage", "full_truth", "normalized_mae"]:
                assert_close(
                    stored[column],
                    expected[column],
                    f"model metric {endpoint} {model} {column}",
                )

    risk_lookup = {scenario_key(row, True): row for _, row in risk.iterrows()}
    max_row_difference = 0.0
    for endpoint in cfg["truth_endpoints"]:
        raw_y = pd.to_numeric(test[ENDPOINT_COLUMNS[endpoint]], errors="coerce").to_numpy(float)
        indices = np.flatnonzero(np.isfinite(raw_y))
        y = transformed(raw_y[indices], endpoint)
        losses = {
            model: np.abs(y - prediction_lookup[(endpoint, model)][indices])
            for model in cfg["primary_models"]
        }
        endpoint_keys = [key for key in expected_keys(cfg, False) if key[1] == endpoint]
        for key in endpoint_keys:
            r, generator_mean, shift = regenerate_mask(key, indices, y, row_audit, cfg)
            for model in cfg["primary_models"]:
                row = risk_lookup[(*key, model)]
                expected_values = {
                    "n_total": len(y),
                    "n_observed": int(r.sum()),
                    "observed_fraction": float(r.mean()),
                    "full_truth": float(losses[model].mean()),
                    "cc": float(losses[model][r.astype(bool)].mean()),
                }
                for column, expected in expected_values.items():
                    difference = abs(float(row[column]) - float(expected))
                    max_row_difference = max(max_row_difference, difference)
                    if difference > 1e-12:
                        raise AssertionError(f"scenario replay mismatch {(*key, model)} {column}")
                if generator_mean is None:
                    if not pd.isna(row["generator_mean_probability"]):
                        raise AssertionError("actual row has generator probability")
                else:
                    assert_close(row["generator_mean_probability"], generator_mean, "generator mean")
                if shift is None:
                    if not pd.isna(row["mnar_intercept_shift"]):
                        raise AssertionError("non-MNAR row has shift")
                else:
                    assert_close(row["mnar_intercept_shift"], shift, "MNAR shift")

    pair_lookup = {scenario_key(row, False): row for _, row in pair.iterrows()}
    max_pair_difference = 0.0
    for key, pair_row in pair_lookup.items():
        row_a = risk_lookup[(*key, cfg["primary_models"][0])]
        row_b = risk_lookup[(*key, cfg["primary_models"][1])]
        for pair_column, risk_column in [
            ("full_truth_delta", "full_truth"),
            ("cc_delta", "cc"),
            ("ipw_delta", "ipw"),
            ("aipw_delta", "aipw"),
        ]:
            expected = float(row_a[risk_column] - row_b[risk_column])
            difference = abs(float(pair_row[pair_column]) - expected)
            max_pair_difference = max(max_pair_difference, difference)
            if difference > 1e-12:
                raise AssertionError(f"paired estimate mismatch {key} {pair_column}")
        assert_close(pair_row["aipw_delta_lcl"], pair_row["aipw_delta"] - 1.96 * pair_row["aipw_delta_se"], "pair lcl")
        assert_close(pair_row["aipw_delta_ucl"], pair_row["aipw_delta"] + 1.96 * pair_row["aipw_delta_se"], "pair ucl")

    for _, row in risk.iterrows():
        assert_close(row["lcl"], row["aipw"] - 1.96 * row["se"], "risk lcl")
        assert_close(row["ucl"], row["aipw"] + 1.96 * row["se"], "risk ucl")
        assert_close(row["cc_abs_bias"], abs(row["cc"] - row["full_truth"]), "CC bias")
        assert_close(row["aipw_abs_bias"], abs(row["aipw"] - row["full_truth"]), "AIPW bias")
        assert_close(row["ci_width"], row["ucl"] - row["lcl"], "CI width")

    recomputed_diagnostics = []
    for mask in cfg["mask_endpoints"]:
        r = row_audit[f"available_{mask}"].to_numpy(int)
        p_raw = row_audit[f"propensity_raw_{mask}"].to_numpy(float)
        p = row_audit[f"propensity_clipped_{mask}"].to_numpy(float)
        weights = r / p
        recomputed_diagnostics.append(
            {
                "mask": mask,
                "coverage": float(r.mean()),
                "observed_n": int(r.sum()),
                "auc": roc_auc_score(r, p_raw),
                "clip_fraction": np.mean((p_raw < cfg["propensity_clip"][0]) | (p_raw > cfg["propensity_clip"][1])),
                "ess": weights.sum() ** 2 / np.square(weights).sum(),
                "ess_fraction": (weights.sum() ** 2 / np.square(weights).sum()) / r.sum(),
            }
        )
    if (
        len(diagnostics) != len(cfg["mask_endpoints"])
        or diagnostics["mask"].duplicated().any()
        or set(diagnostics["mask"]) != set(cfg["mask_endpoints"])
    ):
        raise AssertionError("diagnostic grid mismatch")
    diag_lookup = diagnostics.set_index("mask")
    for item in recomputed_diagnostics:
        if int(diag_lookup.loc[item["mask"], "observed_n"]) != item["observed_n"]:
            raise AssertionError(f"diagnostic observed_n mismatch: {item['mask']}")
        for column in ["coverage", "auc", "clip_fraction", "ess", "ess_fraction"]:
            assert_close(diag_lookup.loc[item["mask"], column], item[column], f"diagnostic {item['mask']} {column}")
    gate4_qualifying = sum(
        item["auc"] >= cfg["gates"]["minimum_propensity_auc"]
        and item["ess_fraction"] >= cfg["gates"]["minimum_ess_fraction"]
        for item in recomputed_diagnostics
    )
    gate4_expected_pass = gate4_qualifying >= cfg["gates"]["minimum_predictable_masks"]
    stored_gate4 = decision["gates"]["gate4"]
    if (
        stored_gate4["pass"] != gate4_expected_pass
        or int(stored_gate4["qualifying_masks"]) != gate4_qualifying
        or int(stored_gate4["required_masks"])
        != int(cfg["gates"]["minimum_predictable_masks"])
        or bool(stored_gate4.get("integrity_failure", True))
    ):
        raise AssertionError("gate4 pass mismatch")
    if len(stored_gate4["details"]) != len(cfg["mask_endpoints"]):
        raise AssertionError("gate4 detail row count mismatch")
    stored_gate4_details = {
        item["mask"]: item for item in stored_gate4["details"]
    }
    if len(stored_gate4_details) != len(stored_gate4["details"]) or set(
        stored_gate4_details
    ) != set(cfg["mask_endpoints"]):
        raise AssertionError("gate4 detail grid mismatch")
    for item in recomputed_diagnostics:
        stored = stored_gate4_details[item["mask"]]
        qualifies = bool(
            item["auc"] >= cfg["gates"]["minimum_propensity_auc"]
            and item["ess_fraction"] >= cfg["gates"]["minimum_ess_fraction"]
        )
        if bool(stored["qualifies"]) != qualifies:
            raise AssertionError(f"gate4 qualifier mismatch: {item['mask']}")
        if int(stored["observed_n"]) != item["observed_n"]:
            raise AssertionError(f"gate4 observed_n mismatch: {item['mask']}")
        for column in ["coverage", "auc", "clip_fraction", "ess", "ess_fraction"]:
            assert_close(
                stored[column], item[column], f"gate4 detail {item['mask']} {column}"
            )

    full_primary = pd.Series(
        {
            (endpoint, model): recomputed_metrics[(endpoint, model)]["full_truth"]
            for endpoint in cfg["truth_endpoints"]
            for model in cfg["primary_models"]
        }
    )
    gate5_actual = []
    gate5_mar = []
    for endpoint in cfg["truth_endpoints"]:
        truth = {model: full_primary.loc[(endpoint, model)] for model in cfg["primary_models"]}
        truth_values = list(truth.values())
        truth_winner = None if truth_values[0] == truth_values[1] else cfg["primary_models"][int(truth_values[1] < truth_values[0])]
        for mask in cfg["mask_endpoints"]:
            part = risk[(risk["mechanism"] == "actual") & (risk["truth_endpoint"] == endpoint) & (risk["mask"] == mask)].set_index("model")
            relative_bias = max(abs(part.loc[model, "cc"] - truth[model]) / truth[model] for model in cfg["primary_models"])
            cc_values = [part.loc[model, "cc"] for model in cfg["primary_models"]]
            cc_winner = None if cc_values[0] == cc_values[1] else cfg["primary_models"][int(cc_values[1] < cc_values[0])]
            rank_change = bool(
                truth_winner is not None
                and cc_winner is not None
                and cc_winner != truth_winner
            )
            gate5_actual.append(
                {
                    "truth_endpoint": endpoint,
                    "mask": mask,
                    "max_relative_bias": float(relative_bias),
                    "rank_change": rank_change,
                    "triggers": bool(
                        relative_bias >= cfg["gates"]["minimum_cc_relative_bias"]
                        or rank_change
                    ),
                }
            )
        for anchor in cfg["mar_anchors"]:
            part = risk[(risk["mechanism"] == "mar") & (risk["truth_endpoint"] == endpoint) & (risk["anchor"] == anchor)]
            median_bias = max(
                np.median(abs(model_part["cc"] - truth[model]) / truth[model])
                for model, model_part in part.groupby("model")
            )
            changes = 0
            if truth_winner is not None:
                for _, seed_part in part.groupby("seed"):
                    values = seed_part.set_index("model")["cc"]
                    cc_pair = [values.loc[model] for model in cfg["primary_models"]]
                    cc_winner = None if cc_pair[0] == cc_pair[1] else cfg["primary_models"][int(cc_pair[1] < cc_pair[0])]
                    changes += int(cc_winner is not None and cc_winner != truth_winner)
            gate5_mar.append(
                {
                    "truth_endpoint": endpoint,
                    "anchor": anchor,
                    "max_model_median_relative_bias": float(median_bias),
                    "winner_changes": int(changes),
                    "triggers": bool(
                        median_bias >= cfg["gates"]["minimum_cc_relative_bias"]
                        or changes >= 21
                    ),
                }
            )
    gate5_expected_pass = any(
        item["triggers"] for item in gate5_actual + gate5_mar
    )
    stored_gate5 = decision["gates"]["gate5"]
    if (
        stored_gate5["pass"] != gate5_expected_pass
        or bool(stored_gate5.get("integrity_failure", True))
    ):
        raise AssertionError("gate5 pass mismatch")
    if len(stored_gate5["actual"]) != len(gate5_actual) or len(
        stored_gate5["mar"]
    ) != len(gate5_mar):
        raise AssertionError("gate5 detail row count mismatch")
    stored_actual = {
        (item["truth_endpoint"], item["mask"]): item
        for item in stored_gate5["actual"]
    }
    stored_mar = {
        (item["truth_endpoint"], item["anchor"]): item
        for item in stored_gate5["mar"]
    }
    if len(stored_actual) != len(stored_gate5["actual"]) or len(
        stored_mar
    ) != len(stored_gate5["mar"]):
        raise AssertionError("duplicate gate5 detail key")
    if set(stored_actual) != {
        (item["truth_endpoint"], item["mask"]) for item in gate5_actual
    } or set(stored_mar) != {
        (item["truth_endpoint"], item["anchor"]) for item in gate5_mar
    }:
        raise AssertionError("gate5 detail grid mismatch")
    for item in gate5_actual:
        stored = stored_actual[(item["truth_endpoint"], item["mask"])]
        assert_close(
            stored["max_relative_bias"], item["max_relative_bias"], "gate5 actual bias"
        )
        if bool(stored["rank_change"]) != item["rank_change"] or bool(
            stored["triggers"]
        ) != item["triggers"]:
            raise AssertionError("gate5 actual detail mismatch")
    for item in gate5_mar:
        stored = stored_mar[(item["truth_endpoint"], item["anchor"])]
        assert_close(
            stored["max_model_median_relative_bias"],
            item["max_model_median_relative_bias"],
            "gate5 MAR bias",
        )
        if int(stored["winner_changes"]) != item["winner_changes"] or bool(
            stored["triggers"]
        ) != item["triggers"]:
            raise AssertionError("gate5 MAR detail mismatch")

    mar = risk[risk["mechanism"] == "mar"]
    cc_bias = np.abs(mar["cc"] - mar["full_truth"])
    aipw_bias = np.abs(mar["aipw"] - mar["full_truth"])
    gate6_reduction = float(1 - np.median(aipw_bias) / np.median(cc_bias))
    gate6_count = int(((mar["lcl"] <= mar["full_truth"]) & (mar["ucl"] >= mar["full_truth"])).sum())
    gate6_coverage = float(gate6_count / 480)
    gate6_width = mar["ucl"] - mar["lcl"]
    gate6_ess_ok = bool(
        (mar["ess_fraction"] >= cfg["gates"]["minimum_ess_fraction"]).all()
    )
    stored_gate6 = decision["gates"]["gate6"]
    assert_close(stored_gate6["bias_reduction"], gate6_reduction, "gate6 reduction")
    assert_close(stored_gate6["coverage"], gate6_coverage, "gate6 coverage")
    assert_close(
        stored_gate6["minimum_ess_fraction"],
        mar["ess_fraction"].min(),
        "gate6 minimum ESS",
    )
    assert_close(
        stored_gate6["interval_width_median"], gate6_width.median(), "gate6 width median"
    )
    assert_close(
        stored_gate6["interval_width_q05"],
        gate6_width.quantile(0.05),
        "gate6 width q05",
    )
    assert_close(
        stored_gate6["interval_width_q95"],
        gate6_width.quantile(0.95),
        "gate6 width q95",
    )
    if (
        int(stored_gate6["row_count"]) != 480
        or int(stored_gate6["coverage_count"]) != gate6_count
        or bool(stored_gate6["all_ess_acceptable"]) != gate6_ess_ok
        or bool(stored_gate6.get("integrity_failure", True))
    ):
        raise AssertionError("gate6 coverage count mismatch")
    gate6_expected_pass = bool(
        gate6_reduction >= cfg["gates"]["minimum_aipw_bias_reduction"]
        and gate6_coverage >= cfg["gates"]["minimum_aipw_ci_coverage"]
        and gate6_ess_ok
    )
    if stored_gate6["pass"] != gate6_expected_pass:
        raise AssertionError("gate6 pass mismatch")

    expected_counts = {
        "train_rows": 5326,
        "test_rows": 2282,
        "ml_ready_rows": 7608,
        "raw_rows": 7618,
        "raw_only_rows": 10,
        "identifier_overlap": 0,
    }
    state_partition_ok = bool(
        (
            recomputed_cohort["numeric"]
            + recomputed_cohort["censored"]
            + recomputed_cohort["blank"]
            == recomputed_cohort["rows"]
        ).all()
        and (recomputed_cohort["other"] == 0).all()
        and (
            recomputed_cohort["numeric"] == recomputed_cohort["clean_numeric"]
        ).all()
    )
    gate1_expected_pass = bool(
        all(recomputed_summary[key] == value for key, value in expected_counts.items())
        and recomputed_summary["ml_ready_rows"]
        >= cfg["gates"]["minimum_parent_rows"]
        and recomputed_summary["raw_numeric_equals_clean_numeric"]
        and len(recomputed_cohort) == 18
        and state_partition_ok
    )
    stored_gate1 = decision["gates"]["gate1"]
    if (
        bool(stored_gate1["pass"]) != gate1_expected_pass
        or bool(stored_gate1["state_partition_ok"]) != state_partition_ok
        or not bool(stored_gate1["cohort_grid_ok"])
    ):
        raise AssertionError("gate1 pass mismatch")
    for key, value in recomputed_summary.items():
        if stored_gate1[key] != value:
            raise AssertionError(f"gate1 field mismatch: {key}")

    expected_feature_shape = [
        len(test),
        int(cfg["fingerprint"]["n_bits"])
        + int(cfg["rdkit_descriptors"]["expected_count"])
        + 2,
    ]
    expected_descriptor_names = observation_contract["descriptor_names"]
    feature_contract_ok = bool(
        manifest["features"]["observation_feature_shape"] == expected_feature_shape
        and manifest["features"]["rdkit_descriptor_count"]
        == cfg["rdkit_descriptors"]["expected_count"]
        and manifest["features"]["rdkit_descriptor_names"]
        == expected_descriptor_names
        and manifest["features"]["fingerprint"] == cfg["fingerprint"]
    )
    if not feature_contract_ok:
        raise AssertionError("feature contract mismatch")
    gate2_expected_pass = bool(
        all(
            coverage >= cfg["gates"]["minimum_truth_coverage"]
            for coverage in endpoint_coverage.values()
        )
    )
    qualifying_models = [
        model
        for model in ["ecfp_ridge", "ecfp_rf", "ecfp_desc_xgb"]
        if all(
            recomputed_metrics[(endpoint, model)]["normalized_mae"]
            < cfg["gates"]["maximum_normalized_mae"]
            for endpoint in cfg["truth_endpoints"]
        )
    ]
    gate3_expected_pass = len(qualifying_models) >= 2
    stored_gate2 = decision["gates"]["gate2"]
    stored_gate3 = decision["gates"]["gate3"]
    if stored_gate2["pass"] != gate2_expected_pass or stored_gate3["pass"] != gate3_expected_pass:
        raise AssertionError("gate2/3 pass mismatch")
    if set(stored_gate2["coverage"]) != set(endpoint_coverage):
        raise AssertionError("gate2 coverage grid mismatch")
    for endpoint, coverage in endpoint_coverage.items():
        assert_close(stored_gate2["coverage"][endpoint], coverage, "gate2 coverage")
    if (
        set(stored_gate3["qualifying_non_dummy_models"]) != set(qualifying_models)
        or int(stored_gate3["required"]) != 2
    ):
        raise AssertionError("gate3 detail mismatch")

    risk_numeric = ["cc", "ipw", "aipw", "se", "lcl", "ucl", "ess", "ess_fraction", "full_truth"]
    pair_numeric = ["full_truth_delta", "cc_delta", "ipw_delta", "aipw_delta", "aipw_delta_se", "aipw_delta_lcl", "aipw_delta_ucl"]
    finite_risk = bool(np.isfinite(risk[risk_numeric].to_numpy(float)).all())
    finite_pair = bool(np.isfinite(pair[pair_numeric].to_numpy(float)).all())
    structural_ok = bool(
        recomputed_structure["train_test_isomeric_overlap"] == 0
        and recomputed_structure["train_test_nonisomeric_overlap"] == 0
        and recomputed_structure["train_exact_duplicates"] == 0
        and recomputed_structure["test_exact_duplicates"] == 0
        and recomputed_structure["valid_train"] == len(train)
        and recomputed_structure["valid_test"] == len(test)
    )
    gate7_expected_pass = bool(
        structural_ok
        and finite_risk
        and finite_pair
        and len(risk) == 2260
        and len(pair) == 1130
        and feature_contract_ok
        and cfg["primary_models"] == ["ecfp_rf", "ecfp_desc_xgb"]
        and rdkit.__version__ == cfg["rdkit_descriptors"]["rdkit_version"]
        and sha256(root / "protocol.md") == cfg["protocol_sha256"]
    )
    stored_gate7 = decision["gates"]["gate7"]
    if bool(stored_gate7["pass"]) != gate7_expected_pass:
        raise AssertionError("gate7 pass mismatch")
    expected_gate7_fields = {
        "structural_ok": structural_ok,
        "finite_required_metrics": finite_risk,
        "finite_pairwise_metrics": finite_pair,
        "expected_risk_rows": 2260,
        "observed_risk_rows": len(risk),
        "observed_pairwise_rows": len(pair),
        "risk_grid_ok": True,
        "pair_grid_ok": True,
        "full_truth_matches_canonical_table": True,
        "pair_truth_matches_canonical_table": True,
        "gate_input_integrity": True,
        "outcome_features_used": False,
        "observation_feature_shape": expected_feature_shape,
        "expected_observation_feature_shape": expected_feature_shape,
        "primary_model_order_matches_config": True,
    }
    for key, value in expected_gate7_fields.items():
        if stored_gate7[key] != value:
            raise AssertionError(f"gate7 field mismatch: {key}")
    if stored_gate7.get("risk_grid_reason") is not None or stored_gate7.get(
        "pair_grid_reason"
    ) is not None:
        raise AssertionError("gate7 grid reason mismatch")

    recomputed_gate_passes = {
        "gate1": gate1_expected_pass,
        "gate2": gate2_expected_pass,
        "gate3": gate3_expected_pass,
        "gate4": gate4_expected_pass,
        "gate5": gate5_expected_pass,
        "gate6": gate6_expected_pass,
        "gate7": gate7_expected_pass,
    }
    critical = all(recomputed_gate_passes[name] for name in ["gate1", "gate2", "gate3", "gate7"])
    expected_decision = (
        "GO_FULL"
        if critical and all(recomputed_gate_passes[name] for name in ["gate4", "gate5", "gate6"])
        else (
            "STOP_PILOT_NEGATIVE"
            if critical
            else "STOP_CRITICAL"
        )
    )
    if decision["decision"] != expected_decision:
        raise AssertionError("overall decision mismatch")

    return {
        "status": "PASS",
        "verification_scope": "full-run artifact hashes and grids, raw/clean row alignment, structure/observation-feature reconstruction, five-fold actual-propensity refit, deterministic mask replay, model-metric replay, and gate aggregation; excludes scenario-specific three-fold nuisance refits",
        "full_nuisance_estimators_refit": False,
        "actual_propensity_refit_max_abs_difference": propensity_refit_max_difference,
        "decision": decision["decision"],
        "risk_rows": len(risk),
        "pairwise_rows": len(pair),
        "maximum_scenario_replay_difference": max_row_difference,
        "maximum_pair_identity_difference": max_pair_difference,
        "gate6_bias_reduction": gate6_reduction,
        "gate6_coverage_count": gate6_count,
        "gate6_minimum_ess_fraction": float(mar["ess_fraction"].min()),
        "gate4_qualifying_masks": int(gate4_qualifying),
        "gate5_trigger": bool(gate5_expected_pass),
        "output_hashes_verified": len(manifest["output_hashes_before_manifest"]),
        "source_hashes_verified": len(manifest["hashes"]),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.output), indent=2))


if __name__ == "__main__":
    main()
