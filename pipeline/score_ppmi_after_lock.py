#!/usr/bin/env python3
"""Score fixed PPMI controls and models only after the public second lock.

This executor is intentionally not run during lock preparation. It prints only
aggregate results; row-ordered probabilities stay in owner-only Mac storage.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import time
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
from sklearn.metrics import balanced_accuracy_score, log_loss, roc_auc_score

from freeze_ppmi_input_binding import freeze, geo_array_hash
from frozen_geo_scorer import FrozenGEOScorer, MODELS
from paired_bootstrap import paired_bootstrap


SEX_BOOTSTRAP_SEED = 3222988182
PRIMARY_BOOTSTRAP_SEED = 20261014
N_BOOTSTRAPS = 5000


def progress(done: float, total: int, started: float, label: str) -> None:
    elapsed = time.monotonic() - started
    fraction = min(1.0, done / total)
    filled = round(24 * fraction)
    eta = elapsed * (1 - fraction) / fraction if fraction else 0.0
    print(f"Progress [{'#' * filled}{'-' * (24 - filled)}] {100*fraction:5.1f}% "
          f"| elapsed {elapsed/60:.1f}m | ETA ~{eta/60:.1f}m | {label}", flush=True)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_public_zenodo_record(lock: dict, manifest_path: Path) -> None:
    doi = lock["zenodo_doi"]
    record_id = doi.rsplit(".", 1)[1]
    url = f"https://zenodo.org/api/records/{record_id}"
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            record = json.load(response)
    except Exception as exc:
        raise ValueError("Cannot verify the public, timestamped Zenodo record") from exc
    if (record.get("status") != "published" or record.get("submitted") is not True
            or record.get("doi") != doi or record.get("metadata", {}).get("doi") != doi):
        raise ValueError("Zenodo DOI is not a published version of this lock")
    matching = [entry for entry in record.get("files", [])
                if entry.get("key") == manifest_path.name]
    if len(matching) != 1:
        raise ValueError("Published Zenodo record lacks the exact lock manifest")
    content_url = f"https://zenodo.org/api/records/{record_id}/files/{manifest_path.name}/content"
    try:
        with urllib.request.urlopen(content_url, timeout=30) as response:
            published_bytes = response.read()
    except Exception as exc:
        raise ValueError("Cannot download the published Zenodo lock manifest") from exc
    if hashlib.sha256(published_bytes).digest() != hashlib.sha256(manifest_path.read_bytes()).digest():
        raise ValueError("Local lock manifest differs from the public Zenodo copy")


def verify_lock(lock: dict, manifest_path: Path, bundle: Path, final_dir: Path, marker_manifest: Path,
                sex_artifact: Path) -> None:
    if lock.get("status") != "published_second_lock":
        raise ValueError("Published second-lock status is absent")
    if not re.fullmatch(r"10\.5281/zenodo\.[0-9]+", str(lock.get("zenodo_doi", ""))):
        raise ValueError("Published Zenodo DOI is absent")
    if not re.fullmatch(r"[0-9a-f]{40}", str(lock.get("git_commit", ""))):
        raise ValueError("Second-lock Git commit is absent")
    verify_public_zenodo_record(lock, manifest_path)
    expected = lock.get("sha256", {})
    inputs = {
        "bundle": bundle,
        "geo_preprocessing": final_dir / "preprocessing.npz",
        "geo_decisions": final_dir / "decisions.json",
        "marker_manifest": marker_manifest,
        "sex_artifact": sex_artifact,
        "frozen_geo_scorer_code": Path(__file__).with_name("frozen_geo_scorer.py"),
        "ppmi_executor_code": Path(__file__),
        "population_audit_code": Path(__file__).with_name("freeze_ppmi_input_binding.py"),
        "candidate_audit_code": Path(__file__).with_name("audit_ppmi_candidate_bundle.py"),
        "paired_bootstrap_code": Path(__file__).with_name("paired_bootstrap.py"),
    }
    inputs.update({f"{model}_model": final_dir / f"{model}_model.npz" for model in MODELS})
    if set(expected) != set(inputs):
        raise ValueError("Second-lock file-hash manifest is incomplete or unexpected")
    for name, path in inputs.items():
        if digest(path) != expected[name]:
            raise ValueError(f"Second-lock SHA-256 mismatch: {name}")
    for key in ("ordered_population_sha256", "raw_matrix_sha256",
                "scored_matrix_geo_array_hash_sha256"):
        if not re.fullmatch(r"[0-9a-f]{64}", str(lock.get(key, ""))):
            raise ValueError(f"Second-lock input digest absent: {key}")
    if lock.get("final_population_counts") != {"PD": 311, "Control": 161}:
        raise ValueError("Second-lock final population counts differ")
    geo_sex = lock.get("geo_sex_transfer_auroc", {})
    if (set(geo_sex) != {"PROBE_to_GENEPARK", "GENEPARK_to_PROBE"}
            or any(not isinstance(value, (int, float)) or value < .95 or value > 1
                   for value in geo_sex.values())):
        raise ValueError("GEO sex-control development gate is not documented as passed")


def load_locked_input(bundle: Path, binding: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with zipfile.ZipFile(bundle) as archive:
        raw = archive.read("rank_features.f32")
        metadata = archive.read("candidate_metadata.jsonl")
        reference = archive.read("common_genes_gencode_v29.tsv")
    if hashlib.sha256(raw).hexdigest() != binding["raw_matrix_sha256"]:
        raise ValueError("Input matrix changed after binding")
    n, g = binding["raw_matrix_shape"]
    lattice = np.frombuffer(raw, dtype="<f4").reshape(n, g).astype(np.float64) * (2 * g)
    if np.max(np.abs(lattice - np.rint(lattice))) >= .01:
        raise ValueError("Input ranks no longer match frozen GEO lattice")
    ranks = np.rint(lattice) / (2 * g)
    if geo_array_hash(ranks) != binding["scored_matrix_geo_array_hash_sha256"]:
        raise ValueError("Recovered scored matrix changed after binding")
    rows = [json.loads(line) for line in metadata.decode("utf-8").splitlines()]
    if len(rows) != n:
        raise ValueError("Candidate metadata row count changed")
    genes = np.array([
        row["gene_id"] for row in csv.DictReader(io.StringIO(reference.decode("utf-8")), delimiter="\t")
    ])
    try:
        age = np.array([float(row["age_at_BL_visit"]) for row in rows], dtype=np.float64)
        sex = np.array([{"female": 0, "male": 1}[str(row["sex"]).strip().lower()] for row in rows], dtype=np.int64)
        y = np.array([{"Control": 0, "PD": 1}[row["historical_label"]] for row in rows], dtype=np.int64)
    except (KeyError, ValueError, TypeError) as exc:
        raise ValueError("Locked candidate metadata schema is invalid") from exc
    donors = np.array([str(row["PATNO"]) for row in rows])
    if len(set(donors.tolist())) != len(donors):
        raise ValueError("Locked candidate donor IDs are not unique")
    return ranks, age, sex, y, donors, genes


def sex_gate(y: np.ndarray, sex: np.ndarray, predictions: np.ndarray,
             tick=None) -> dict:
    control = y == 0
    labels = sex[control]
    scores = predictions[control]
    if not np.isfinite(scores).all():
        raise ValueError("Control sex predictions must be finite")
    if set(np.unique(labels)) != {0, 1}:
        return {"status": "unevaluable_one_sex_class", "pass": False, "control_count": int(control.sum())}
    point = float(roc_auc_score(labels, scores))
    rng = np.random.default_rng(SEX_BOOTSTRAP_SEED)
    strata = [np.flatnonzero(labels == value) for value in (0, 1)]
    draws = np.empty(N_BOOTSTRAPS, dtype=float)
    for i in range(N_BOOTSTRAPS):
        index = np.concatenate([rng.choice(ix, size=len(ix), replace=True) for ix in strata])
        draws[i] = roc_auc_score(labels[index], scores[index])
        if tick is not None and (i + 1) % 500 == 0:
            tick(i + 1)
    ci = np.quantile(draws, [.025, .975]).tolist()
    passed = point >= .95 and ci[0] >= .90
    return {
        "status": "passed" if passed else "failed",
        "pass": bool(passed),
        "control_count": int(control.sum()),
        "female_controls": int((labels == 0).sum()),
        "male_controls": int((labels == 1).sum()),
        "AUROC": point,
        "CI95": ci,
        "BA_at_0_5_descriptive": float(balanced_accuracy_score(labels, scores >= .5)),
        "bootstrap_repetitions": N_BOOTSTRAPS,
        "bootstrap_seed": SEX_BOOTSTRAP_SEED,
    }


def disease_summary(y: np.ndarray, donors: np.ndarray, predictions: dict[str, np.ndarray]) -> dict:
    metrics = {
        name: {
            "AUROC": float(roc_auc_score(y, predictions[name])),
            "log_loss": float(log_loss(y, predictions[name], labels=[0, 1])),
        } for name in MODELS
    }
    if np.array_equal(predictions["M5"], predictions["M4"]):
        primary = {
            "branch": "M5_identical_to_M4",
            "delta_AUROC": 0.0,
            "CI95": [0.0, 0.0],
            "confirmatory_margin_interpretation": "not_applicable_identical_model_branch",
        }
    else:
        primary, _draws = paired_bootstrap(
            y, predictions["M4"], predictions["M5"], donors,
            repetitions=N_BOOTSTRAPS, seed=PRIMARY_BOOTSTRAP_SEED,
        )
        primary["branch"] = "distinct_models"
        primary["confirmatory_margin_interpretation"] = (
            "upper_bound_below_0_05" if primary["CI95"][1] < .05 else "0_05_not_excluded"
        )
    return {"metrics": metrics, "primary": primary}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--published-lock-manifest", required=True, type=Path)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--final-dir", required=True, type=Path)
    parser.add_argument("--marker-manifest", required=True, type=Path)
    parser.add_argument("--sex-artifact", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    progress(0, 6, started, "checking published lock and frozen files")
    private_root = (Path.home() / "PPMI_private").resolve()
    output_dir = args.output_dir.resolve()
    if output_dir != private_root and private_root not in output_dir.parents:
        raise ValueError("PPMI outputs must remain under ~/PPMI_private")
    if output_dir.exists():
        raise ValueError("Private output directory already exists; use a fresh path")

    lock = json.loads(args.published_lock_manifest.read_text(encoding="utf-8"))
    verify_lock(lock, args.published_lock_manifest, args.bundle, args.final_dir,
                args.marker_manifest, args.sex_artifact)
    progress(1, 6, started, "published-lock manifest and file hashes checked")
    binding = freeze(args.bundle, args.final_dir / "preprocessing.npz")
    for key in ("ordered_population_sha256", "raw_matrix_sha256",
                "scored_matrix_geo_array_hash_sha256"):
        if binding[key] != lock[key]:
            raise ValueError(f"PPMI input differs from published second lock: {key}")
    if binding["group_counts"] != lock["final_population_counts"]:
        raise ValueError("PPMI group counts differ from published second lock")
    progress(2, 6, started, "population and both matrix hashes checked")
    ranks, age, sex, y, donors, genes = load_locked_input(args.bundle, binding)
    scorer = FrozenGEOScorer(args.final_dir, args.marker_manifest, args.sex_artifact)
    progress(3, 6, started, "locked input and GEO scorer loaded")

    os.umask(0o077)
    output_dir.mkdir(mode=0o700, parents=False)
    control_predictions = np.full(len(y), np.nan, dtype=np.float64)
    control_predictions[y == 0] = scorer.predict_sex_control(ranks[y == 0], genes)
    private_sex_file = output_dir / "PRIVATE_ROW_ORDERED_SEX_CONTROL.npz"
    np.savez_compressed(private_sex_file, recorded_sex=sex, historical_group=y,
                        sex_control=control_predictions)
    private_sex_file.chmod(0o600)
    gate = sex_gate(y, sex, control_predictions,
                    tick=lambda count: progress(3 + count / N_BOOTSTRAPS, 6, started,
                                                f"sex-control bootstrap {count}/{N_BOOTSTRAPS}"))
    progress(4, 6, started, "sex-control gate complete")
    aggregate = {
        "status": "sex_gate_passed" if gate["pass"] else "sex_gate_failed_or_unevaluable",
        "zenodo_doi": lock["zenodo_doi"],
        "ordered_population_sha256": binding["ordered_population_sha256"],
        "raw_matrix_sha256": binding["raw_matrix_sha256"],
        "scored_matrix_geo_array_hash_sha256": binding["scored_matrix_geo_array_hash_sha256"],
        "sex_gate": gate,
        "PD_model_scored": False,
    }
    if gate["pass"]:
        model_predictions = scorer.predict_models(ranks, age, sex, genes)
        aggregate["disease"] = disease_summary(y, donors, model_predictions)
        aggregate["PD_model_scored"] = True
        aggregate["status"] = "completed_after_sex_gate"
        private_file = output_dir / "PRIVATE_ROW_ORDERED_PREDICTIONS.npz"
        np.savez_compressed(private_file, y=y, sex=sex, sex_control=control_predictions,
                            **model_predictions)
        private_file.chmod(0o600)
        progress(5, 6, started, "fixed M1-M5 and paired bootstrap complete")
    else:
        progress(5, 6, started, "PD scoring suspended by sex-control gate")
    aggregate_file = output_dir / "AGGREGATE_RESULTS.json"
    aggregate_file.write_text(json.dumps(aggregate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    aggregate_file.chmod(0o600)
    progress(6, 6, started, "owner-only outputs finalized")
    print(json.dumps({"status": aggregate["status"], "PD_model_scored": aggregate["PD_model_scored"],
                      "aggregate_output": str(aggregate_file)}))


if __name__ == "__main__":
    main()
