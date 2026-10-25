"""Shared safeguards and numerical helpers for the draft E1–E10 analyses.

Importing this module never reads PPMI or GEO data. Real-data entry points must
call ``require_execution_authorization`` before loading any analysis input.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess

import numpy as np

MASTER_SEED = 20261014
PPMI_PRIVATE = Path.home() / "PPMI_private"
FREETSA_SIGNER_SHA256 = "8bfb0305bb64e2571ca507552ef3245cb1c2fee8728e0ff8689225081ea13467"
FREETSA_CA_SHA256 = "2151b61137ffa86bf664691ba67e7da0b19f98c758e3d228d5d8ebf27e044438"


def rng_for(number: int, repetition: int = 0) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([MASTER_SEED, number, repetition]))


def require_execution_authorization(record: Path, *, analysis: str) -> None:
    """Require an independently timestamped amendment and pre-run lock audit.

    A private hash-only timestamp is an alternative provenance anchor while
    public methods-only posting awaits DPC advice. It is not a DOI or approval.
    E3 and E9 retain their separate access/eligibility evidence gate.
    """
    obj = json.loads(Path(record).read_text(encoding="utf-8"))
    if (obj.get("phase") != 3 or obj.get("real_data_execution_allowed") is not True
            or analysis not in obj.get("approved_analyses", [])):
        raise PermissionError(f"Real-data {analysis} execution gate is not satisfied")
    if analysis in {"E3", "E9"} and not obj.get("approved_application_evidence_sha256"):
        raise PermissionError(f"Real-data {analysis} needs separate scope evidence")

    receipt = Path(obj.get("private_lock_audit_receipt", ""))
    if (not receipt.is_file() or sha256(receipt) != obj.get("private_lock_audit_sha256")):
        raise PermissionError("Full private-artifact lock audit receipt is missing or changed")
    audit = json.loads(receipt.read_text(encoding="utf-8"))
    if (audit.get("result") != "PASS" or audit.get("reference_hash_match") is not True
            or not isinstance(audit.get("verified_file_count"), int)
            or audit["verified_file_count"] < 1 or audit.get("mismatches")):
        raise PermissionError("Full private-artifact lock audit has not passed")

    if obj.get("amendment_registration_doi"):
        return
    manifest = Path(obj.get("timestamped_manifest", ""))
    response = Path(obj.get("timestamp_response", ""))
    ca = Path(obj.get("timestamp_ca_certificate", ""))
    tsa = Path(obj.get("timestamp_tsa_certificate", ""))
    if not all(path.is_file() for path in (manifest, response, ca, tsa)):
        raise PermissionError("Independent hash-only timestamp receipt is incomplete")
    if sha256(ca) != FREETSA_CA_SHA256 or sha256(tsa) != FREETSA_SIGNER_SHA256:
        raise PermissionError("Timestamp trust certificates differ from the pinned issuer")
    manifest_obj = json.loads(manifest.read_text(encoding="utf-8"))
    if manifest_obj.get("kind") != "post_primary_methods_and_code_hashes_v1":
        raise PermissionError("Timestamp manifest has an unexpected format")
    root = Path(__file__).resolve().parents[1]
    files = manifest_obj.get("files", {})
    if (not isinstance(files, dict) or
            "revision/PHASE_2_AMENDMENT_DRAFT.md" not in files or
            "revision/postlock_common.py" not in files):
        raise PermissionError("Timestamp manifest omits the amendment or execution gate")
    for name, expected in files.items():
        path = (root / name).resolve()
        if not path.is_relative_to(root) or not path.is_file() or sha256(path) != expected:
            raise PermissionError(f"Timestamped file changed: {name}")
    try:
        verified = subprocess.run(
            ["openssl", "ts", "-verify", "-data", str(manifest), "-in", str(response),
             "-CAfile", str(ca), "-untrusted", str(tsa)],
            capture_output=True, text=True, timeout=30, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise PermissionError("Timestamp certificate could not be verified") from exc
    if verified.returncode != 0:
        raise PermissionError("Timestamp certificate verification failed")


def private_output(path: Path) -> Path:
    path = Path(path).expanduser().resolve()
    private = PPMI_PRIVATE.resolve()
    if not path.is_relative_to(private):
        raise ValueError("PPMI diagnostic outputs must remain under ~/PPMI_private")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def write_aggregate(path: Path, obj: dict) -> None:
    path = private_output(path)
    data = json.dumps(obj, indent=2, sort_keys=True, allow_nan=False).encode("utf-8")
    temp = path.with_name(path.name + ".tmp")
    fd = os.open(temp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        temp.replace(path)
    finally:
        if temp.exists():
            temp.unlink()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def summary(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError("A one-dimensional finite vector is required")
    return {"n": int(len(values)),
            "mean": float(values.mean()) if len(values) else None,
            "sample_sd": float(values.std(ddof=1)) if len(values) > 1 else None}


def percentile_interval(values: np.ndarray) -> list[float] | None:
    values = np.asarray(values, dtype=float)
    if len(values) == 0 or not np.isfinite(values).all():
        return None
    return np.quantile(values, [0.025, 0.975]).astype(float).tolist()


def stratified_indices(y: np.ndarray, draws: int, number: int) -> np.ndarray:
    y = np.asarray(y)
    if len(np.unique(y)) != 2:
        raise ValueError("Both outcome classes are required")
    groups = [np.flatnonzero(y == level) for level in sorted(np.unique(y))]
    rng = rng_for(number)
    return np.vstack([np.concatenate([rng.choice(g, len(g), replace=True)
                                      for g in groups]) for _ in range(draws)])


def synthetic_ranks(n: int = 20, g: int = 30, seed: int = 7) -> np.ndarray:
    rng = np.random.default_rng(seed)
    counts = rng.normal(size=(n, g))
    # Distinct continuous values: deterministic rank lattice with row mean 0.5.
    return (np.argsort(np.argsort(counts, axis=1), axis=1) + .5) / g
