"""Shared safeguards and numerical helpers for the draft E1–E10 analyses.

Importing this module never reads PPMI or GEO data. Real-data entry points must
call ``require_execution_authorization`` before loading any analysis input.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np

MASTER_SEED = 20261014
PPMI_PRIVATE = Path.home() / "PPMI_private"


def rng_for(number: int, repetition: int = 0) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([MASTER_SEED, number, repetition]))


def require_execution_authorization(record: Path, *, analysis: str) -> None:
    """Fail closed until a separate, explicit PPMI application/DPC gate record exists.

    A record is not created by the analysis scripts. Its SHA-256 and the
    underlying approval wording must be verified separately before phase 3.
    """
    obj = json.loads(Path(record).read_text(encoding="utf-8"))
    if (obj.get("phase") != 3 or obj.get("real_data_execution_allowed") is not True
            or analysis not in obj.get("approved_analyses", [])
            or not obj.get("approved_application_evidence_sha256")
            or not obj.get("amendment_registration_doi")):
        raise PermissionError(f"Real-data {analysis} execution gate is not satisfied")


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
