"""One-scope leakage/runtime pilot for registered train-only age/sex OLS residualization."""
from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import geo_pipeline as geo
from pipeline.prediction_warning_capture import instrument_predict
_ORIGINAL_CHECKED_FIT = geo.checked_fit

SCOPE = "GSE57475/outer/0"
OUT = ROOT / "private_data" / "geo_age_sex_resid_pilot_20261019"
RELEASE = ROOT / "results/geo_development_20261018/GEO_CV_RELEASE_MANIFEST_20261018.json"
AUDIT = ROOT / "data/manifests/input_audit.json"
WARNING_LOG = OUT / "prediction_warnings.jsonl"
_PREP_FIT_AUDIT = []
_PREP_MAX_NORMAL_EQ_AUDIT = []


class TrainOLSResidualPrep(geo.Prep):
    """M1-basis OLS residualization, fitted separately inside every training partition."""

    def __init__(self, train, max_genes=1000):
        # Retain the primary age/sex and marker baseline preprocessing exactly.
        super().__init__(train, max_genes)
        design = _design(train, self.age_mean, self.age_sd)
        self.resid_beta = np.linalg.lstsq(design, train.X, rcond=None)[0]
        residual = train.X - design @ self.resid_beta
        variability = []
        for source in np.unique(train.source):
            x = residual[train.source == source]
            variability.append(np.median(np.abs(x - np.median(x, axis=0)), axis=0))
        mad = np.min(variability, axis=0)
        eligible = np.flatnonzero(mad > 0)
        order = np.lexsort((train.genes[eligible], -mad[eligible]))
        self.selected = eligible[order[:max_genes]]
        selected = residual[:, self.selected]
        self.gene_mean, self.gene_sd = selected.mean(axis=0), selected.std(axis=0)
        if not len(self.selected) or np.any(self.gene_sd <= 0):
            raise ValueError("Residualized positive-MAD selection is empty or constant")
        self.fit_samples_hash = geo.array_hash(train.sample)
        _PREP_FIT_AUDIT.append(self.fit_samples_hash)
        self.train_ols_max_abs_normal_equation = float(
            np.max(np.abs(design.T @ residual)))
        _PREP_MAX_NORMAL_EQ_AUDIT.append(
            self.train_ols_max_abs_normal_equation)

    def transform(self, data, model):
        raw = self.baseline_raw(data)
        z = (raw - self.base_mean) / np.where(self.base_sd > 0, self.base_sd, 1.)
        keep = [i for i in self.keep if i < 4] if model in ("M1", "M3") else self.keep
        baseline = z[:, keep] if model != "M2" else np.empty((len(z), 0))
        residual = data.X[:, self.selected] - _design(
            data, self.age_mean, self.age_sd) @ self.resid_beta[:, self.selected]
        gene = (residual - self.gene_mean) / self.gene_sd
        if model in ("M1", "M4"):
            return np.ascontiguousarray(baseline), baseline.shape[1]
        return np.ascontiguousarray(np.column_stack((baseline, gene))), baseline.shape[1]

    def fingerprint(self):
        return geo.array_hash(np.array([self.age_mean, self.age_sd]), self.base_mean,
            self.base_sd, np.array(self.keep), self.selected, self.gene_mean,
            self.gene_sd, self.resid_beta[:, self.selected])

    def save(self, path):
        np.savez_compressed(path, age_mean=self.age_mean, age_sd=self.age_sd,
            base_mean=self.base_mean, base_sd=self.base_sd, keep=self.keep,
            selected=self.selected, gene_mean=self.gene_mean, gene_sd=self.gene_sd,
            genes=self.genes, fit_samples_hash=self.fit_samples_hash,
            resid_beta=self.resid_beta[:, self.selected],
            residualization_columns=np.array(
                ["intercept", "standardized_age", "standardized_age_squared",
                 "sex", "standardized_age_by_sex"]))


def _design(data, age_mean, age_sd):
    age = (data.age - age_mean) / age_sd
    return np.column_stack((np.ones(len(data.y)), age, age ** 2,
                            data.sex, age * data.sex))


def _instrumented_checked_fit(*args, **kwargs):
    fit, records = _ORIGINAL_CHECKED_FIT(*args, **kwargs)
    context = dict(args[5] if len(args) > 5 else kwargs["context"])
    context.setdefault("model", args[2] if len(args) > 2 else kwargs["model"])
    context.setdefault("lambda", args[3] if len(args) > 3 else kwargs["lam"])
    context.setdefault("rho", args[4] if len(args) > 4 else kwargs["rho"])
    return instrument_predict(fit, context, WARNING_LOG), records


def main():
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite prior run: {OUT}")
    OUT.mkdir(parents=True)
    release = json.loads(RELEASE.read_text())
    audit = json.loads(AUDIT.read_text())
    pipeline_hash = geo.digest(Path(geo.__file__))
    if pipeline_hash != release["pipeline_code_sha256"]:
        raise ValueError("Locked primary pipeline hash differs from the published release")
    data, hashes = geo.load_geo()
    if hashes != audit["input_sha256"]:
        raise ValueError("Prepared GEO inputs differ from their frozen audit")
    splits = geo.split_manifest(data, hashes)
    split_hash = hashlib.sha256(json.dumps(
        splits, indent=2, ensure_ascii=False, allow_nan=False).encode()).hexdigest()
    if split_hash != release["split_manifest_sha256"]:
        raise ValueError("Reconstructed split manifest differs from the published release")
    found = [s for s in splits["splits"] if s["scope"] == SCOPE]
    if len(found) != 1:
        raise ValueError("Frozen outer scope was not unique")
    split = found[0]
    lookup = {sample: i for i, sample in enumerate(data.sample)}
    train = data.take(np.array([lookup[x] for x in split["train"]], dtype=int))
    heldout = data.take(np.array([lookup[x] for x in split["test"]], dtype=int))
    local = {sample: i for i, sample in enumerate(train.sample)}
    inner = [(np.array([local[x] for x in f["train"]], dtype=int),
              np.array([local[x] for x in f["validation"]], dtype=int))
             for f in split["inner"]]

    primary_prep = geo.Prep(train)
    residual_prep = TrainOLSResidualPrep(train)
    baseline_checks = {}
    for label, subset in (("train", train), ("heldout", heldout)):
        for model in ("M1", "M4"):
            baseline_checks[f"{label}_{model}_identical"] = np.array_equal(
                primary_prep.transform(subset, model)[0],
                residual_prep.transform(subset, model)[0])
    if not all(baseline_checks.values()):
        raise AssertionError("Residualization changed the fixed baseline block")
    if residual_prep.train_ols_max_abs_normal_equation > 1e-7:
        raise AssertionError("Training residuals do not satisfy the OLS normal equations")

    allowed_fit_hashes = {geo.array_hash(train.sample)}
    expected_prep_hashes = [geo.array_hash(train.sample[ti]) for ti, _ in inner]
    expected_prep_hashes.append(geo.array_hash(train.sample))
    for ti, _ in inner:
        allowed_fit_hashes.add(geo.array_hash(train.sample[ti]))
        if set(heldout.sample) & set(train.sample[ti]):
            raise AssertionError("Outer holdout appears in an inner training partition")
    if set(heldout.sample) & set(train.sample):
        raise AssertionError("Outer holdout appears in the outer training partition")

    started, cpu_started = time.monotonic(), time.process_time()
    original_prep = geo.Prep
    original_checked_fit = geo.checked_fit
    _PREP_FIT_AUDIT.clear()
    geo.Prep = TrainOLSResidualPrep
    geo.checked_fit = _instrumented_checked_fit
    try:
        _, _, _ = geo.train_suite(train, inner, OUT, SCOPE)
    finally:
        geo.Prep = original_prep
        geo.checked_fit = original_checked_fit
    seconds = time.monotonic() - started
    cpu_seconds = time.process_time() - cpu_started
    records = [json.loads(line) for line in (OUT / "fits.jsonl").read_text().splitlines()]
    if not {r["train_sample_sha256"] for r in records} <= allowed_fit_hashes:
        raise AssertionError("A fit used samples outside a frozen training partition")
    if _PREP_FIT_AUDIT != expected_prep_hashes:
        raise AssertionError(
            "Residualization was not fit exactly once per inner-training partition "
            "and once on the outer training partition")
    warning_count = (sum(1 for _ in WARNING_LOG.open(encoding="utf-8"))
                     if WARNING_LOG.exists() else 0)

    result = dict(
        status="timing_pilot_complete", sensitivity="train_only_age_sex_OLS_residualization",
        scope=SCOPE, seconds=seconds, process_CPU_seconds=cpu_seconds,
        inner_fold_count=len(inner), checked_fit_count=len(records),
        checked_fit_training_hashes_match_frozen_train_partitions=True,
        outer_holdout_used_for_performance=False,
        baseline_transform_exactly_unchanged=baseline_checks,
        residualization_predictors=[
            "intercept", "training-standardized age", "age squared", "sex",
            "age-by-sex"],
        residualization_precedes_source_stratified_MAD_selection=True,
        residualization_fit_partition_count=len(_PREP_FIT_AUDIT),
        residualization_fit_sample_hashes=_PREP_FIT_AUDIT,
        expected_residualization_fit_sample_hashes=expected_prep_hashes,
        prediction_warning_log=str(WARNING_LOG),
        prediction_warning_count=warning_count,
        train_ols_max_abs_normal_equation=residual_prep.train_ols_max_abs_normal_equation,
        selected_gene_count=int(len(residual_prep.selected)),
        original_pipeline_sha256=pipeline_hash,
        published_split_manifest_sha256=split_hash,
        input_sha256=hashes, python=sys.version, platform=platform.platform(),
        PPMI_loaded=False, PPMI_scores_computed=False, metrics_reported=False,
        interpretation="One-scope runtime/leakage pilot only; not a sensitivity effect estimate.")
    geo.write_json(OUT / "RESULTS.json", result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
