"""Run the registered GEO age/sex residualization sensitivity in all frozen scopes.

This wrapper uses the locked nested-CV/transfer pipeline. It changes only the
added gene block, instruments model.predict warnings with fit context, and
scores each outer holdout only after inner model selection is complete.
"""
from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import geo_pipeline as geo
from pipeline.pilot_geo_age_sex_residualization import (
    TrainOLSResidualPrep,
    _PREP_FIT_AUDIT,
    _PREP_MAX_NORMAL_EQ_AUDIT,
)
from pipeline.prediction_warning_capture import instrument_predict

FULL_OUT = ROOT / "private_data" / "geo_age_sex_resid_full_20261019"
RELEASE = ROOT / "results/geo_development_20261018/GEO_CV_RELEASE_MANIFEST_20261018.json"
AUDIT = ROOT / "data/manifests/input_audit.json"
DECISION = ROOT / "docs/GEO_AGE_SEX_RESIDUALIZATION_DECISION_20261019_EN.md"
WARNING_LOG = FULL_OUT / "prediction_warnings.jsonl"
ORIGINAL_PREP = geo.Prep
ORIGINAL_CHECKED_FIT = geo.checked_fit


def _instrumented_checked_fit(*args, **kwargs):
    fit, records = ORIGINAL_CHECKED_FIT(*args, **kwargs)
    context = dict(args[5] if len(args) > 5 else kwargs["context"])
    context.setdefault("model", args[2] if len(args) > 2 else kwargs["model"])
    context.setdefault("lambda", args[3] if len(args) > 3 else kwargs["lam"])
    context.setdefault("rho", args[4] if len(args) > 4 else kwargs["rho"])
    return instrument_predict(fit, context, WARNING_LOG), records


def _scope_inputs(data, split):
    lookup = {sample: i for i, sample in enumerate(data.sample)}
    train = data.take(np.array([lookup[s] for s in split["train"]], dtype=int))
    heldout = data.take(np.array([lookup[s] for s in split["test"]], dtype=int))
    local = {sample: i for i, sample in enumerate(train.sample)}
    inner = [(np.array([local[s] for s in fold["train"]], dtype=int),
              np.array([local[s] for s in fold["validation"]], dtype=int))
             for fold in split["inner"]]
    return train, heldout, inner


def _warning_count(scope):
    if not WARNING_LOG.exists():
        return 0
    return sum(json.loads(line)["scope"] == scope
               for line in WARNING_LOG.read_text(encoding="utf-8").splitlines())


def _numerical_gate_summary(records):
    """Require the final permitted attempt for every candidate to pass KKT."""
    groups = defaultdict(list)
    for record in records:
        context = tuple((key, json.dumps(record.get(key), sort_keys=True)) for key in (
            "scope", "stage", "fold", "model", "lambda", "rho",
            "expansion_round", "baseline_only", "preprocessing_sha256",
            "train_sample_sha256"))
        groups[context].append(record)
    final_records = [max(group, key=lambda record: record["attempt"])
                     for group in groups.values()]
    if not all(record.get("passed") is True for record in final_records):
        raise AssertionError("At least one candidate exhausted its KKT attempts")
    return len(final_records)


def _run_scope(data, split):
    scope = split["scope"]
    scope_dir = FULL_OUT / scope.replace("/", "__")
    scope_dir.mkdir()
    train, heldout, inner = _scope_inputs(data, split)
    heldout_samples = set(heldout.sample.tolist())
    heldout_donors = set(heldout.donor.tolist())

    primary_prep = ORIGINAL_PREP(train)
    residual_prep = TrainOLSResidualPrep(train)
    baseline_checks = {}
    for label, subset in (("train", train), ("heldout", heldout)):
        for model in ("M1", "M4"):
            baseline_checks[f"{label}_{model}_identical"] = bool(np.array_equal(
                primary_prep.transform(subset, model)[0],
                residual_prep.transform(subset, model)[0]))
    if not all(baseline_checks.values()):
        raise AssertionError(f"Fixed baseline changed in {scope}")

    allowed_fit_hashes = {geo.array_hash(train.sample)}
    expected_prep_hashes = []
    for train_indices, valid_indices in inner:
        fold_train = train.take(train_indices)
        fold_valid = train.take(valid_indices)
        expected_prep_hashes.append(geo.array_hash(fold_train.sample))
        allowed_fit_hashes.add(geo.array_hash(fold_train.sample))
        if set(fold_train.donor.tolist()) & set(fold_valid.donor.tolist()):
            raise AssertionError(f"Inner donor overlap in {scope}")
        if heldout_samples & set(fold_train.sample.tolist()):
            raise AssertionError(f"Outer holdout sample entered inner training in {scope}")
        if heldout_donors & set(fold_train.donor.tolist()):
            raise AssertionError(f"Outer holdout donor entered inner training in {scope}")
    expected_prep_hashes.append(geo.array_hash(train.sample))
    if heldout_samples & set(train.sample.tolist()):
        raise AssertionError(f"Outer holdout sample entered outer training in {scope}")
    if heldout_donors & set(train.donor.tolist()):
        raise AssertionError(f"Outer holdout donor entered outer training in {scope}")

    _PREP_FIT_AUDIT.clear()
    _PREP_MAX_NORMAL_EQ_AUDIT.clear()
    started, cpu_started = time.monotonic(), time.process_time()
    prep, fits, decisions = geo.train_suite(train, inner, scope_dir, scope)
    # The outer holdout is first transformed/scored after inner selection and
    # all final refits have completed.
    probabilities = (geo.predict_suite(prep, fits, decisions, heldout)
                     if len(heldout.y) else {})
    if probabilities:
        np.savez_compressed(scope_dir / "heldout_predictions.npz",
                            samples=heldout.sample, y=heldout.y, **probabilities)

    fit_records = [json.loads(line) for line in
                   (scope_dir / "fits.jsonl").read_text(encoding="utf-8").splitlines()]
    fit_hashes = {record["train_sample_sha256"] for record in fit_records}
    if not fit_hashes <= allowed_fit_hashes:
        raise AssertionError(f"A fit used samples outside frozen training partitions: {scope}")
    successful_fit_count = _numerical_gate_summary(fit_records)
    if _PREP_FIT_AUDIT != expected_prep_hashes:
        raise AssertionError(
            f"OLS preprocessing was not fit once in each inner-training partition "
            f"and once in outer training for {scope}")
    if len(_PREP_MAX_NORMAL_EQ_AUDIT) != len(expected_prep_hashes):
        raise AssertionError(f"Missing OLS normal-equation audit in {scope}")
    max_normal_equation = max(_PREP_MAX_NORMAL_EQ_AUDIT)
    if max_normal_equation > 1e-7:
        raise AssertionError(f"OLS normal-equation residual too large in {scope}")

    warning_count = _warning_count(scope)
    result = dict(
        status="complete", scope=scope,
        seconds=time.monotonic() - started,
        process_CPU_seconds=time.process_time() - cpu_started,
        inner_fold_count=len(inner), fit_attempt_record_count=len(fit_records),
        successful_candidate_fit_count=successful_fit_count,
        all_candidate_fits_passed_numerical_gate=True,
        fit_training_hashes_match_frozen_training_partitions=True,
        residualization_fit_count=len(_PREP_FIT_AUDIT),
        residualization_fit_hashes_match_inner_and_outer_training=True,
        residualization_fit_sample_hashes=list(_PREP_FIT_AUDIT),
        max_training_OLS_normal_equation_residual=max_normal_equation,
        baseline_transform_exactly_unchanged=baseline_checks,
        residualization_precedes_source_stratified_MAD_selection=True,
        residualization_predictors=[
            "intercept", "training-standardized age", "age squared", "sex",
            "age-by-sex"],
        selected_gene_count=int(len(prep.selected)),
        pooled_OLS_uses_one_source_free_training_fit=True,
        heldout_predictions_generated_after_inner_selection=bool(probabilities),
        heldout_metrics={name: geo.metrics(heldout.y, pred)
                         for name, pred in probabilities.items()},
        prediction_warning_count=warning_count,
        PPMI_loaded=False, PPMI_scores_computed=False)
    geo.write_json(scope_dir / "RESULTS.json", result)
    print(f"{scope} complete; fits: {len(fit_records)}; warnings: {warning_count}",
          flush=True)
    return result


def main():
    if FULL_OUT.exists():
        raise FileExistsError(f"Refusing to overwrite prior run: {FULL_OUT}")
    release = json.loads(RELEASE.read_text(encoding="utf-8"))
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    decision_hash = geo.digest(DECISION)
    pipeline_hash = geo.digest(Path(geo.__file__))
    if pipeline_hash != release["pipeline_code_sha256"]:
        raise ValueError("Locked primary pipeline hash differs from published release")
    data, input_hashes = geo.load_geo()
    if input_hashes != audit["input_sha256"]:
        raise ValueError("Prepared GEO inputs differ from their frozen audit")
    splits = geo.split_manifest(data, input_hashes)
    split_hash = hashlib.sha256(json.dumps(
        splits, indent=2, ensure_ascii=False, allow_nan=False).encode()).hexdigest()
    if split_hash != release["split_manifest_sha256"]:
        raise ValueError("Reconstructed split manifest differs from published release")
    all_splits = splits["splits"]
    if len(all_splits) != 13:
        raise ValueError(f"Expected 13 registered GEO scopes, found {len(all_splits)}")

    FULL_OUT.mkdir(parents=True)
    total_started, total_cpu_started = time.monotonic(), time.process_time()
    scope_results = {}
    original_prep, original_checked_fit = geo.Prep, geo.checked_fit
    geo.Prep = TrainOLSResidualPrep
    geo.checked_fit = _instrumented_checked_fit
    try:
        for index, split in enumerate(all_splits, start=1):
            print(f"[{index}/{len(all_splits)}] starting {split['scope']}", flush=True)
            scope_result = _run_scope(data, split)
            scope_results[split["scope"]] = scope_result
    finally:
        geo.Prep = original_prep
        geo.checked_fit = original_checked_fit

    warning_count = (sum(1 for _ in WARNING_LOG.open(encoding="utf-8"))
                     if WARNING_LOG.exists() else 0)
    result = dict(
        status="complete", sensitivity="registered_train_only_age_sex_OLS_residualization",
        scope_count=len(scope_results), scopes=scope_results,
        seconds=time.monotonic() - total_started,
        process_CPU_seconds=time.process_time() - total_cpu_started,
        fit_attempt_record_count=sum(item["fit_attempt_record_count"]
                                     for item in scope_results.values()),
        successful_candidate_fit_count=sum(item["successful_candidate_fit_count"]
                                           for item in scope_results.values()),
        all_candidate_fits_passed_numerical_gate=all(
            item["all_candidate_fits_passed_numerical_gate"]
            for item in scope_results.values()),
        all_training_hashes_match_frozen_partitions=all(
            item["fit_training_hashes_match_frozen_training_partitions"]
            for item in scope_results.values()),
        all_residualization_partitions_match_expected=all(
            item["residualization_fit_hashes_match_inner_and_outer_training"]
            for item in scope_results.values()),
        all_baseline_transforms_exactly_unchanged=all(
            all(item["baseline_transform_exactly_unchanged"].values())
            for item in scope_results.values()),
        total_prediction_warning_count=warning_count,
        original_pipeline_sha256=pipeline_hash,
        published_split_manifest_sha256=split_hash,
        implementation_decision_sha256=decision_hash,
        input_sha256=input_hashes, python=sys.version, platform=platform.platform(),
        PPMI_loaded=False, PPMI_scores_computed=False,
        interpretation="Registered GEO residualization sensitivity; no PPMI data were accessed.")
    geo.write_json(FULL_OUT / "RESULTS.json", result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
