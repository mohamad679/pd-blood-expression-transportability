"""One-scope runtime pilot for the registered GEO X/Y-removal sensitivity.

This wrapper leaves the locked primary pipeline source unchanged. It runs the
registered training suite on a frozen outer scope, restricting only the added
gene block to autosomes. It deliberately does not score the outer holdout.
"""
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


SCOPE = "GSE57475/outer/0"
OUT = ROOT / "private_data" / "geo_xy_sensitivity_timing_20261019"
FULL_OUT = ROOT / "private_data" / "geo_xy_sensitivity_full_20261019"
RELEASE = ROOT / "results/geo_development_20261018/GEO_CV_RELEASE_MANIFEST_20261018.json"
AUDIT = ROOT / "data/manifests/input_audit.json"
WARNING_LOG = FULL_OUT / "prediction_warnings.jsonl"
_ORIGINAL_CHECKED_FIT = geo.checked_fit


def _instrumented_checked_fit(*args, **kwargs):
    fit, records = _ORIGINAL_CHECKED_FIT(*args, **kwargs)
    context = dict(args[5] if len(args) > 5 else kwargs["context"])
    context.setdefault("model", args[2] if len(args) > 2 else kwargs["model"])
    context.setdefault("lambda", args[3] if len(args) > 3 else kwargs["lam"])
    context.setdefault("rho", args[4] if len(args) > 4 else kwargs["rho"])
    return instrument_predict(fit, context, WARNING_LOG), records


class AutosomalPrep(geo.Prep):
    """Keep primary baseline features; select/scalе the added block on train autosomes."""

    def __init__(self, train, max_genes=1000):
        super().__init__(train, max_genes)
        autosomal = np.isin(_CHROMOSOMES, [f"chr{i}" for i in range(1, 23)])
        if len(autosomal) != len(train.genes):
            raise ValueError("Chromosome annotation is not aligned to the frozen gene matrix")
        eligible = np.flatnonzero(autosomal)
        variability = []
        for source in np.unique(train.source):
            x = train.X[train.source == source][:, eligible]
            variability.append(np.median(np.abs(x - np.median(x, axis=0)), axis=0))
        mad = np.min(variability, axis=0)
        eligible = eligible[mad > 0]
        mad = mad[mad > 0]
        order = np.lexsort((train.genes[eligible], -mad))
        self.selected = eligible[order[:max_genes]]
        genes = train.X[:, self.selected]
        self.gene_mean, self.gene_sd = genes.mean(axis=0), genes.std(axis=0)
        if np.any(self.gene_sd <= 0):
            raise ValueError("Positive-MAD autosomal selection included a constant gene")
        self.fit_samples_hash = geo.array_hash(train.sample)


def scope_inputs(data, split):
    lookup = {sample: i for i, sample in enumerate(data.sample)}
    train = data.take(np.array([lookup[x] for x in split["train"]], dtype=int))
    heldout = data.take(np.array([lookup[x] for x in split["test"]], dtype=int))
    local = {sample: i for i, sample in enumerate(train.sample)}
    inner = [(np.array([local[x] for x in f["train"]], dtype=int),
              np.array([local[x] for x in f["validation"]], dtype=int))
             for f in split["inner"]]
    return train, heldout, inner


def main():
    global _CHROMOSOMES
    if FULL_OUT.exists():
        raise FileExistsError(f"Refusing to overwrite prior run: {FULL_OUT}")
    FULL_OUT.mkdir(parents=True)
    release = json.loads(RELEASE.read_text())
    audit = json.loads(AUDIT.read_text())
    original_pipeline_hash = geo.digest(Path(geo.__file__))
    if original_pipeline_hash != release["pipeline_code_sha256"]:
        raise ValueError("Locked primary pipeline hash differs from the published release")

    data, hashes = geo.load_geo()
    if hashes != audit["input_sha256"]:
        raise ValueError("Prepared GEO inputs differ from their frozen audit")
    splits = geo.split_manifest(data, hashes)
    split_hash = hashlib.sha256(json.dumps(
        splits, indent=2, ensure_ascii=False, allow_nan=False).encode()).hexdigest()
    if split_hash != release["split_manifest_sha256"]:
        raise ValueError("Reconstructed split manifest differs from the published release")
    with np.load(ROOT / "prepared/GSE57475_development.npz", allow_pickle=False) as z:
        _CHROMOSOMES = z["chromosomes"]
    if len(_CHROMOSOMES) != len(data.genes):
        raise ValueError("Chromosome annotation length mismatch")
    if int(np.isin(_CHROMOSOMES, [f"chr{i}" for i in range(1, 23)]).sum()) != 15720:
        raise ValueError("Expected exactly 15,720 autosomal genes")

    all_splits = splits["splits"]
    if len(all_splits) != 13:
        raise ValueError(f"Expected 13 registered GEO scopes; found {len(all_splits)}")
    original_class = geo.Prep
    geo.Prep = AutosomalPrep
    original_checked_fit = geo.checked_fit
    geo.checked_fit = _instrumented_checked_fit
    total_started = time.monotonic()
    total_cpu_started = time.process_time()
    scope_results = {}
    try:
        for split in all_splits:
            scope = split["scope"]
            scope_dir = FULL_OUT / scope.replace("/", "__")
            scope_dir.mkdir()
            train, heldout, inner = scope_inputs(data, split)
            heldout_samples = set(heldout.sample.tolist())

            # The X/Y sensitivity changes only the added gene block; baseline
            # marker/age features must remain bit-identical in every scope.
            original_prep = original_class(train)
            sensitivity_prep = AutosomalPrep(train)
            baseline_checks = {}
            for label, subset in (("train", train), ("heldout", heldout)):
                for model in ("M1", "M4"):
                    before = original_prep.transform(subset, model)[0]
                    after = sensitivity_prep.transform(subset, model)[0]
                    baseline_checks[f"{label}_{model}_identical"] = bool(
                        np.array_equal(before, after))
            if not all(baseline_checks.values()):
                raise AssertionError(f"Baseline changed in {scope}")
            if not np.all(np.isin(_CHROMOSOMES[sensitivity_prep.selected],
                                  [f"chr{i}" for i in range(1, 23)])):
                raise AssertionError(f"Sex-chromosome gene selected in {scope}")

            allowed_fit_hashes = {geo.array_hash(train.sample)}
            for ti, _ in inner:
                allowed_fit_hashes.add(geo.array_hash(train.sample[ti]))
                if heldout_samples & set(train.sample[ti].tolist()):
                    raise AssertionError(f"Outer holdout appears in inner training: {scope}")
            if heldout_samples & set(train.sample.tolist()):
                raise AssertionError(f"Outer holdout appears in outer training: {scope}")

            scope_started = time.monotonic()
            cpu_started = time.process_time()
            prep, fits, decisions = geo.train_suite(train, inner, scope_dir, scope)
            probabilities = geo.predict_suite(prep, fits, decisions, heldout) if len(heldout.y) else {}
            if probabilities:
                np.savez_compressed(scope_dir / "heldout_predictions.npz",
                                    samples=heldout.sample, y=heldout.y, **probabilities)
            records = [json.loads(line) for line in (scope_dir / "fits.jsonl").read_text().splitlines()]
            fit_hashes = {r["train_sample_sha256"] for r in records}
            if not fit_hashes <= allowed_fit_hashes:
                raise AssertionError(f"Fit samples differ from frozen training partitions: {scope}")
            scope_result = dict(
                status="complete", scope=scope,
                seconds=time.monotonic() - scope_started,
                process_CPU_seconds=time.process_time() - cpu_started,
                inner_fold_count=len(inner), checked_fit_count=len(records),
                checked_fit_training_hashes_match_frozen_train_partitions=True,
                baseline_transform_exactly_unchanged=baseline_checks,
                selected_gene_count=int(len(prep.selected)),
                selected_autosomal_only=True,
                heldout_metrics={name: geo.metrics(heldout.y, pred)
                                 for name, pred in probabilities.items()},
                heldout_samples_scored=bool(probabilities),
                PPMI_loaded=False, PPMI_scores_computed=False)
            scope_results[scope] = scope_result
            geo.write_json(scope_dir / "RESULTS.json", scope_result)
            print(scope, "complete; fits:", len(records), flush=True)
    finally:
        geo.Prep = original_class
        geo.checked_fit = original_checked_fit

    result = dict(
        status="complete", sensitivity="registered_XY_gene_block_removal",
        scope_count=len(scope_results), scopes=scope_results,
        seconds=time.monotonic() - total_started,
        process_CPU_seconds=time.process_time() - total_cpu_started,
        checked_fit_count=sum(x["checked_fit_count"] for x in scope_results.values()),
        all_fit_training_hashes_match_frozen_partitions=all(
            x["checked_fit_training_hashes_match_frozen_train_partitions"]
            for x in scope_results.values()),
        all_baseline_transforms_exactly_unchanged=all(
            all(x["baseline_transform_exactly_unchanged"].values())
            for x in scope_results.values()),
        all_selected_genes_autosomal=all(
            x["selected_autosomal_only"] for x in scope_results.values()),
        original_pipeline_sha256=original_pipeline_hash,
        published_split_manifest_sha256=split_hash,
        input_sha256=hashes, python=sys.version, platform=platform.platform(),
        PPMI_loaded=False, PPMI_scores_computed=False,
        interpretation="Registered GEO X/Y-removal sensitivity; no PPMI data were accessed.")
    geo.write_json(FULL_OUT / "RESULTS.json", result)
    print(json.dumps(result, indent=2), flush=True)


_CHROMOSOMES = None

if __name__ == "__main__":
    main()
