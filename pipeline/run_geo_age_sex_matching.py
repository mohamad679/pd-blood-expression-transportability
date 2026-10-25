"""Run the preregistered full GEO age/sex-matching sensitivity on a local Mac.

This runner performs only metadata matching and GEO model fits. It never loads
PPMI. Run with --plan for a lightweight plan/integrity check before the full
run. Output directories are immutable and a fresh one is required.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import threading
import time
import traceback
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import geo_pipeline as geo
from pipeline.audit_geo_age_sex_matching import match_ids
from pipeline.prediction_warning_capture import instrument_predict

RELEASE = ROOT / "results/geo_development_20261018/GEO_CV_RELEASE_MANIFEST_20261018.json"
INPUT_AUDIT = ROOT / "data/manifests/input_audit.json"
MATCH_AUDIT = ROOT / "results/geo_age_sex_matching_audit_20261019.json"
DECISION = ROOT / "docs/GEO_AGE_SEX_MATCHING_DECISION_20261019_EN.md"
PILOT_DIR = ROOT / "private_data/geo_age_sex_matching_pilot_20261019_GSE99039_outer_1"
OUT = ROOT / "private_data/geo_age_sex_matching_full_20261019"
PILOT_SCOPE = "GSE99039/outer/1"
ELIGIBLE_SCOPES = [
    "GSE99039/outer/1", "GSE99039/outer/2", "GSE99039/outer/3",
    "GSE99039/outer/4", "GSE99039/transfer", "pooled/final",
]
PILOT_SECONDS = 873.1997913669984
PILOT_PAIRS = 70
MIN_CLASS_N = 5
TICK_SECONDS = 30
STOP_MULTIPLE = 2.5

_ORIGINAL_CHECKED_FIT = geo.checked_fit
_ORIGINAL_PREP = geo.Prep
_PREP_FIT_SAMPLE_HASHES: list[str] = []
_STATE_LOCK = threading.Lock()
_STATE: dict = {}
_STOP_EVENT = threading.Event()
_WARNING_LOG: Path | None = None
_DEADLINE = float("inf")


class AuditedPrep(_ORIGINAL_PREP):
    """Record exact sample membership for every fitted preprocessing object."""

    def __init__(self, train, max_genes=1000):
        super().__init__(train, max_genes)
        _PREP_FIT_SAMPLE_HASHES.append(geo.array_hash(train.sample))


def _write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                   encoding="utf-8")
    tmp.replace(path)


def _split_hash(splits: dict) -> str:
    payload = json.dumps(splits, indent=2, ensure_ascii=False, allow_nan=False).encode()
    return hashlib.sha256(payload).hexdigest()


def _indices_for_pairs(data, pair_rows) -> np.ndarray:
    paired_ids = {sample for pair in pair_rows
                  for sample in (pair["case_sample"], pair["control_sample"])}
    return np.asarray([i for i, sample in enumerate(data.sample)
                       if str(sample) in paired_ids], dtype=int)


def _status() -> dict:
    with _STATE_LOCK:
        return dict(_STATE)


def _progress_thread() -> None:
    while not _STOP_EVENT.wait(TICK_SECONDS):
        state = _status()
        now = time.monotonic()
        elapsed = now - state.get("started_monotonic", now)
        expected_total = state.get("expected_total_seconds", 1.0)
        completed_estimate = state.get("completed_estimated_seconds", 0.0)
        current_estimate = state.get("current_estimated_seconds", 0.0)
        current_start = state.get("current_started_monotonic", now)
        current_elapsed = max(0.0, now - current_start)
        completed_actual = state.get("completed_actual_seconds", 0.0)
        completed_expected = max(completed_estimate, 1.0)
        scale = max(0.65, min(2.5, completed_actual / completed_expected))
        estimated_done = completed_estimate + min(current_elapsed, current_estimate) \
            + max(0.0, current_elapsed - current_estimate) * 0.15
        fraction = min(0.995, estimated_done / max(expected_total, 1.0))
        remaining = max(0.0, expected_total - completed_estimate - current_elapsed) * scale
        bar_width = 20
        filled = int(bar_width * fraction)
        bar = "#" * filled + "-" * (bar_width - filled)
        scope = state.get("scope", "preparing")
        idx = state.get("scope_index", 0)
        total = len(ELIGIBLE_SCOPES)
        fit_count = state.get("fit_count", 0)
        label = (f"{elapsed/60:.1f} min | {fraction*100:5.1f}% [{bar}] "
                 f"scope {idx}/{total} {scope} | fits {fit_count} | "
                 f"ETA ~{remaining/60:.1f} min")
        print(label, flush=True)
        with _STATE_LOCK:
            snapshot = dict(_STATE)
        snapshot.update({"status": "running", "elapsed_seconds": elapsed,
                         "approx_progress_fraction": fraction,
                         "approx_remaining_seconds": remaining,
                         "updated_unix_time": time.time()})
        _write_json(OUT / "RUN_PROGRESS.json", snapshot)


def _instrumented_checked_fit(*args, **kwargs):
    if time.monotonic() >= _DEADLINE:
        raise TimeoutError("Elapsed time passed the preregistered 2.5x runtime ceiling")
    fit, records = _ORIGINAL_CHECKED_FIT(*args, **kwargs)
    if time.monotonic() >= _DEADLINE:
        raise TimeoutError("Elapsed time crossed the 2.5x runtime ceiling after a completed fit")
    context = dict(args[5] if len(args) > 5 else kwargs["context"])
    context.setdefault("model", args[2] if len(args) > 2 else kwargs["model"])
    context.setdefault("lambda", args[3] if len(args) > 3 else kwargs["lam"])
    context.setdefault("rho", args[4] if len(args) > 4 else kwargs["rho"])
    with _STATE_LOCK:
        _STATE["fit_count"] = _STATE.get("fit_count", 0) + 1
    return instrument_predict(fit, context, _WARNING_LOG), records


def _prepare_scope(data, split, scope):
    train_match = match_ids(data, split["train"])
    if not train_match["passes_minimum_class_n"]:
        raise ValueError(f"Outer training matching falls below minimum for {scope}")
    train = data.take(_indices_for_pairs(data, train_match["pairs"]))
    train_sample_set = set(map(str, train.sample))
    if scope == "pooled/final":
        test_match = None
        test = None
    else:
        test_match = match_ids(data, split["test"])
        if not test_match["passes_minimum_class_n"]:
            raise ValueError(f"Outer test matching falls below minimum for {scope}")
        test = data.take(_indices_for_pairs(data, test_match["pairs"]))

    local = {str(sample): i for i, sample in enumerate(train.sample)}
    inner, expected_hashes, partition_counts = [], [], []
    for fold, frozen in enumerate(split["inner"]):
        restricted_train = [sample for sample in frozen["train"] if sample in train_sample_set]
        restricted_valid = [sample for sample in frozen["validation"] if sample in train_sample_set]
        a_match, b_match = match_ids(data, restricted_train), match_ids(data, restricted_valid)
        if not b_match["passes_minimum_class_n"]:
            raise ValueError(f"Inner validation {fold} falls below minimum for {scope}")
        ai = _indices_for_pairs(train, a_match["pairs"])
        bi = _indices_for_pairs(train, b_match["pairs"])
        a, b = train.take(ai), train.take(bi)
        a_donors, b_donors = set(map(str, a.donor)), set(map(str, b.donor))
        if a_donors & b_donors:
            raise AssertionError(f"Donor leakage in {scope} inner fold {fold}")
        a_ids, b_ids = set(map(str, a.sample)), set(map(str, b.sample))
        if a_ids & b_ids or not (a_ids | b_ids) <= train_sample_set:
            raise AssertionError(f"Matched membership differs from frozen S in {scope} fold {fold}")
        inner.append((np.asarray([local[str(s)] for s in a.sample], dtype=int),
                      np.asarray([local[str(s)] for s in b.sample], dtype=int)))
        expected_hashes.append(geo.array_hash(a.sample))
        partition_counts.append({
            "fold": fold,
            "train_restricted_samples": len(restricted_train),
            "train_pairs_after_rematching": a_match["matched_pairs"],
            "train_samples_after_rematching": a_match["matched_n"],
            "validation_restricted_samples": len(restricted_valid),
            "validation_pairs_after_rematching": b_match["matched_pairs"],
            "validation_samples_after_rematching": b_match["matched_n"],
            "validation_age_gap_years": b_match["age_gap_years"],
        })
    expected_hashes.append(geo.array_hash(train.sample))
    return train, inner, test, train_match, test_match, expected_hashes, partition_counts


def _plan(count_audit: dict) -> dict:
    rows = {row["scope"]: row for row in count_audit["scopes"]}
    if set(rows) != {row["scope"] for row in count_audit["scopes"]}:
        raise ValueError("Duplicate matching audit scope")
    eligible = [scope for scope in ELIGIBLE_SCOPES if rows.get(scope, {}).get("tuning_scope_evaluable")]
    if eligible != ELIGIBLE_SCOPES:
        raise ValueError(f"Frozen eligibility mismatch. Expected {ELIGIBLE_SCOPES}; audit says {eligible}")
    plans = []
    for scope in ELIGIBLE_SCOPES:
        pairs = int(rows[scope]["outer_train_S_sample_n"] // 2)
        seconds = PILOT_SECONDS * pairs / PILOT_PAIRS
        plans.append({"scope": scope, "matched_training_pairs_from_metadata_audit": pairs,
                      "estimated_seconds": seconds,
                      "reporting_role": "descriptive hyperparameter selection only; no PPMI evaluation"
                      if scope == "pooled/final" else
                      "matched-population held-out performance" if "outer/" not in scope
                      else "descriptive fold only; not an aggregate matched GENEPARK CV estimate"})
    est = sum(x["estimated_seconds"] for x in plans)
    return {"eligible_scopes": plans,
            "excluded_scopes": [row["scope"] for row in count_audit["scopes"]
                                if row["scope"] not in ELIGIBLE_SCOPES],
            "pilot_seconds": PILOT_SECONDS,
            "estimated_total_seconds": est,
            "estimated_total_minutes": est / 60,
            "runtime_stop_ceiling_seconds": est * STOP_MULTIPLE,
            "runtime_stop_ceiling_minutes": est * STOP_MULTIPLE / 60,
            "progress_update_seconds": TICK_SECONDS,
            "progress_note": "Approximate; runtime scales pilot seconds linearly by matched outer-training pairs.",
            "PPMI_loaded": False}


def _preflight(load_data: bool):
    release = json.loads(RELEASE.read_text(encoding="utf-8"))
    input_audit = json.loads(INPUT_AUDIT.read_text(encoding="utf-8"))
    count_audit = json.loads(MATCH_AUDIT.read_text(encoding="utf-8"))
    pilot = json.loads((PILOT_DIR / "RESULTS.json").read_text(encoding="utf-8"))
    if count_audit.get("status") != "metadata_count_audit_complete":
        raise ValueError("Metadata matching audit is incomplete")
    if count_audit.get("PPMI_loaded") is not False:
        raise ValueError("Matching audit indicates PPMI was loaded")
    if geo.digest(DECISION) != pilot.get("decision_sha256"):
        raise ValueError("Decision document has changed since the correctness pilot")
    if pilot.get("published_pipeline_sha256") != geo.digest(Path(geo.__file__)):
        raise ValueError("GEO pipeline differs from the pipeline used in the pilot")
    if pilot.get("input_sha256") != input_audit.get("input_sha256"):
        raise ValueError("Pilot prepared-input hashes differ from frozen input audit")
    if not pilot.get("preprocessing_fit_membership_matches"):
        raise ValueError("Pilot preprocessing membership/hash audit did not pass")
    plan = _plan(count_audit)
    if not load_data:
        return plan, None, None, None
    data, input_hashes = geo.load_geo()
    if input_hashes != input_audit["input_sha256"]:
        raise ValueError("Prepared GEO inputs differ from frozen input audit")
    splits = geo.split_manifest(data, input_hashes)
    split_hash = _split_hash(splits)
    if split_hash != release.get("split_manifest_sha256"):
        raise ValueError("Reconstructed frozen split manifest hash differs from release")
    split_by_scope = {entry["scope"]: entry for entry in splits["splits"]}
    if set(split_by_scope) != set(row["scope"] for row in count_audit["scopes"]):
        raise ValueError("Frozen split scope inventory differs from metadata audit")
    return plan, data, split_by_scope, {"input_sha256": input_hashes,
                                        "split_manifest_sha256": split_hash,
                                        "pipeline_sha256": geo.digest(Path(geo.__file__)),
                                        "decision_sha256": geo.digest(DECISION)}


def _run_scope(data, split, scope, scope_out, pilot_decisions=None):
    global _WARNING_LOG
    _PREP_FIT_SAMPLE_HASHES.clear()
    train, inner, test, train_match, test_match, expected_hashes, partition_counts = \
        _prepare_scope(data, split, scope)
    scope_out.mkdir(parents=True)
    _WARNING_LOG = scope_out / "prediction_warnings.jsonl"
    started, cpu_started = time.monotonic(), time.process_time()
    geo.checked_fit, geo.Prep = _instrumented_checked_fit, AuditedPrep
    try:
        prep, fits, decisions = geo.train_suite(train, inner, scope_out, scope)
    finally:
        geo.checked_fit, geo.Prep = _ORIGINAL_CHECKED_FIT, _ORIGINAL_PREP
    prep_matches = _PREP_FIT_SAMPLE_HASHES == expected_hashes
    if not prep_matches:
        raise AssertionError(f"Preprocessing train membership/hash mismatch for {scope}")
    if pilot_decisions is not None and decisions != pilot_decisions:
        raise AssertionError("Pilot-versus-integrated tuning decisions differ in GSE99039/outer/1")

    metrics, prediction_file = None, None
    if test is not None:
        predictions = geo.predict_suite(prep, fits, decisions, test)
        metrics = {name: float(geo.macro_auc(test, predictions[name]))
                   for name in geo.MODELS}
        prediction_file = "outer_test_predictions.npz"
        np.savez_compressed(scope_out / prediction_file,
                            sample=np.asarray(test.sample), y=np.asarray(test.y),
                            **{name: np.asarray(predictions[name]) for name in geo.MODELS})

    boundary = {}
    for model, dec in decisions.items():
        selected = dec.get("selected", dec.get("hyperparameters"))
        if selected is not None:
            boundary[model] = {
                "selected": selected,
                "persistent_boundary": bool(dec.get("persistent_boundary", False)),
                "expansion_rounds": dec.get("expansion_rounds"),
            }
    warning_count = (sum(bool(line.strip()) for line in _WARNING_LOG.read_text(encoding="utf-8").splitlines())
                     if _WARNING_LOG.exists() else 0)
    result = {
        "status": "scope_complete", "scope": scope,
        "seconds": time.monotonic() - started,
        "process_CPU_seconds": time.process_time() - cpu_started,
        "outer_train_matched_pairs": train_match["matched_pairs"],
        "outer_train_matched_samples": len(train.y),
        "outer_test_matched_pairs": None if test_match is None else test_match["matched_pairs"],
        "outer_test_scored": test is not None,
        "outer_test_metrics_matched_population_only": metrics,
        "outer_test_predictions_file": prediction_file,
        "inner_fold_count": len(inner), "inner_partitions": partition_counts,
        "preprocessing_fit_sample_hashes": list(_PREP_FIT_SAMPLE_HASHES),
        "expected_inner_train_plus_final_fit_hashes": expected_hashes,
        "preprocessing_fit_membership_matches": prep_matches,
        "pilot_tuning_decisions_identical": None if pilot_decisions is None else True,
        "selected_hyperparameters": {name: dec.get("selected", dec.get("hyperparameters"))
                                     for name, dec in decisions.items()},
        "boundary_selection_records": boundary,
        "prediction_warning_count": warning_count,
        "prediction_warning_log": str(_WARNING_LOG),
        "checked_fit_call_count": sum(1 for _ in (scope_out / "fits.jsonl").open(encoding="utf-8")),
        "pooled_final_used_for_PPMI_evaluation": False,
        "PPMI_loaded": False, "PPMI_scores_computed": False,
        "interpretation": "Held-out metrics, where present, describe performance in the matched population and are not directly comparable to the primary full-population estimate. Four eligible GSE99039 within-cohort folds remain descriptive individually; do not report their mean as matched GENEPARK CV.",
        "python": sys.version, "platform": platform.platform(),
    }
    _write_json(scope_out / "RESULTS.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", action="store_true", help="metadata/run integrity plan only; no model fits")
    args = parser.parse_args()
    # Even --plan reconstructs and verifies the frozen GEO inputs/splits, but it
    # never calls train_suite and therefore performs no model fits or scoring.
    plan, data, split_by_scope, integrity = _preflight(load_data=True)
    print(json.dumps({"status": "plan_ready" if args.plan else "preflight_passed",
                      "plan": plan, "integrity": integrity}, indent=2), flush=True)
    if args.plan:
        return 0
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite a prior full run: {OUT}")

    OUT.mkdir(parents=True)
    expected_total = float(plan["estimated_total_seconds"])
    started = time.monotonic()
    global _DEADLINE
    _DEADLINE = started + expected_total * STOP_MULTIPLE
    with _STATE_LOCK:
        _STATE.update({"started_monotonic": started,
                       "expected_total_seconds": expected_total,
                       "completed_estimated_seconds": 0.0,
                       "completed_actual_seconds": 0.0,
                       "scope_index": 0, "scope": "preparing",
                       "fit_count": 0})
    ticker = threading.Thread(target=_progress_thread, daemon=True)
    ticker.start()
    _write_json(OUT / "RUN_PLAN.json", {"plan": plan, "integrity": integrity,
                                        "stop_rule": "stop after a failed KKT retry, hash/fold-separation mismatch, exception/incomplete fit, or 2.5x runtime ceiling; boundary lambda selection is recorded, not a stop condition"})
    rows, boundaries = [], []
    pilot_decisions = json.loads((PILOT_DIR / "decisions.json").read_text(encoding="utf-8"))
    try:
        for idx, item in enumerate(plan["eligible_scopes"], start=1):
            scope = item["scope"]
            if time.monotonic() >= _DEADLINE:
                raise TimeoutError("Run crossed its 2.5x runtime ceiling between scopes")
            expected_scope = float(item["estimated_seconds"])
            with _STATE_LOCK:
                _STATE.update({"scope_index": idx, "scope": scope,
                               "current_started_monotonic": time.monotonic(),
                               "current_estimated_seconds": expected_scope,
                               "fit_count": 0})
            decisions = pilot_decisions if scope == PILOT_SCOPE else None
            result = _run_scope(data, split_by_scope[scope], scope,
                                OUT / scope.replace("/", "__"), decisions)
            rows.append(result)
            for model, record in result["boundary_selection_records"].items():
                if record["persistent_boundary"]:
                    boundaries.append({"scope": scope, "model": model, **record})
            _write_json(OUT / "BOUNDARY_SELECTIONS.json", {
                "scope": "this integrated age/sex-matching run only",
                "persistent_boundary_count": len(boundaries), "selections": boundaries})
            with _STATE_LOCK:
                _STATE["completed_estimated_seconds"] += expected_scope
                _STATE["completed_actual_seconds"] += result["seconds"]
                _STATE["scope_index"] = idx
                _STATE["scope"] = scope + " complete"
            print(f"COMPLETED {idx}/{len(ELIGIBLE_SCOPES)} {scope} "
                  f"({result['seconds']/60:.1f} min); warnings={result['prediction_warning_count']}; "
                  f"test_scored={result['outer_test_scored']}", flush=True)
        if time.monotonic() - started > expected_total * STOP_MULTIPLE:
            raise TimeoutError("Full run exceeded the 2.5x runtime ceiling")
        summary = {
            "status": "complete", "scope_count": len(rows), "scopes": rows,
            "persistent_boundary_selection_count": len(boundaries),
            "persistent_boundary_selections": boundaries,
            "pilot_decisions_reproduced_exactly": True,
            "four_GSE99039_inner_outer_folds_must_remain_descriptive_individually": True,
            "pooled_final_for_PPMI_evaluation": False,
            "PPMI_loaded": False,
            "elapsed_seconds": time.monotonic() - started,
            "plan": plan, "integrity": integrity,
        }
        _write_json(OUT / "SUMMARY.json", summary)
        print("FULL_MATCHING_RUN_COMPLETE", flush=True)
        print(json.dumps({"status": summary["status"],
                          "elapsed_minutes": summary["elapsed_seconds"] / 60,
                          "scope_count": summary["scope_count"],
                          "persistent_boundary_selection_count": len(boundaries),
                          "output": str(OUT / "SUMMARY.json")}, indent=2), flush=True)
        return 0
    except Exception as exc:
        failure = {"status": "failed_or_stopped", "error_type": type(exc).__name__,
                   "error": str(exc), "elapsed_seconds": time.monotonic() - started,
                   "completed_scopes": [row["scope"] for row in rows],
                   "traceback": traceback.format_exc(), "integrity": integrity}
        _write_json(OUT / "FAILURE.json", failure)
        print(f"STOPPED: {type(exc).__name__}: {exc}", flush=True)
        print(f"Diagnostic saved to {OUT / 'FAILURE.json'}", flush=True)
        return 2
    finally:
        _STOP_EVENT.set()
        ticker.join(timeout=2)


if __name__ == "__main__":
    raise SystemExit(main())
