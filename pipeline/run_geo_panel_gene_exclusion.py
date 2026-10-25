"""Run the registered 13-scope GEO panel-gene exclusion sensitivity.

Only M5 is tuned/refit. Frozen M4 inner scores and held-out predictions are
reused from the primary run; the locked primary pipeline is not modified.
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

OUT = ROOT / "private_data" / "geo_panel_gene_exclusion_full_20261019"
RELEASE = ROOT / "results/geo_development_20261018/GEO_CV_RELEASE_MANIFEST_20261018.json"
AUDIT = ROOT / "data/manifests/input_audit.json"
DECISION = ROOT / "docs/GEO_PANEL_GENE_EXCLUSION_DECISION_20261019_EN.md"
MARKERS = ROOT / "protocol/v3_marker_manifest.json"
ORIGINAL_RUN = ROOT / "private_data/geo_runs_20261018"
WARNING_LOG = OUT / "prediction_warnings.jsonl"
GRID = (.0001, .001, .01, .1, 1.)
RHOS = (0., .5, 1.)
EXPANSION_ROUNDS = 2
MAX_GENES = 1000
EXPECTED_MARKER_HASH = "729a9c450deadd6c65d0949c1e8f09ee76b455989a91580a60d41a1a758a8bc1"


class PanelExcludedPrep(geo.Prep):
    """Filter the fixed panel union before the frozen positive-MAD top-k."""

    def __init__(self, train, excluded_ids, max_genes=MAX_GENES):
        # Prep returns all eligible genes in the locked MAD/tie-break order.
        super().__init__(train, max_genes=len(train.genes))
        eligible = self.selected
        keep = np.array([str(self.genes[index]) not in excluded_ids for index in eligible])
        positions = np.flatnonzero(keep)[:max_genes]
        self.selected = eligible[positions]
        self.gene_mean = self.gene_mean[positions]
        self.gene_sd = self.gene_sd[positions]
        self.panel_gene_universe_count = int(sum(str(g) in excluded_ids for g in self.genes))
        self.panel_gene_eligible_count = int(sum(str(self.genes[i]) in excluded_ids
                                                 for i in eligible))
        if np.any(self.gene_sd <= 0):
            raise AssertionError("Panel-exclusion selected a constant gene")


def _split_scope(data, split):
    lookup = {sample: i for i, sample in enumerate(data.sample)}
    train = data.take(np.array([lookup[s] for s in split["train"]], dtype=int))
    heldout = data.take(np.array([lookup[s] for s in split["test"]], dtype=int))
    local = {sample: i for i, sample in enumerate(train.sample)}
    inner = [(np.array([local[s] for s in fold["train"]], dtype=int),
              np.array([local[s] for s in fold["validation"]], dtype=int))
             for fold in split["inner"]]
    return train, heldout, inner


def _original_m4(scope):
    original_dir = ORIGINAL_RUN / scope.replace("/", "_")
    decisions = json.loads((original_dir / "decisions.json").read_text(encoding="utf-8"))
    score = float(decisions["M4"]["inner_macro_AUROC"])
    return original_dir, score


def _tune_refit_m5(train, inner, excluded_ids, scope, scope_dir, m4_score):
    log_path = scope_dir / "fits.jsonl"
    if log_path.exists():
        raise ValueError(f"Refusing to append to existing fit log: {log_path}")
    fold_data = []
    for fold, (train_ix, valid_ix) in enumerate(inner):
        a, b = train.take(train_ix), train.take(valid_ix)
        if set(a.donor.tolist()) & set(b.donor.tolist()):
            raise AssertionError(f"Inner donor leakage in {scope}, fold {fold}")
        fold_data.append((PanelExcludedPrep(a, excluded_ids), a, b))
    final_prep = PanelExcludedPrep(train, excluded_ids)

    candidates = [dict(baseline=True, **{"lambda": 0.}, rho=0.,
                       score=float(m4_score), nonzero=False)]
    lambdas = list(GRID)
    for round_i in range(EXPANSION_ROUNDS + 1):
        existing = {(c["lambda"], c["rho"]) for c in candidates if not c["baseline"]}
        for lam in sorted(lambdas, reverse=True):
            for rho in RHOS:
                if (lam, rho) in existing:
                    continue
                scores, any_nonzero = [], False
                for fold, (prep, a, b) in enumerate(fold_data):
                    fit, records = geo.checked_fit(
                        prep, a, "M5", lam, rho,
                        dict(scope=scope, stage="inner", fold=fold,
                             baseline_only=False, expansion_round=round_i), log_path)
                    x, _ = prep.transform(b, "M5")
                    context = dict(scope=scope, stage="inner_validation", fold=fold,
                                   model="M5", **{"lambda": lam}, rho=rho)
                    wrapped = instrument_predict(fit, context, WARNING_LOG)
                    scores.append(geo.macro_auc(b, wrapped.predict(x)))
                    any_nonzero |= records[-1]["diagnostics"][
                        "nonzero_penalized_coefficients"] > 0
                candidates.append(dict(baseline=False, **{"lambda": float(lam)},
                    rho=float(rho), score=float(np.mean(scores)), nonzero=bool(any_nonzero)))
            geo.write_json(scope_dir / "M5_tuning.inprogress.json", candidates)
        selected = geo.choose(candidates)
        edge = (not selected["baseline"] and selected["nonzero"]
                and selected["lambda"] in (min(lambdas), max(lambdas)))
        if not edge or round_i == EXPANSION_ROUNDS:
            break
        lambdas.append(min(lambdas) / 10 if selected["lambda"] == min(lambdas)
                       else max(lambdas) * 10)

    result = dict(selected=dict(selected), candidate_count=len(candidates),
                  expansion_rounds=round_i, selected_gene_count=len(final_prep.selected),
                  panel_gene_universe_count=final_prep.panel_gene_universe_count,
                  panel_gene_eligible_count=final_prep.panel_gene_eligible_count,
                  final_gene_block_nonzero_coefficients=None)
    final_fit = None
    if selected["baseline"]:
        result["final_gene_block_status"] = "gene_block_not_selected_baseline_candidate"
    else:
        final_fit, records = geo.checked_fit(
            final_prep, train, "M5", selected["lambda"], selected["rho"],
            dict(scope=scope, stage="outer_or_final_refit", fold=None,
                 baseline_only=False), log_path)
        n_nonzero = int(records[-1]["diagnostics"]["nonzero_penalized_coefficients"])
        result["final_gene_block_nonzero_coefficients"] = n_nonzero
        if n_nonzero == 0:
            final_fit = None
            result["final_gene_block_status"] = "gene_block_not_selected_final_gene_block_zero"
        else:
            result["final_gene_block_status"] = "gene_block_selected"
            np.savez_compressed(scope_dir / "M5_model.npz",
                                intercept=final_fit.intercept_, coef=final_fit.coef_)
    fit_records = [json.loads(line) for line in log_path.read_text(encoding="utf-8").splitlines()]
    result["checked_fit_count"] = len(fit_records)
    result["fit_training_hashes"] = sorted({r["train_sample_sha256"] for r in fit_records})
    geo.write_json(scope_dir / "M5_decision.json", dict(selected=selected,
                    candidates=candidates, persistent_boundary=bool(edge),
                    expansion_rounds=round_i))
    return final_prep, final_fit, result


def _check_primary_inputs(release, audit):
    pipeline_hash = geo.digest(Path(geo.__file__))
    if pipeline_hash != release["pipeline_code_sha256"]:
        raise ValueError("Locked primary pipeline hash differs from the published release")
    data, hashes = geo.load_geo()
    if hashes != audit["input_sha256"]:
        raise ValueError("Prepared GEO inputs differ from the frozen input audit")
    splits = geo.split_manifest(data, hashes)
    split_hash = hashlib.sha256(json.dumps(
        splits, indent=2, ensure_ascii=False, allow_nan=False).encode()).hexdigest()
    if split_hash != release["split_manifest_sha256"]:
        raise ValueError("Reconstructed split manifest differs from the published release")
    return data, hashes, splits, pipeline_hash, split_hash


def main():
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite prior run: {OUT}")
    release = json.loads(RELEASE.read_text(encoding="utf-8"))
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    marker_manifest = json.loads(MARKERS.read_text(encoding="utf-8"))
    excluded_ids = {gene for panel in marker_manifest["cell_marker_panels"].values()
                    for gene in panel["gene_ids"]}
    marker_hash = hashlib.sha256(MARKERS.read_bytes()).hexdigest()
    if len(excluded_ids) != 172 or marker_hash != EXPECTED_MARKER_HASH:
        raise ValueError("Frozen six-panel gene union/hash does not match the decision")
    data, hashes, splits, pipeline_hash, split_hash = _check_primary_inputs(release, audit)
    all_splits = splits["splits"]
    if len(all_splits) != 13:
        raise ValueError(f"Expected 13 registered GEO scopes; found {len(all_splits)}")

    OUT.mkdir(parents=True)
    (OUT / "RUNNING").write_text("Full registered 13-scope sensitivity in progress\n",
                                  encoding="utf-8")
    total_started, total_cpu = time.monotonic(), time.process_time()
    results = {}
    for index, split in enumerate(all_splits, start=1):
        scope = split["scope"]
        scope_dir = OUT / scope.replace("/", "__")
        scope_dir.mkdir()
        train, heldout, inner = _split_scope(data, split)
        heldout_samples = set(heldout.sample.tolist())
        if heldout_samples & set(train.sample.tolist()):
            raise AssertionError(f"Outer heldout appears in training: {scope}")
        allowed_hashes = {geo.array_hash(train.sample)}
        for fold, (train_ix, valid_ix) in enumerate(inner):
            a, b = train.take(train_ix), train.take(valid_ix)
            if heldout_samples & set(a.sample.tolist()):
                raise AssertionError(f"Outer heldout appears in inner training: {scope}")
            if set(a.donor.tolist()) & set(b.donor.tolist()):
                raise AssertionError(f"Inner donor leakage in {scope}, fold {fold}")
            allowed_hashes.add(geo.array_hash(a.sample))

        # Verify the change leaves all frozen baseline features exactly intact.
        primary_prep = geo.Prep(train)
        excluded_prep = PanelExcludedPrep(train, excluded_ids)
        baseline_checks = {}
        for label, subset in (("train", train), ("heldout", heldout)):
            for model in ("M1", "M4"):
                baseline_checks[f"{label}_{model}_identical"] = bool(np.array_equal(
                    primary_prep.transform(subset, model)[0],
                    excluded_prep.transform(subset, model)[0]))
        if not all(baseline_checks.values()):
            raise AssertionError(f"M1/M4 transform changed in {scope}: {baseline_checks}")

        original_dir, m4_inner_score = _original_m4(scope)
        started, cpu_started = time.monotonic(), time.process_time()
        final_prep, final_fit, m5_result = _tune_refit_m5(
            train, inner, excluded_ids, scope, scope_dir, m4_inner_score)

        m4_probability = m5_probability = None
        heldout_metrics = None
        if len(heldout.y):
            primary = np.load(original_dir / "heldout_predictions.npz", allow_pickle=False)
            if not np.array_equal(primary["samples"], heldout.sample):
                raise AssertionError(f"Frozen M4 prediction sample order differs: {scope}")
            if not np.array_equal(primary["y"], heldout.y):
                raise AssertionError(f"Frozen heldout labels differ from split: {scope}")
            m4_probability = np.asarray(primary["M4"], dtype=float)
            if final_fit is None:
                # Protocol-mandated no-gene branch: reuse M4 exactly.
                m5_probability = m4_probability.copy()
            else:
                x, _ = final_prep.transform(heldout, "M5")
                context = dict(scope=scope, stage="heldout_scoring", fold=None,
                               model="M5", **{"lambda": m5_result["selected"]["lambda"]},
                               rho=m5_result["selected"]["rho"])
                wrapped = instrument_predict(final_fit, context, WARNING_LOG)
                m5_probability = wrapped.predict(x)
            heldout_metrics = dict(
                M4=geo.metrics(heldout.y, m4_probability),
                M5_panel_gene_excluded=geo.metrics(heldout.y, m5_probability))
            heldout_metrics["M5_minus_M4"] = {
                key: float(heldout_metrics["M5_panel_gene_excluded"][key]
                           - heldout_metrics["M4"][key])
                for key in ("AUROC", "AUPRC", "Brier", "log_loss")}
            np.savez_compressed(scope_dir / "heldout_predictions.npz",
                                samples=heldout.sample, y=heldout.y,
                                M4=m4_probability, M5_panel_gene_excluded=m5_probability)

        fit_hashes = set(m5_result["fit_training_hashes"])
        if not fit_hashes <= allowed_hashes:
            raise AssertionError(f"M5 fit used non-frozen training partition in {scope}")
        scope_result = dict(
            status="complete", scope=scope,
            seconds=time.monotonic() - started,
            process_CPU_seconds=time.process_time() - cpu_started,
            inner_fold_count=len(inner), baseline_transform_exactly_unchanged=baseline_checks,
            M4_inner_macro_AUROC_reused=m4_inner_score,
            M1_to_M4_refit=False, M5_only_refit=True,
            selected_gene_count=m5_result["selected_gene_count"],
            panel_gene_union_count=len(excluded_ids),
            panel_gene_universe_count=m5_result["panel_gene_universe_count"],
            panel_gene_eligible_count=m5_result["panel_gene_eligible_count"],
            M5=m5_result, heldout_metrics=heldout_metrics,
            heldout_samples_scored=bool(len(heldout.y)),
            checked_fit_training_hashes_match_frozen_partitions=True,
            PPMI_loaded=False, PPMI_scores_computed=False)
        geo.write_json(scope_dir / "RESULTS.json", scope_result)
        results[scope] = scope_result
        geo.write_json(OUT / "progress.json", dict(completed_scope_count=len(results),
                         expected_scope_count=len(all_splits), completed_scopes=list(results)))
        print(f"[{index}/13] {scope} complete; M5 checked fits: "
              f"{m5_result['checked_fit_count']}", flush=True)

    warning_count = sum(1 for _ in WARNING_LOG.open(encoding="utf-8")) if WARNING_LOG.exists() else 0
    result = dict(
        status="complete", sensitivity="six_panel_gene_exclusion_from_M5_added_block",
        scope_count=len(results), scopes=results,
        seconds=time.monotonic() - total_started,
        process_CPU_seconds=time.process_time() - total_cpu,
        checked_fit_count=sum(s["M5"]["checked_fit_count"] for s in results.values()),
        all_baseline_transforms_exactly_unchanged=all(
            all(s["baseline_transform_exactly_unchanged"].values()) for s in results.values()),
        all_fit_training_hashes_match_frozen_partitions=all(
            s["checked_fit_training_hashes_match_frozen_partitions"] for s in results.values()),
        marker_panel_union_count=len(excluded_ids), marker_manifest_sha256=marker_hash,
        original_pipeline_sha256=pipeline_hash,
        published_split_manifest_sha256=split_hash, input_sha256=hashes,
        decision_sha256=geo.digest(DECISION), prediction_warning_log=str(WARNING_LOG),
        prediction_warning_count=warning_count,
        python=sys.version, platform=platform.platform(),
        PPMI_loaded=False, PPMI_scores_computed=False,
        interpretation="Registered GEO sensitivity; per-scope heldout scores are descriptive, not an independent pooled estimate.")
    geo.write_json(OUT / "RESULTS.json", result)
    (OUT / "RUNNING").unlink()
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
