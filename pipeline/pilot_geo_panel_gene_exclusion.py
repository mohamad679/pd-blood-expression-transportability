"""One-scope correctness/runtime pilot for the registered M5 panel-gene exclusion."""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import geo_pipeline as geo
from pipeline.prediction_warning_capture import instrument_predict

OUT = ROOT / "private_data" / "geo_panel_gene_exclusion_pilot_20261019"
RELEASE = ROOT / "results/geo_development_20261018/GEO_CV_RELEASE_MANIFEST_20261018.json"
AUDIT = ROOT / "data/manifests/input_audit.json"
DECISION = ROOT / "docs/GEO_PANEL_GENE_EXCLUSION_DECISION_20261019_EN.md"
MARKERS = ROOT / "protocol/v3_marker_manifest.json"
ORIGINAL_RUN = ROOT / "private_data/geo_runs_20261018"
PILOT_SCOPE = "GSE57475/outer/0"
GRID = (.0001, .001, .01, .1, 1.)
RHOS = (0., .5, 1.)
EXPANSION_ROUNDS = 2
MAX_GENES = 1000


class PanelExcludedPrep(geo.Prep):
    """Apply the frozen Prep, then filter panel IDs before its MAD top-k."""

    def __init__(self, train, excluded_ids, max_genes=MAX_GENES):
        # Ask the locked Prep for the complete eligible positive-MAD ordering;
        # filtering this ordering is equivalent to excluding before ranking.
        super().__init__(train, max_genes=len(train.genes))
        eligible = self.selected
        allowed = np.array([str(self.genes[index]) not in excluded_ids
                            for index in eligible], dtype=bool)
        kept_positions = np.flatnonzero(allowed)[:max_genes]
        self.selected = eligible[kept_positions]
        self.gene_mean = self.gene_mean[kept_positions]
        self.gene_sd = self.gene_sd[kept_positions]
        self.panel_gene_universe_count = int(sum(
            str(gene) in excluded_ids for gene in self.genes))
        self.panel_gene_eligible_count = int(sum(
            str(self.genes[index]) in excluded_ids for index in eligible))
        if np.any(self.gene_sd <= 0):
            raise AssertionError("Panel-exclusion selected a constant gene")


def _split_scope(data, split):
    lookup = {sample: i for i, sample in enumerate(data.sample)}
    train = data.take(np.array([lookup[s] for s in split["train"]], dtype=int))
    local = {sample: i for i, sample in enumerate(train.sample)}
    inner = [(np.array([local[s] for s in fold["train"]], dtype=int),
              np.array([local[s] for s in fold["validation"]], dtype=int))
             for fold in split["inner"]]
    return train, inner


def _m5_only_tune(train, inner, prep_factory, scope, out, baseline_m4_score):
    """Mirror only the locked M5 tuning/refit branch; never fit M1–M4."""
    log_path = out / "fits.jsonl"
    if log_path.exists():
        raise ValueError(f"Refusing to append to existing fit log: {log_path}")
    fold_data = []
    for fold, (train_ix, valid_ix) in enumerate(inner):
        a, b = train.take(train_ix), train.take(valid_ix)
        if set(a.donor.tolist()) & set(b.donor.tolist()):
            raise AssertionError(f"Inner donor leakage in fold {fold}")
        fold_data.append((prep_factory(a), a, b))
    final_prep = prep_factory(train)
    candidates = [dict(baseline=True, **{"lambda": 0.}, rho=0.,
                       score=float(baseline_m4_score), nonzero=False)]
    lambdas = list(GRID)
    for round_i in range(EXPANSION_ROUNDS + 1):
        existing = {(candidate["lambda"], candidate["rho"])
                    for candidate in candidates if not candidate["baseline"]}
        for lam in sorted(lambdas, reverse=True):
            for rho in RHOS:
                if (lam, rho) in existing:
                    continue
                scores, any_nonzero = [], False
                for fold, (prep, a, b) in enumerate(fold_data):
                    fit, records = geo.checked_fit(
                        prep, a, "M5", lam, rho,
                        dict(scope=scope, stage="inner", fold=fold,
                             baseline_only=False, expansion_round=round_i),
                        log_path)
                    validation_x, _ = prep.transform(b, "M5")
                    context = dict(scope=scope, stage="inner_validation",
                                   fold=fold, model="M5", **{"lambda": lam}, rho=rho)
                    wrapped = instrument_predict(fit, context,
                                                 out / "prediction_warnings.jsonl")
                    scores.append(geo.macro_auc(b, wrapped.predict(validation_x)))
                    any_nonzero |= (records[-1]["diagnostics"]
                                    ["nonzero_penalized_coefficients"] > 0)
                candidates.append(dict(baseline=False, **{"lambda": float(lam)},
                    rho=float(rho), score=float(np.mean(scores)),
                    nonzero=bool(any_nonzero)))
        selected = geo.choose(candidates)
        edge = (not selected["baseline"] and selected["nonzero"]
                and selected["lambda"] in (min(lambdas), max(lambdas)))
        if not edge or round_i == EXPANSION_ROUNDS:
            break
        lambdas.append(min(lambdas) / 10 if selected["lambda"] == min(lambdas)
                       else max(lambdas) * 10)

    outcome = dict(selected=dict(selected), selected_gene_count=int(len(final_prep.selected)),
                   panel_gene_universe_count=final_prep.panel_gene_universe_count,
                   panel_gene_eligible_count=final_prep.panel_gene_eligible_count,
                   final_gene_block_nonzero_coefficients=None,
                   final_gene_block_status="not_refit_in_pilot")
    if selected["baseline"]:
        outcome["final_gene_block_status"] = "gene_block_not_selected_baseline_candidate"
    else:
        fit, records = geo.checked_fit(
            final_prep, train, "M5", selected["lambda"], selected["rho"],
            dict(scope=scope, stage="outer_or_final_refit", fold=None,
                 baseline_only=False), log_path)
        n_nonzero = int(records[-1]["diagnostics"]["nonzero_penalized_coefficients"])
        outcome["final_gene_block_nonzero_coefficients"] = n_nonzero
        outcome["final_gene_block_status"] = (
            "gene_block_not_selected_final_fit_zero" if n_nonzero == 0
            else "gene_block_selected")
    outcome["candidate_count"] = len(candidates)
    outcome["expansion_rounds"] = round_i
    outcome["fit_attempt_records"] = sum(
        1 for _ in log_path.read_text(encoding="utf-8").splitlines())
    return outcome


def main():
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite prior pilot: {OUT}")
    release = json.loads(RELEASE.read_text(encoding="utf-8"))
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    marker_manifest = json.loads(MARKERS.read_text(encoding="utf-8"))
    excluded_ids = {gene for panel in marker_manifest["cell_marker_panels"].values()
                    for gene in panel["gene_ids"]}
    if len(excluded_ids) != 172:
        raise ValueError(f"Expected fixed 172-gene panel union, found {len(excluded_ids)}")
    marker_hash = hashlib.sha256(MARKERS.read_bytes()).hexdigest()
    expected_marker_hash = "729a9c450deadd6c65d0949c1e8f09ee76b455989a91580a60d41a1a758a8bc1"
    if marker_hash != expected_marker_hash:
        raise ValueError("Marker panel manifest differs from the decision")
    pipeline_hash = geo.digest(Path(geo.__file__))
    if pipeline_hash != release["pipeline_code_sha256"]:
        raise ValueError("Locked primary pipeline hash differs from published release")
    data, input_hashes = geo.load_geo()
    if input_hashes != audit["input_sha256"]:
        raise ValueError("Prepared GEO inputs differ from frozen input audit")
    splits = geo.split_manifest(data, input_hashes)
    split_hash = hashlib.sha256(json.dumps(
        splits, indent=2, ensure_ascii=False, allow_nan=False).encode()).hexdigest()
    if split_hash != release["split_manifest_sha256"]:
        raise ValueError("Reconstructed split manifest differs from published release")
    split = next(item for item in splits["splits"] if item["scope"] == PILOT_SCOPE)
    original_dir = ORIGINAL_RUN / "GSE57475_outer_0"
    original_decisions = json.loads((original_dir / "decisions.json").read_text())
    m4_score = original_decisions["M4"]["inner_macro_AUROC"]
    train, inner = _split_scope(data, split)
    for fold, (train_ix, valid_ix) in enumerate(inner):
        a, b = train.take(train_ix), train.take(valid_ix)
        if set(a.donor.tolist()) & set(b.donor.tolist()):
            raise AssertionError(f"Inner donor leakage in fold {fold}")
        primary_prep = geo.Prep(a)
        panel_prep = PanelExcludedPrep(a, excluded_ids)
        for subset in (a, b):
            if not np.array_equal(primary_prep.transform(subset, "M4")[0],
                                  panel_prep.transform(subset, "M4")[0]):
                raise AssertionError(f"M4 transform changed in inner fold {fold}")

    OUT.mkdir(parents=True)
    started, cpu_started = time.monotonic(), time.process_time()
    outcome = _m5_only_tune(
        train, inner, lambda partition: PanelExcludedPrep(partition, excluded_ids),
        PILOT_SCOPE, OUT, m4_score)
    result = dict(
        status="timing_pilot_complete", sensitivity="six_panel_gene_exclusion_from_M5_added_block",
        scope=PILOT_SCOPE, seconds=time.monotonic()-started,
        process_CPU_seconds=time.process_time()-cpu_started,
        inner_fold_count=len(inner), M5_only_refit=True, M1_to_M4_refit=False,
        M4_baseline_score_reused_from_original_run=True,
        marker_panel_union_count=len(excluded_ids), marker_manifest_sha256=marker_hash,
        selected_panel_gene_count_after_refit=outcome["selected_gene_count"],
        panel_gene_universe_count=outcome["panel_gene_universe_count"],
        panel_gene_eligible_count=outcome["panel_gene_eligible_count"],
        result=outcome, published_pipeline_sha256=pipeline_hash,
        published_split_manifest_sha256=split_hash, input_sha256=input_hashes,
        decision_sha256=geo.digest(DECISION),
        outer_holdout_loaded=False, heldout_performance_scored=False,
        PPMI_loaded=False, PPMI_scores_computed=False,
        metrics_reported=False,
        interpretation="Single-scope correctness/runtime pilot only; not an effect estimate.")
    geo.write_json(OUT / "RESULTS.json", result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
