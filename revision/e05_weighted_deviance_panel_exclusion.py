#!/usr/bin/env python3
"""E5 V1: weighted-deviance tuning plus panel-gene exclusion, staged transfer."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
import geo_pipeline as geo
from e03_ppmi_nested_cv import fit_pair, score_model, donor_folds
from e04_transductive_universe import (EligiblePrep, metrics, ppmi_data,
                                       restore_prediction, save_fit)
from postlock_common import (private_output, require_execution_authorization,
                             sha256, synthetic_ranks, write_aggregate)


def panel_excluded_indices(data: geo.Data) -> np.ndarray:
    excluded = np.unique(np.concatenate(data.marker_indices))
    eligible = np.setdiff1d(np.arange(len(data.genes)), excluded, assume_unique=True)
    if len(eligible) == 0:
        raise ValueError("No gene remains after marker-panel exclusion")
    return eligible


def code_hashes() -> dict:
    folder = Path(__file__).parent
    return {name: sha256(folder / name) for name in (
        "e05_weighted_deviance_panel_exclusion.py", "e04_transductive_universe.py",
        "e03_ppmi_nested_cv.py")}


def geo_stage(output: Path) -> None:
    data, input_hashes = geo.load_geo()
    eligible = panel_excluded_indices(data)
    np.save(output / "eligible_indices.npy", eligible, allow_pickle=False)
    report = {"analysis": "E5", "stage": "GEO_complete_before_PPMI_scoring",
              "exploratory_two_factor_variant": True,
              "tuning_objective": "negative_train_weighted_validation_log_loss",
              "excluded_panel_gene_count": int(len(data.genes) - len(eligible)),
              "eligible_gene_ids_sha256": geo.array_hash(data.genes[eligible]),
              "code_sha256": code_hashes(), "geo_input_sha256": input_hashes,
              "transfer": {}}
    factory = lambda train, max_genes: EligiblePrep(train, eligible, max_genes)
    for source in sorted(set(map(str, data.source))):
        train = data.take(np.flatnonzero(data.source == source))
        test = data.take(np.flatnonzero(data.source != source))
        inner = geo.donor_folds(train, 5, f"E5/{source}/transfer/inner")
        out = output / "transfer" / source
        prep, m4, m5, decision = fit_pair(train, inner, out,
            f"E5/{source}/transfer", prep_factory=factory,
            objective="weighted_deviance")
        save_fit(out, prep, m4, m5, decision)
        p4 = score_model(prep, m4, test, "M4")
        p5 = p4 if m5 is None else score_model(prep, m5, test, "M5")
        report["transfer"][f"{source}_to_{str(test.source[0])}"] = {
            "M4": metrics(test.y, p4), "M5": metrics(test.y, p5),
            "M5_minus_M4_AUROC": float(roc_auc_score(test.y, p5)-roc_auc_score(test.y, p4)),
            "M5_minus_M4_AP": float(metrics(test.y, p5)["AP"]-metrics(test.y, p4)["AP"]),
            "M5_minus_M4_log_loss": float(metrics(test.y, p5)["log_loss"]-metrics(test.y, p4)["log_loss"]),
            "selection": decision}
        print(f"E5 GEO transfer {source} complete", flush=True)
    inner = geo.donor_folds(data, 5, "E5/pooled/final/inner")
    out = output / "pooled"
    prep, m4, m5, decision = fit_pair(data, inner, out, "E5/pooled/final",
        prep_factory=factory, objective="weighted_deviance")
    save_fit(out, prep, m4, m5, decision)
    report["pooled_selection"] = decision
    write_aggregate(output / "E5_GEO_AGGREGATE.json", report)
    print("E5 GEO stage frozen; PPMI scoring is a separate command", flush=True)


def ppmi_stage(output: Path, bundle: Path, binding: Path) -> None:
    report = json.loads((output / "E5_GEO_AGGREGATE.json").read_text(encoding="utf-8"))
    if (report.get("stage") != "GEO_complete_before_PPMI_scoring"
            or report.get("code_sha256") != code_hashes()):
        raise ValueError("E5 GEO stage or code changed before PPMI scoring")
    target = ppmi_data(bundle, binding)
    p4 = restore_prediction(output / "pooled", target, "M4")
    p5 = restore_prediction(output / "pooled", target, "M5")
    m4, m5 = metrics(target.y, p4), metrics(target.y, p5)
    result = {"analysis": "E5", "stage": "PPMI_after_frozen_GEO",
              "exploratory_two_factor_variant": True, "M4": m4, "M5": m5,
              "M5_minus_M4": {key: m5[key] - m4[key]
                              for key in ("AUROC", "AP", "log_loss")}}
    write_aggregate(output / "E5_PPMI_AGGREGATE.json", result)
    print("E5_COMPLETE: private aggregate PPMI output written")


def synthetic_self_test() -> None:
    rng = np.random.default_rng(51)
    n, g = 48, 40
    ranks = synthetic_ranks(n, g, 51)
    data = geo.Data(ranks, rng.uniform(50, 75, n), rng.integers(0, 2, n),
                    np.array([0, 1] * (n//2)), np.repeat("SYNTHETIC", n),
                    np.array([f"s{i}" for i in range(n)]),
                    np.array([f"s{i}" for i in range(n)]),
                    np.array([f"G{i:04d}" for i in range(g)]),
                    [np.arange(i*3, i*3+3) for i in range(6)])
    eligible = panel_excluded_indices(data)
    assert len(eligible) == g - 18
    import tempfile
    with tempfile.TemporaryDirectory() as temporary:
        out = Path(temporary) / "fit"
        prep, m4, m5, decision = fit_pair(data, list(donor_folds(data, 2, 900)),
            out, "E5/synthetic", max_genes=8, grid=(.1,), rhos=(0.,),
            expansion_rounds=0,
            prep_factory=lambda train, max_genes: EligiblePrep(train, eligible, max_genes),
            objective="weighted_deviance")
        assert set(prep.selected).isdisjoint(set(np.concatenate(data.marker_indices)))
        save_fit(out, prep, m4, m5, decision)
        direct = score_model(prep, m4, data, "M4")
        assert np.allclose(direct, restore_prediction(out, data, "M4"), atol=1e-10)
    print("E5_SYNTHETIC_TEST_PASSED")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic-self-test", action="store_true")
    parser.add_argument("--stage", choices=("geo", "ppmi"))
    parser.add_argument("--authorization-record", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--binding", type=Path)
    args = parser.parse_args()
    if args.synthetic_self_test:
        synthetic_self_test()
        return
    if any(value is None for value in (args.stage, args.authorization_record, args.output_dir)):
        parser.error("Real E5 requires stage, authorization and output")
    require_execution_authorization(args.authorization_record, analysis="E5")
    os.umask(0o077)
    output = private_output(args.output_dir)
    if args.stage == "geo":
        if output.exists():
            raise ValueError("Use a fresh E5 GEO output directory")
        output.mkdir(mode=0o700, parents=True)
        geo_stage(output)
    else:
        if args.bundle is None or args.binding is None:
            parser.error("PPMI stage requires bundle and binding")
        ppmi_stage(output, args.bundle, args.binding)


if __name__ == "__main__":
    main()
