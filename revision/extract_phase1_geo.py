"""Phase-1 read-only extraction from saved GEO artifacts (no fitting)."""

from __future__ import annotations

import csv
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "revision/outputs"
POOLED = ROOT / "private_data/geo_runs_20261018/pooled_final"
BASE_NAMES = ["age", "age_squared", "sex", "age_times_sex", "T_cells", "B_lineage",
              "Monocytic_lineage", "Neutrophils", "Erythroid", "Platelet"]


def write_csv(name: str, rows: list[dict], fields: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"Empty output: {name}")
    with (OUT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def model_specs() -> None:
    with np.load(POOLED / "preprocessing.npz", allow_pickle=False) as prep:
        genes = prep["genes"].astype(str)
        selected = prep["selected"].astype(int)
        keep = prep["keep"].astype(int)
    annotation = {}
    with (ROOT / "data/manifests/common_genes_gencode_v29.tsv").open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            annotation[row["gene_id"]] = row
    marker = json.loads((ROOT / "protocol/v3_marker_manifest.json").read_text())
    memberships = defaultdict(list)
    for panel, details in marker["cell_marker_panels"].items():
        for gene_id in details["gene_ids"]:
            memberships[gene_id].append(panel)
    sex_gate = {entry["gene_id"] for entry in marker["sex_markers"]}
    overlaps = {gene_id: panels for gene_id, panels in memberships.items() if len(panels) > 1}
    (OUT / "marker_union_audit.json").write_text(json.dumps({
        "panel_memberships": sum(len(x["gene_ids"]) for x in marker["cell_marker_panels"].values()),
        "unique_genes": len(memberships), "shared_gene_ids_and_panels": overlaps,
        "explanation": "173 panel memberships but 172 unique IDs when one gene belongs to two panels."
    }, indent=2) + "\n")

    decisions = json.loads((POOLED / "decisions.json").read_text())
    specs = []
    gene_rows = []
    for model in ("M1", "M2", "M3", "M4", "M5"):
        with np.load(POOLED / f"{model}_model.npz", allow_pickle=False) as data:
            coef = data["coef"].astype(float)
            intercept = float(data["intercept"])
        base_indices = [] if model == "M2" else ([i for i in keep if i < 4] if model in ("M1", "M3") else list(keep))
        n_base = len(base_indices)
        gene_coef = coef[n_base:]
        if len(gene_coef) not in (0, len(selected)):
            raise ValueError(f"Unexpected coefficient length in {model}")
        decision = decisions[model]
        selected_hyper = decision.get("selected", {})
        specs.append({
            "model": model, "lambda": selected_hyper.get("lambda", ""),
            "rho": selected_hyper.get("rho", ""),
            "grid_expansions": decision.get("expansion_rounds", ""),
            "no_gene_chosen": selected_hyper.get("baseline", False),
            "nonzero_genes": int(np.count_nonzero(gene_coef)),
            "intercept": intercept,
            "baseline_coefficients_standardized": json.dumps({BASE_NAMES[i]: float(coef[j]) for j, i in enumerate(base_indices)}),
            "source_file": str(POOLED / f"{model}_model.npz"),
        })
        if model in ("M2", "M3", "M5"):
            for j in np.flatnonzero(gene_coef):
                gene_id = str(genes[selected[j]])
                ref = annotation.get(gene_id, {})
                sym = ref.get("symbol", "")
                gene_rows.append({
                    "model": model, "gene_id": gene_id, "symbol": sym,
                    "coefficient_standardized": float(gene_coef[j]),
                    "panel_membership": ";".join(memberships.get(gene_id, [])),
                    "chromosome": ref.get("chromosome", ""),
                    "HB_prefix": sym.upper().startswith("HB"),
                    "erythroid_panel": "Erythroid" in memberships.get(gene_id, []),
                    "sex_gate_gene": gene_id in sex_gate,
                })
    write_csv("final_model_specs.csv", specs)
    write_csv("final_model_genes.csv", gene_rows)


def boundary_audit() -> None:
    table = (ROOT / "docs/GEO_SUPPLEMENTARY_BOUNDARY_TABLE_20261020_EN.md").read_text()
    rows = []
    in_detail = False
    for line in table.splitlines():
        if line.startswith("| Analysis | Scope or ComBat arm / scope"):
            in_detail = True
            continue
        if not in_detail or not line.startswith("|") or line.startswith("|---"):
            continue
        parts = [x.strip() for x in line.strip("|").split("|")]
        if len(parts) != 6:
            break
        analysis, scope, model, lam, rho, extensions = parts
        value = float(lam)
        rows.append({"analysis": analysis, "scope": scope, "model": model,
                     "lambda": value, "rho": float(rho), "grid_expansions": int(extensions),
                     "boundary_direction": "lower" if value < 1e-4 else "upper",
                     "source_document": "docs/GEO_SUPPLEMENTARY_BOUNDARY_TABLE_20261020_EN.md"})
    if len(rows) != 26:
        raise ValueError(f"Expected 26 documented boundary decisions; got {len(rows)}")
    write_csv("lambda_boundary_audit.csv", rows)


def batch_fold_auc() -> None:
    metadata = json.loads((ROOT / "prepared/GSE99039_sample_metadata.json").read_text())
    by_id = {r["gsm"]: r for r in metadata}
    rows = []
    for fold in range(5):
        path = ROOT / f"private_data/geo_runs_20261018/GSE99039_outer_{fold}/heldout_predictions.npz"
        with np.load(path, allow_pickle=False) as data:
            samples = data["samples"].astype(str)
            labels = data["y"].astype(int)
            scores = {m: data[m].astype(float) for m in ("M4", "M5")}
        batches = np.array([by_id[s]["batch"] for s in samples])
        if any(int(by_id[s]["case"]) != int(y) for s, y in zip(samples, labels)):
            raise ValueError(f"GEO fold {fold} prediction labels differ from metadata")
        for batch in np.unique(batches):
            ix = np.flatnonzero(batches == batch)
            n0 = int(np.sum(labels[ix] == 0))
            n1 = int(np.sum(labels[ix] == 1))
            eligible = n0 >= 3 and n1 >= 3
            row = {"fold": fold, "batch": batch, "n_control": n0, "n_PD": n1,
                   "n_pairs": n0 * n1, "eligible_ge3_each_class": eligible}
            for model in ("M4", "M5"):
                row[f"{model}_AUROC"] = float(roc_auc_score(labels[ix], scores[model][ix])) if eligible else ""
            rows.append(row)
    write_csv("within_batch_fold_cells.csv", rows)
    eligible = [r for r in rows if r["eligible_ge3_each_class"]]
    summary = {"definition": "Descriptive within-batch×outer-fold AUROC; eligible cells have >=3 per class. Weighted mean weights each cell by n_PD*n_control.",
               "eligible_cells": len(eligible), "total_cells": len(rows)}
    for model in ("M4", "M5"):
        values = np.array([r[f"{model}_AUROC"] for r in eligible], float)
        weights = np.array([r["n_pairs"] for r in eligible], float)
        summary[model] = {"unweighted": float(values.mean()) if len(values) else None,
                          "pair_weighted": float(np.average(values, weights=weights)) if len(values) else None}
    batch_rows = list(csv.DictReader((ROOT / "results/geo_within_batch_auc_20261020/batch_results.csv").open()))
    for omit in (False, True):
        selected = [r for r in batch_rows if r["evaluable"] == "True" and
                    (not omit or r["batch"] != "09/02/09_TRO")]
        summary["pooled_oof_without_09_02_09" if omit else "pooled_oof_all_eligible"] = {
            model: float(np.mean([float(r[f"{model}_within_batch_AUROC"]) for r in selected]))
            for model in ("M4", "M5")
        }
    (OUT / "within_batch_fold_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def probe_chip_proxy() -> None:
    """GENEPARK batch-only rule reapplied to recorded PROBE chip labels."""
    metadata = json.loads((ROOT / "prepared/GSE57475_sample_metadata.json").read_text())
    by_id = {row["gsm"]: row for row in metadata}
    manifest = json.loads((ROOT / "data/manifests/splits.json").read_text())
    rows = []
    for split in manifest["splits"]:
        if not split["scope"].startswith("GSE57475/outer/"):
            continue
        train = [by_id[x] for x in split["train"]]
        test = [by_id[x] for x in split["test"]]
        overall = float(np.mean([r["case"] for r in train]))
        rates = {batch: float(np.mean([r["case"] for r in train if r["batch"] == batch]))
                 for batch in {r["batch"] for r in train}}
        scores = [rates.get(r["batch"], overall) for r in test]
        rows.append({"scope": split["scope"], "n_train": len(train), "n_test": len(test),
                     "unseen_chip_fallback_count": sum(r["batch"] not in rates for r in test),
                     "chip_proxy_AUROC": float(roc_auc_score([r["case"] for r in test], scores))})
    if len(rows) != 5:
        raise ValueError("Expected 5 PROBE outer folds")
    write_csv("PROBE_chip_proxy_fold_AUROC.csv", rows)
    (OUT / "PROBE_chip_proxy_summary.json").write_text(json.dumps({
        "mean_fold_AUROC": float(np.mean([r["chip_proxy_AUROC"] for r in rows])),
        "n_folds": len(rows),
        "unseen_chip_fallback_total": sum(r["unseen_chip_fallback_count"] for r in rows),
        "label": "post hoc, GEO-only, metadata-only",
        "definition": "For each frozen outer fold, test samples receive the PD prevalence of their chip in the outer training subset; unseen chips use overall training PD prevalence. No smoothing."
    }, indent=2) + "\n")


def geo_characteristics() -> None:
    rows = []
    for source, label in (("GSE57475", "PROBE"), ("GSE99039", "GENEPARK-TRO")):
        metadata = json.loads((ROOT / f"prepared/{source}_sample_metadata.json").read_text())
        for diagnosis in (0, 1):
            group = [r for r in metadata if r["case"] == diagnosis]
            ages = [float(r["age"]) for r in group]
            rows.append({"source": label, "diagnosis": "PD" if diagnosis else "Control",
                         "n": len(group), "age_n": len(ages), "age_mean": statistics.mean(ages),
                         "age_sd": statistics.stdev(ages),
                         "female_n": sum(r["sex"] == 0 for r in group),
                         "male_n": sum(r["sex"] == 1 for r in group),
                         "disease_duration": "not available in prepared metadata",
                         "Hoehn_Yahr": "not available in prepared metadata",
                         "MDS_UPDRS_III": "not available in prepared metadata",
                         "dopaminergic_treatment": "not available in prepared metadata",
                         "LEDD": "not available in prepared metadata"})
    write_csv("GEO_candidate_characteristics_P1.csv", rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    model_specs()
    boundary_audit()
    batch_fold_auc()
    probe_chip_proxy()
    geo_characteristics()
    print("GEO_PHASE1_EXTRACTION_COMPLETE")


if __name__ == "__main__":
    main()
