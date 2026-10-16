"""Prepare public GEO development inputs without fitting any model."""
import collections
import csv
import gzip
import hashlib
import json
import math
import re
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "prepared"


def rank_percentiles(matrix):
    result = np.empty_like(matrix, dtype=np.float64)
    for index, row in enumerate(matrix):
        _, inverse, counts = np.unique(row, return_inverse=True, return_counts=True)
        midranks = np.cumsum(counts) - (counts - 1) / 2
        result[index] = (midranks[inverse] - 0.5) / row.size
    return result


def metadata(accession):
    result, current = [], None
    for line in (ROOT / "metadata" / f"{accession}.soft").read_text().splitlines():
        if line.startswith("^SAMPLE ="):
            if current is not None:
                result.append(current)
            current = {"gsm": line.split("=", 1)[1].strip()}
        elif current is not None and line.startswith("!Sample_title ="):
            current["title"] = line.split("=", 1)[1].strip()
        elif current is not None and line.startswith("!Sample_characteristics_ch1 ="):
            value = line.split("=", 1)[1].strip()
            if ":" in value:
                key, value = value.split(":", 1)
                current[key.strip().lower()] = value.strip()
    if current is not None:
        result.append(current)
    return {row["gsm"]: row for row in result}


def eligible(accession, rows):
    selected = []
    for row in rows.values():
        label = row.get("disease state") if accession == "GSE57475" else row.get("disease label")
        labels = ({"Dopamine Transporter Imaging (DAT)-confirmed PD": 1, "healthy control": 0}
                  if accession == "GSE57475" else {"IPD": 1, "CONTROL": 0})
        if label not in labels:
            continue
        age_string = row.get("age") if accession == "GSE57475" else row.get("age at exam")
        try:
            age = float(age_string)
            if not math.isfinite(age) or not 0 < age < 120:
                continue
        except (ValueError, TypeError):
            continue
        sex = row.get("gender", "") if accession == "GSE57475" else row.get("sex", "")
        coding = {"M": 1, "F": 0, "Male": 1, "Female": 0}
        if sex not in coding:
            continue
        batch, kind, position = row.get("batch", ""), "provided_batch", ""
        if accession == "GSE57475":
            match = re.fullmatch(r"(\d+)_([A-Z])", row.get("title", ""))
            if match is None:
                raise ValueError("Unrecognized Illumina title; no invented chip proxy")
            batch, position, kind = match.group(1), match.group(2), "title_derived_chip_proxy"
        selected.append({"gsm": row["gsm"], "case": labels[label], "age": age, "sex": coding[sex],
                         "subject_id": row.get("subject_id", ""), "batch": batch,
                         "batch_kind": kind, "array_position": position})
    return selected


def probe_expression(accession):
    path = ROOT / "expression" / f"{accession}_series_matrix.txt.gz"
    samples, probes, values = [], [], []
    inside, header = False, False
    with gzip.open(path, "rt", newline="") as handle:
        for line in handle:
            if line.startswith("!Sample_geo_accession"):
                samples = next(csv.reader([line], delimiter="\t"))[1:]
            if line.startswith("!series_matrix_table_begin"):
                inside, header = True, True
                continue
            if line.startswith("!series_matrix_table_end"):
                inside = False
            if inside:
                cells = next(csv.reader([line], delimiter="\t"))
                if header:
                    if cells[1:] != samples:
                        raise ValueError("Matrix sample header mismatch")
                    header = False
                    continue
                probes.append(cells[0])
                values.append(np.asarray(cells[1:], dtype=np.float64))
    expression = np.stack(values)
    if not np.isfinite(expression).all():
        raise ValueError("Nonfinite source values; investigate before preparing features")
    return samples, probes, expression


def main():
    OUT.mkdir(exist_ok=True)
    with (ROOT / "mapping" / "common_genes_gencode_v29.tsv").open(newline="") as handle:
        gene_rows = list(csv.DictReader(handle, delimiter="\t"))
    genes = [row["gene_id"] for row in gene_rows]
    if genes != sorted(set(genes)):
        raise ValueError("Common gene list not unique and sorted")
    audit = {"transformation": "Median of uniquely mapped probes per gene, then per-sample midrank percentile (rank-0.5)/G",
             "universe": "common_genes_gencode_v29.tsv", "genes": len(genes),
             "no_model_fitting": True, "no_cross_sample_normalization_added": True,
             "inherited_source_processing": {"GSE57475": "GenomeStudio average normalization", "GSE99039": "RMA/log2/quantile normalization"},
             "datasets": {}}
    for accession, platform in [("GSE57475", "GPL6947"), ("GSE99039", "GPL570")]:
        selected = eligible(accession, metadata(accession))
        samples, probes, expression = probe_expression(accession)
        by_sample = {sample: i for i, sample in enumerate(samples)}
        by_probe = {probe: i for i, probe in enumerate(probes)}
        if len(by_sample) != len(samples) or len(by_probe) != len(probes):
            raise ValueError("Duplicate source keys")
        indices = [by_sample[row["gsm"]] for row in selected]
        gene_probes = collections.defaultdict(list)
        with (ROOT / "mapping" / f"{platform}_probe_to_gene.tsv").open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                gene_probes[row["gene_id"]].append(by_probe[row["probe_id"]])
        gene_expression = np.column_stack([
            np.median(expression[gene_probes[gene]][:, indices], axis=0) for gene in genes
        ])
        ranks = rank_percentiles(gene_expression)
        if not np.isfinite(gene_expression).all() or not np.isfinite(ranks).all():
            raise ValueError("Nonfinite prepared values")
        if not np.allclose(ranks.mean(axis=1), 0.5):
            raise ValueError("Per-sample rank invariant failed")
        output = OUT / f"{accession}_development.npz"
        np.savez_compressed(output, gene_ids=np.asarray(genes), sample_ids=np.asarray([row["gsm"] for row in selected]),
                            gene_expression=gene_expression.astype(np.float32), rank_features=ranks.astype(np.float32),
                            age=np.asarray([row["age"] for row in selected]), sex=np.asarray([row["sex"] for row in selected]),
                            y=np.asarray([row["case"] for row in selected]), chromosomes=np.asarray([row["chromosome"] for row in gene_rows]))
        (OUT / f"{accession}_sample_metadata.json").write_text(json.dumps(selected, indent=2) + "\n")
        reloaded = np.load(output, allow_pickle=False)
        if not np.array_equal(reloaded["gene_ids"], np.asarray(genes)) or reloaded["rank_features"].shape != (len(selected), len(genes)):
            raise ValueError("Saved input alignment failed")
        audit["datasets"][accession] = {
            "source_samples": len(samples), "source_probes": len(probes), "source_nonfinite_values": 0,
            "selected_samples": len(selected), "cases": sum(row["case"] for row in selected),
            "controls": sum(not row["case"] for row in selected), "female": sum(not row["sex"] for row in selected),
            "male": sum(row["sex"] for row in selected), "genes": len(genes),
            "source_value_min": float(expression.min()), "source_value_max": float(expression.max()),
            "range_of_sample_medians": [float(np.median(expression, axis=0).min()), float(np.median(expression, axis=0).max())],
            "prepared_bytes": output.stat().st_size, "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
            "explicit_unique_subject_ids": len({row["subject_id"] for row in selected if row["subject_id"]}),
            "batch_suffix_counts": dict(collections.Counter(row["batch"].rsplit("_", 1)[-1] for row in selected)) if accession == "GSE99039" else {},
            "distinct_batches_or_chip_proxies": len({row["batch"] for row in selected}),
            "batch_kind": selected[0]["batch_kind"],
        }
        del expression
    (OUT / "preprocessing_audit_summary.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
