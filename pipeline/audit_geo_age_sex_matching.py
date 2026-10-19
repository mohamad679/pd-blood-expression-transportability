#!/usr/bin/env python3
"""Metadata-only, fold-safe counts for the registered GEO matching sensitivity."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import csv

import numpy as np

from pipeline import geo_pipeline as geo


ROOT = geo.ROOT
CALIPER_YEARS = 3.0
MIN_CLASS_N = 5


class SplitData:
    """Minimal metadata view required by the frozen GEO split builder."""

    def __init__(self, age, sex, y, source, donor, sample):
        self.age, self.sex, self.y = age, sex, y
        self.source, self.donor, self.sample = source, donor, sample

    def take(self, indices):
        ix = np.asarray(indices, dtype=int)
        return SplitData(*(getattr(self, key)[ix] for key in
                           ("age", "sex", "y", "source", "donor", "sample")))


def metadata_only():
    parts = []
    hashes = {}
    for accession in ("GSE57475", "GSE99039"):
        npz_path = ROOT / "prepared" / f"{accession}_development.npz"
        meta_path = ROOT / "prepared" / f"{accession}_sample_metadata.json"
        hashes[npz_path.name] = geo.digest(npz_path)
        hashes[meta_path.name] = geo.digest(meta_path)
        with np.load(npz_path, allow_pickle=False) as archive:
            sample_ids = archive["sample_ids"].astype(str)
            age = archive["age"].astype(float)
            sex = archive["sex"].astype(int)
            y = archive["y"].astype(int)
        rows = json.loads(meta_path.read_text())
        by_id = {row["gsm"]: row for row in rows}
        if len(by_id) != len(rows) or set(by_id) != set(sample_ids):
            raise ValueError(f"metadata/sample mismatch for {accession}")
        ordered = [by_id[sample] for sample in sample_ids]
        expected_age = np.asarray([float(row["age"]) for row in ordered])
        expected_sex = np.asarray([
            {1: 1, 0: 0, "M": 1, "Male": 1, "F": 0, "Female": 0}[row["sex"]]
            for row in ordered
        ])
        expected_y = np.asarray([int(row["case"]) for row in ordered])
        if not (np.array_equal(age, expected_age)
                and np.array_equal(sex, expected_sex)
                and np.array_equal(y, expected_y)):
            raise ValueError(f"metadata/array disagreement for {accession}")
        donor = np.asarray([
            accession + ":" + str(row["subject_id"] or row["gsm"])
            for row in ordered
        ])
        parts.append(SplitData(age, sex, y, np.repeat(accession, len(y)), donor,
                               sample_ids))

    data = SplitData(*(np.concatenate([getattr(part, field) for part in parts])
                       for field in ("age", "sex", "y", "source", "donor", "sample")))
    return data, hashes


def match_ids(data, ids):
    """Return deterministic greedy 1:1 pairs and unmatched sample IDs."""
    ids = set(map(str, ids))
    ix = np.flatnonzero(np.isin(data.sample, list(ids)))
    edges = []
    for source in sorted(set(data.source[ix])):
        source_ix = ix[data.source[ix] == source]
        for sex in sorted(set(data.sex[source_ix])):
            cell = source_ix[data.sex[source_ix] == sex]
            cases = [int(i) for i in cell[data.y[cell] == 1]]
            controls = [int(i) for i in cell[data.y[cell] == 0]]
            for case in cases:
                for control in controls:
                    gap = abs(float(data.age[case]) - float(data.age[control]))
                    if gap <= CALIPER_YEARS:
                        edges.append((gap, str(data.donor[case]),
                                      str(data.donor[control]),
                                      str(data.sample[case]), str(data.sample[control]),
                                      case, control))
    edges.sort(key=lambda edge: edge[:5])
    used_donors, matched_ix, pair_gaps, pairs = set(), set(), [], []
    for gap, case_donor, ctrl_donor, case_id, ctrl_id, case_ix, ctrl_ix in edges:
        if case_donor in used_donors or ctrl_donor in used_donors:
            continue
        used_donors.update((case_donor, ctrl_donor))
        matched_ix.update((case_ix, ctrl_ix))
        pair_gaps.append(gap)
        pairs.append({"case_sample": case_id, "control_sample": ctrl_id,
                      "case_donor": case_donor, "control_donor": ctrl_donor,
                      "absolute_age_gap_years": gap})
    selected = set(map(str, data.sample[ix]))
    retained = {str(data.sample[i]) for i in matched_ix}
    unmatched = sorted(selected - retained)
    n_cases = sum(int(data.y[i]) == 1 for i in matched_ix)
    n_controls = sum(int(data.y[i]) == 0 for i in matched_ix)
    gaps = np.asarray(pair_gaps, dtype=float)
    return {
        "input_n": int(len(ix)), "matched_pairs": len(pairs),
        "matched_n": int(len(matched_ix)), "matched_cases": n_cases,
        "matched_controls": n_controls, "unmatched_n": len(unmatched),
        "unmatched_sample_ids": unmatched, "age_gap_years": {
            "mean": float(gaps.mean()) if len(gaps) else None,
            "median": float(np.median(gaps)) if len(gaps) else None,
            "max": float(gaps.max()) if len(gaps) else None,
        },
        "minimum_class_n": MIN_CLASS_N,
        "passes_minimum_class_n": n_cases >= MIN_CLASS_N and n_controls >= MIN_CLASS_N,
        "pairs": pairs,
    }


def id_intersection(ids, allowed):
    allowed = set(allowed)
    return [sample for sample in ids if sample in allowed]


def main():
    data, hashes = metadata_only()
    frozen = geo.split_manifest(data, hashes)
    rows = []
    for entry in frozen["splits"]:
        scope = entry["scope"]
        outer_train = match_ids(data, entry["train"])
        train_s = {sample for pair in outer_train["pairs"]
                   for sample in (pair["case_sample"], pair["control_sample"])}
        row = {"scope": scope, "outer_train_before_matching": {
            "n": len(entry["train"]), "matched": outer_train},
            "outer_train_S_sample_n": len(train_s), "inner": []}
        if entry["test"]:
            outer_test = match_ids(data, entry["test"])
            row["outer_test"] = {"n": len(entry["test"]), "matched": outer_test}
            row["outer_test_evaluable"] = outer_test["passes_minimum_class_n"]
        else:
            row["outer_test"] = None
            row["outer_test_evaluable"] = None
        for number, inner in enumerate(entry["inner"]):
            inner_train_ids = id_intersection(inner["train"], train_s)
            inner_val_ids = id_intersection(inner["validation"], train_s)
            inner_train = match_ids(data, inner_train_ids)
            inner_val = match_ids(data, inner_val_ids)
            row["inner"].append({
                "fold": number,
                "train_restricted_to_S": {"n": len(inner_train_ids), "matched": inner_train},
                "validation_restricted_to_S": {"n": len(inner_val_ids), "matched": inner_val},
                "validation_evaluable": inner_val["passes_minimum_class_n"],
            })
        row["tuning_scope_evaluable"] = all(
            fold["validation_evaluable"] for fold in row["inner"])
        rows.append(row)

    manifest_bytes = json.dumps(frozen, sort_keys=True, separators=(",", ":")).encode()
    output = {
        "status": "metadata_count_audit_complete",
        "decision": "docs/GEO_AGE_SEX_MATCHING_DECISION_20261019_FA.md",
        "matching": {"caliper_years": CALIPER_YEARS, "minimum_class_n": MIN_CLASS_N,
                     "algorithm": "greedy_global_edges_no_replacement",
                     "tie_break": ["absolute_age_gap", "case_donor_id", "control_donor_id",
                                   "case_sample_id", "control_sample_id"],
                     "optimization_guarantee": "none; greedy is intentionally deterministic"},
        "frozen_split_manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "input_sha256": hashes,
        "PPMI_loaded": False, "models_fit": False, "scores_computed": False,
        "performance_metrics_computed": False,
        "scope_count": len(rows), "scopes": rows,
    }
    out_path = ROOT / "results" / "geo_age_sex_matching_audit_20261019.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    table_path = ROOT / "results" / "geo_age_sex_matching_audit_20261019.csv"
    with table_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            "scope", "partition", "n_before_matching", "n_after_matching",
            "matched_cases", "matched_controls", "matched_pairs", "unmatched_n",
            "age_gap_mean_years", "age_gap_median_years", "age_gap_max_years",
            "minimum_5_per_class_pass",
        ])
        writer.writeheader()
        for row in rows:
            train = row["outer_train_before_matching"]["matched"]
            writer.writerow({"scope": row["scope"], "partition": "outer_train",
                             "n_before_matching": train["input_n"],
                             "n_after_matching": train["matched_n"],
                             "matched_cases": train["matched_cases"],
                             "matched_controls": train["matched_controls"],
                             "matched_pairs": train["matched_pairs"],
                             "unmatched_n": train["unmatched_n"],
                             "age_gap_mean_years": train["age_gap_years"]["mean"],
                             "age_gap_median_years": train["age_gap_years"]["median"],
                             "age_gap_max_years": train["age_gap_years"]["max"],
                             "minimum_5_per_class_pass": "not_applicable"})
            if row["outer_test"] is not None:
                test = row["outer_test"]["matched"]
                writer.writerow({"scope": row["scope"], "partition": "outer_test",
                                 "n_before_matching": test["input_n"],
                                 "n_after_matching": test["matched_n"],
                                 "matched_cases": test["matched_cases"],
                                 "matched_controls": test["matched_controls"],
                                 "matched_pairs": test["matched_pairs"],
                                 "unmatched_n": test["unmatched_n"],
                                 "age_gap_mean_years": test["age_gap_years"]["mean"],
                                 "age_gap_median_years": test["age_gap_years"]["median"],
                                 "age_gap_max_years": test["age_gap_years"]["max"],
                                 "minimum_5_per_class_pass": test["passes_minimum_class_n"]})
            for inner in row["inner"]:
                for side, key in (("inner_train_restricted_to_S", "train_restricted_to_S"),
                                  ("inner_validation_restricted_to_S", "validation_restricted_to_S")):
                    matched = inner[key]["matched"]
                    writer.writerow({"scope": row["scope"],
                                     "partition": f"{side}_{inner['fold']}",
                                     "n_before_matching": matched["input_n"],
                                     "n_after_matching": matched["matched_n"],
                                     "matched_cases": matched["matched_cases"],
                                     "matched_controls": matched["matched_controls"],
                                     "matched_pairs": matched["matched_pairs"],
                                     "unmatched_n": matched["unmatched_n"],
                                     "age_gap_mean_years": matched["age_gap_years"]["mean"],
                                     "age_gap_median_years": matched["age_gap_years"]["median"],
                                     "age_gap_max_years": matched["age_gap_years"]["max"],
                                     "minimum_5_per_class_pass": (
                                         matched["passes_minimum_class_n"]
                                         if side.startswith("inner_validation")
                                         else "not_applicable")})
    print(json.dumps({
        "status": output["status"], "output": str(out_path), "table": str(table_path),
        "frozen_split_manifest_sha256": output["frozen_split_manifest_sha256"],
        "scope_count": len(rows), "PPMI_loaded": False, "models_fit": False,
        "scores_computed": False, "performance_metrics_computed": False,
        "scopes": [{
            "scope": row["scope"],
            "train_pairs": row["outer_train_before_matching"]["matched"]["matched_pairs"],
            "train_unmatched": row["outer_train_before_matching"]["matched"]["unmatched_n"],
            "test_pairs": row["outer_test"]["matched"]["matched_pairs"] if row["outer_test"] else None,
            "test_minimum_pass": row["outer_test_evaluable"],
            "inner_validation_pairs": [f["validation_restricted_to_S"]["matched"]["matched_pairs"]
                                       for f in row["inner"]],
            "inner_validation_minimum_pass": [f["validation_evaluable"] for f in row["inner"]],
            "tuning_scope_evaluable": row["tuning_scope_evaluable"],
        } for row in rows]
    }, indent=2))


if __name__ == "__main__":
    main()
