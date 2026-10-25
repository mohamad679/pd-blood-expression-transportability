"""Conservative numerical provenance candidates for phase-1 review.

PPMI trace rows remain under ~/PPMI_private; the public CSV covers GEO results
only. A Y means a numeric source candidate matches printed precision, not that
the statistic's scientific provenance has been independently established.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRIVATE = Path.home() / "PPMI_private"
MANUSCRIPT = PRIVATE / "manuscript_20261021_revision1/MANUSCRIPT_PRIVATE_REVISED.md"
NUMBER = re.compile(r"(?<![A-Za-z0-9])(?P<value>[+−-]?(?:\d+\.\d+|\d+))(?![A-Za-z0-9])")
GEO_SOURCES = [
    ROOT / "results/geo_development_20261018/aggregate_summary.json",
    ROOT / "results/geo_transfer_descriptive_bootstrap_20261020.json",
    ROOT / "results/geo_within_batch_auc_20261020/SUMMARY.json",
    ROOT / "results/geo_within_batch_auc_20261020/batch_results.csv",
    ROOT / "results/geo_batch_only_diagnostic_20261019/SUMMARY.json",
    ROOT / "results/geo_controls_20261018/aggregate_summary.json",
    ROOT / "results/geo_boundary_selection_inventory_20261020.json",
    ROOT / "results/geo_age_sex_matching_audit_20261019.json",
    ROOT / "results/geo_technical_metadata_audit_20261019/RESULTS.json",
    ROOT / "results/geo_combat_numeric_audit_20261020/SUMMARY.json",
    ROOT / "revision/outputs/within_batch_fold_summary.json",
    ROOT / "revision/outputs/PROBE_chip_proxy_fold_AUROC.csv",
    ROOT / "revision/outputs/final_model_specs.csv",
    ROOT / "revision/outputs/PROBE_chip_proxy_summary.json",
]
PPMI_SOURCES = [
    PRIVATE / "score_20261021/AGGREGATE_RESULTS.json",
    PRIVATE / "secondary_20261021_retry2/AGGREGATE_SECONDARY.json",
    PRIVATE / "revision_outputs/candidate_characteristics_P1.json",
    PRIVATE / "revision_outputs/prior_descriptive_ci_aggregate_20261021.json",
    ROOT / "private_data/ppmi_candidate_collection_timing_20261020.json",
    ROOT / "private_data/ppmi_ledd_at_collection_20261020.json",
]


def numeric_leaves(value, path=""):
    if isinstance(value, dict):
        for key, item in value.items():
            yield from numeric_leaves(item, f"{path}.{key}" if path else key)
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from numeric_leaves(item, f"{path}[{index}]")
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        yield path, float(value)


def source_index(paths):
    index = []
    for path in paths:
        if not path.is_file():
            continue
        if path.suffix == ".json":
            for key, value in numeric_leaves(json.loads(path.read_text())):
                index.append((str(path), key, value))
        elif path.suffix == ".csv":
            with path.open(newline="") as handle:
                for n, row in enumerate(csv.DictReader(handle), 2):
                    for key, raw in row.items():
                        try:
                            value = float(raw)
                        except (TypeError, ValueError):
                            continue
                        index.append((str(path), f"row{n}.{key}", value))
    return index


def trace(lines, index, eligible):
    rows = []
    for line_no, line in enumerate(lines, 1):
        if not eligible(line_no):
            continue
        # Remove citation numbers and DOI/URL material, which are not results.
        cleaned = re.sub(r"\[\d+(?:\s*,\s*\d+)*\]", "", line)
        cleaned = re.sub(r"https?://\S+", "", cleaned)
        for match in NUMBER.finditer(cleaned):
            token = match.group("value")
            try:
                target = float(token.replace("−", "-"))
            except ValueError:
                continue
            decimals = len(token.split(".", 1)[1]) if "." in token else 0
            tolerance = 0.5 * 10 ** (-decimals) + 1e-9 if decimals else 1e-9
            hits = [(path, key, value) for path, key, value in index if abs(value - target) <= tolerance]
            # Numeric matches only nominate a source; manual statistic-level check follows.
            best = hits[0] if hits else None
            rows.append({"number": token, "location": f"line {line_no}",
                         "source_file": best[0] if best else "",
                         "source_field": best[1] if best else "",
                         "source_value": best[2] if best else "",
                         "match": "Y" if best else "N",
                         "review_note": "precision-level candidate; manual verification required" if best else
                         "no numeric source candidate in indexed aggregate outputs"})
    return rows


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("number", "location", "source_file", "source_field",
                                                   "source_value", "match", "review_note"))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--private-output", type=Path, required=True)
    args = parser.parse_args()
    out = args.private_output.resolve()
    if not out.is_relative_to(PRIVATE.resolve()):
        raise ValueError("Private number trace must remain under ~/PPMI_private")
    lines = MANUSCRIPT.read_text().splitlines()
    # Mixed GEO/PPMI prose and captions remain private as whole lines.
    mixed_private_lines = {87, 134}
    geo = trace(lines, source_index(GEO_SOURCES),
                lambda n: (83 <= n <= 92 or 123 <= n <= 149) and n not in mixed_private_lines)
    ppmi = trace(lines, source_index(PPMI_SOURCES),
                 lambda n: 93 <= n <= 102 or 150 <= n <= 160 or n in mixed_private_lines)
    # Full private census covers abstract, Methods, Results, tables, captions,
    # Discussion and declarations; references are citation metadata, not results.
    full = trace(lines, source_index(GEO_SOURCES + PPMI_SOURCES), lambda n: 9 <= n < 189)
    write(ROOT / "revision/outputs/number_trace_P1.csv", geo)
    write(out, ppmi)
    full_path = out.parent / "number_trace_P1_FULL_PRIVATE.csv"
    write(full_path, full)
    out.chmod(0o600)
    full_path.chmod(0o600)
    print(f"GEO trace: {len(geo)} numbers, {sum(r['match']=='N' for r in geo)} unresolved; "
          f"PPMI trace: {len(ppmi)} numbers, {sum(r['match']=='N' for r in ppmi)} unresolved (private); "
          f"full manuscript: {len(full)} numbers, {sum(r['match']=='N' for r in full)} unresolved (private)")


if __name__ == "__main__":
    main()
