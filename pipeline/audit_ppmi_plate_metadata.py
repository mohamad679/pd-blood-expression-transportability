#!/usr/bin/env python3
"""Aggregate-only PPMI plate crosswalk and descriptive in-sample AUROC."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import zipfile

from audit_ppmi_candidate_bundle import audit as audit_candidate_bundle


EXPECTED_BUNDLE_SHA256 = "1b191f042885eedf2eea00c69ce79a269cc647ef4d5148b37b678b97c7d0333c"


def canonical_plate(value: object, *, portal: bool) -> str:
    pattern = r"p([1-9][0-9]*)" if portal else r"([1-9][0-9]*)"
    match = re.fullmatch(pattern, str(value).strip(), flags=re.IGNORECASE)
    if match is None:
        raise ValueError("A plate value does not match the prespecified release/portal syntax")
    return str(int(match.group(1)))


def calculate(bundle: Path) -> dict:
    core = audit_candidate_bundle(bundle)
    if core["bundle_sha256"] != EXPECTED_BUNDLE_SHA256:
        raise ValueError("Private bundle differs from the prespecified candidate set")
    with zipfile.ZipFile(bundle) as archive:
        rows = [json.loads(line) for line in archive.read("candidate_metadata.jsonl").decode("utf-8").splitlines()]

    counts: dict[str, Counter[str]] = defaultdict(Counter)
    raw_pairs: Counter[tuple[str, str]] = Counter()
    release_to_portal: dict[str, set[str]] = defaultdict(set)
    portal_to_release: dict[str, set[str]] = defaultdict(set)
    labels_and_plates: list[tuple[str, str]] = []
    for row in rows:
        release = str(row["plate"]).strip()
        portal = str(row["portal_plate"]).strip()
        if canonical_plate(release, portal=False) != canonical_plate(portal, portal=True):
            raise ValueError("Release and portal plates disagree after prespecified normalization")
        release_to_portal[release].add(portal)
        portal_to_release[portal].add(release)
        raw_pairs[(release, portal)] += 1
        label = row["historical_label"]
        counts[release][label] += 1
        labels_and_plates.append((label, release))
    if any(len(values) != 1 for values in release_to_portal.values()) or any(
        len(values) != 1 for values in portal_to_release.values()
    ):
        raise ValueError("Observed release-to-portal plate mapping is not one-to-one")
    if len(counts) != core["plate_count"] or sum(sum(value.values()) for value in counts.values()) != core["candidate_count"]:
        raise ValueError("Plate counts do not reconcile with the candidate audit")
    if counts.get("201") != {"PD": 1}:
        raise ValueError("The recorded one-person plate 201 has changed")

    prevalence = {plate: Fraction(value["PD"], sum(value.values())) for plate, value in counts.items()}
    pd_scores = [prevalence[plate] for label, plate in labels_and_plates if label == "PD"]
    control_scores = [prevalence[plate] for label, plate in labels_and_plates if label == "Control"]
    pair_count = len(pd_scores) * len(control_scores)
    twice_favorable = sum(2 if pd > control else 1 if pd == control else 0
                          for pd in pd_scores for control in control_scores)
    if (len(pd_scores), len(control_scores)) != (core["group_counts"]["PD"], core["group_counts"]["Control"]):
        raise ValueError("Class counts do not reconcile with the candidate audit")

    return {
        "status": "metadata_only_plate_diagnostic_complete_provisional_population",
        "candidate_count": len(rows),
        "group_counts": core["group_counts"],
        "bundle_sha256": core["bundle_sha256"],
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "plate_reference": "metaDataIR3_original Plate",
        "plate_count": len(counts),
        "portal_plate_mapping": "strip whitespace; portal pN and release N map to positive integer N",
        "rowwise_normalized_plate_agreement_count": len(rows),
        "one_to_one_mapping_both_directions": True,
        "plate_crosswalk": [
            {"release_plate": release, "portal_plate": portal, "candidate_count": n}
            for (release, portal), n in sorted(raw_pairs.items(), key=lambda item: int(item[0][0]))
        ],
        "plate_201": {"PD": counts["201"]["PD"], "Control": counts["201"]["Control"],
                      "retained": True, "purpose": "unknown"},
        "plate_only_in_sample_auroc": twice_favorable / (2 * pair_count),
        "pd_control_pair_count": pair_count,
        "twice_favorable_pair_count": twice_favorable,
        "interpretation": "Descriptive in-sample metadata-only association; not cross-validated or confirmatory",
        "ppmi_expression_or_model_scores_computed": False,
        "final_population_status": "pending",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Audit output already exists; do not overwrite an earlier attempt")
    result = calculate(args.bundle)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
