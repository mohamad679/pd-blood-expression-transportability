#!/usr/bin/env python3
"""Aggregate-only metadata audit of a private pre-lock PPMI candidate bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import zipfile
from collections import Counter, defaultdict
from pathlib import Path


def clean(value: object) -> str:
    return " ".join(str(value or "").strip().lower().replace("_", " ").split())


def normalized_patno(value: object) -> str:
    text = str(value or "").strip()
    return text[3:] if text.upper().startswith("PP-") else text


def audit(bundle: Path) -> dict:
    if not bundle.is_file():
        raise ValueError("Private bundle was not found")
    with zipfile.ZipFile(bundle) as archive:
        if archive.testzip() is not None:
            raise ValueError("Private bundle ZIP integrity failed")
        expected = {"MANIFEST.json", "rank_features.f32", "candidate_metadata.jsonl",
                    "common_genes_gencode_v29.tsv"}
        if set(archive.namelist()) != expected:
            raise ValueError("Private bundle has unexpected or missing files")
        manifest = json.loads(archive.read("MANIFEST.json"))
        if manifest.get("scores_computed") is not False or manifest.get("models_fitted") is not False:
            raise ValueError("Bundle manifest does not assert score-free preparation")
        if manifest.get("candidate_count") != 472 or manifest.get("gene_count") != 16364:
            raise ValueError("Bundle manifest dimensions differ from the fixed audit")
        if manifest.get("feature_shape") != [472, 16364] or manifest.get("feature_dtype") != "little-endian float32":
            raise ValueError("Bundle feature schema differs from the fixed audit")
        for name, expected_hash in manifest["files_sha256"].items():
            if hashlib.sha256(archive.read(name)).hexdigest() != expected_hash:
                raise ValueError("Internal file hash mismatch")
        if archive.getinfo("rank_features.f32").file_size != 472 * 16364 * 4:
            raise ValueError("Feature byte size does not match manifest shape")
        lines = archive.read("candidate_metadata.jsonl").decode("utf-8").splitlines()
        if len(lines) != 472:
            raise ValueError("Candidate metadata row count differs from feature rows")
        rows = [json.loads(line) for line in lines]
    participants: set[str] = set()
    samples: set[str] = set()
    release_rows: set[int] = set()
    groups = Counter()
    plate_groups: dict[str, Counter[str]] = defaultdict(Counter)
    portal_plate_compare = Counter()
    for row in rows:
        release = row["release_metadata"]
        portal = row["portal_metadata"]
        participant = normalized_patno(row["PATNO"])
        if not participant or participant in participants:
            raise ValueError("Candidate participant identifiers are missing or repeated")
        participants.add(participant)
        sample = str(release.get("Sample") or "").strip()
        if not sample or sample in samples:
            raise ValueError("Candidate release sample identifiers are missing or repeated")
        samples.add(sample)
        release_i = row["release_row_index"]
        if not isinstance(release_i, int) or release_i in release_rows:
            raise ValueError("Candidate release row indexes are invalid or repeated")
        release_rows.add(release_i)
        if any(normalized_patno(source.get("PATNO")) != participant for source in (release, portal)):
            raise ValueError("Participant identifier disagrees across linked metadata")
        if row["EVENT_ID"] != "BL" or clean(release.get("CLINICAL_EVENT")) != "bl" or clean(portal.get("Clinical Event")) != "bl":
            raise ValueError("Candidate is not baseline in all metadata")
        if clean(release.get("QCflagIR3")) != "pass" or clean(portal.get("QC")) != "highqc" or clean(portal.get("Long RNA-seq")) != "yes":
            raise ValueError("Candidate does not meet release and portal Long RNA QC rules")
        age = float(row["age_at_BL_visit"])
        if not math.isfinite(age) or not 0 < age < 120:
            raise ValueError("Candidate baseline age is missing or invalid")
        sex = clean(row["sex"])
        if sex not in {"male", "female"} or sex != clean(release.get("GENDER")) or sex != clean(portal.get("Sex")):
            raise ValueError("Candidate recorded sex disagrees across metadata")
        label = row["historical_label"]
        if label == "PD":
            wanted = ("idiopathic pd", "pd", "case", "pd")
        elif label == "Control":
            wanted = ("healthy control", "control", "control", "control")
        else:
            raise ValueError("Candidate historical label is unrecognized")
        actual = (clean(portal.get("Disease Status")), clean(portal.get("Diagnosis")),
                  clean(portal.get("Case Control")), clean(release.get("DIAGNOSIS")))
        if actual != wanted:
            raise ValueError("Candidate historical labels disagree across metadata")
        groups[label] += 1
        plate = str(release.get("Plate") or "").strip()
        if not plate or plate != row["plate"] or row.get("plate_source") != "metaDataIR3_original":
            raise ValueError("Release plate is absent or not the selected plate")
        plate_groups[plate][label] += 1
        portal_plate = str(portal.get("Plate") or "").strip()
        if portal_plate != row["portal_plate"]:
            raise ValueError("Portal plate was not retained verbatim")
        portal_plate_compare["missing" if not portal_plate else
                             "same_literal_value" if portal_plate == plate else
                             "different_literal_value"] += 1
    if groups != {"PD": 311, "Control": 161}:
        raise ValueError("Candidate group counts differ from prior aggregate audit")
    return {
        "status": "core_candidate_metadata_checks_passed_final_population_pending",
        "scope": "Pre-second-lock metadata audit only; no PPMI scores or model fitting",
        "candidate_count": len(rows), "group_counts": dict(groups),
        "unique_participants": len(participants), "unique_release_samples": len(samples),
        "unique_release_rows": len(release_rows),
        "plate_source": "metaDataIR3_original Plate",
        "plate_count": len(plate_groups),
        "plates_with_both_groups": sum(x["PD"] > 0 and x["Control"] > 0 for x in plate_groups.values()),
        "plate_by_group_counts": [
            {"plate": plate, "PD": counts["PD"], "Control": counts["Control"],
             "total": counts["PD"] + counts["Control"]}
            for plate, counts in sorted(plate_groups.items())
        ],
        "portal_plate_literal_comparison": dict(portal_plate_compare),
        "bundle_sha256": hashlib.sha256(bundle.read_bytes()).hexdigest(),
        "internal_file_hashes_match": True,
        "core_checks": "baseline, Long RNA, release and portal QC, age, sex, historical label agreement, exact one-to-one metadata rows",
        "remaining_before_final_population": [
            "Pre-QC exclusion counts by group require full release and portal metadata, not only selected candidates.",
            "Medication at sampling and blood-draw timing require independent clinical/collection metadata.",
            "Computational one-to-one sample linkage does not independently establish biological donor identity.",
            "The second public Zenodo lock has not been published; no PPMI score is authorized yet."
        ],
        "privacy": "Aggregate-only report without participant or sample identifiers or individual features",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = audit(args.bundle)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise ValueError("Audit output already exists; do not overwrite an earlier attempt")
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
