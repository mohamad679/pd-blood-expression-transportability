"""End-to-end synthetic check of the score-free IR3 cache export."""

import argparse
import csv
import hashlib
import io
import json
import struct
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
sys.path.insert(0, str(ROOT / "pipeline" / "windows_ppmi"))
import prepare_ppmi_remaining_inputs as one  # noqa: E402
import prepare_ppmi_windows_bundle as original  # noqa: E402


def csv_bytes(rows):
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


class OneShotExportTest(unittest.TestCase):
    def test_complete_cache_to_aggregate_and_two_private_bundles(self):
        release = []
        portal = []
        for i, (pid, group, event, sex) in enumerate([
            ("1", "PD", "BL", "Male"),
            ("2", "Control", "BL", "Female"),
            ("1", "PD", "V02", "Male"),
            ("2", "Control", "V02", "Female"),
        ], 1):
            release.append({
                "PATNO": pid, "Sample": f"S{i}", "HudAlphaID": f"S{i}",
                "CLINICAL_EVENT": event, "DIAGNOSIS": group, "QCflagIR3": "Pass",
                "GENDER": sex, "Plate": "1",
            })
            if event == "BL":
                portal.append({
                    "PATNO": pid, "HudAlphaSampleName": f"S{i}", "Clinical Event": "BL",
                    "Disease Status": "Idiopathic PD" if group == "PD" else "Healthy Control",
                    "Diagnosis": group, "Case Control": "Case" if group == "PD" else "Control",
                    "Long RNA-seq": "Yes", "QC": "HighQC", "Sex": sex,
                })
        ages = [{"PATNO": pid, "EVENT_ID": "BL", "AGE_AT_VISIT": "60"} for pid in ("1", "2")]
        statuses = [{"PATNO": pid} for pid in ("1", "2")]
        index = [{"cache_entry": f"counts/{i:06d}.f64", "count_header": f"S{i+1}",
                  "member_name": f"S{i+1}"} for i in range(4)]
        reference = ("gene_id\tsymbol\tchromosome\tbiotype\texon_union_length_bp\n"
                     "G1\tA\t1\tx\t1\nG2\tB\t1\tx\t1\n"
                     "G3\tC\t1\tx\t1\nG4\tD\t1\tx\t1\n").encode()
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            cache = base / "IR3_PRIVATE_CACHE.zip"
            with zipfile.ZipFile(cache, "w") as archive:
                archive.writestr(one.RELEASE_MEMBER, csv_bytes(release))
                archive.writestr("inputs/" + one.PORTAL, csv_bytes(portal))
                archive.writestr("inputs/" + one.AGE, csv_bytes(ages))
                archive.writestr("inputs/" + one.STATUS, csv_bytes(statuses))
                archive.writestr(one.INDEX_MEMBER, json.dumps(index))
                archive.writestr("inputs/" + one.GENES, reference)
                for i in range(4):
                    archive.writestr(f"counts/{i:06d}.f64", struct.pack("<4d", 0, 1, 2, 3))
            output = base / "out"
            with mock.patch.object(one, "EXPECTED_RELEASE_ROWS", 4), \
                    mock.patch.object(one, "EXPECTED_CANDIDATES", 2), \
                    mock.patch.object(one, "EXPECTED_GROUPS", {"PD": 1, "Control": 1}), \
                    mock.patch.object(original, "EXPECTED_COUNT_FILES", 4), \
                    mock.patch.object(original, "GENE_COUNT", 4), \
                    mock.patch.object(original, "REFERENCE_SHA256", hashlib.sha256(reference).hexdigest()):
                one.run(argparse.Namespace(cache=cache, capsule=None, candidate_bundle=None,
                                           output_dir=output))
                one.run(argparse.Namespace(cache=None, capsule=output / one.CAPSULE_NAME,
                                           candidate_bundle=None, output_dir=base / "mac_verification"))
            report = json.loads((output / one.REPORT_NAME).read_text())
            mac_report = json.loads((base / "mac_verification" / one.REPORT_NAME).read_text())
            self.assertEqual(mac_report["source_kind"], "metadata_capsule_verification")
            self.assertEqual(mac_report["preqc"], report["preqc"])
            self.assertEqual(report["selection"]["candidate_historical_groups"],
                             {"PD": 1, "Control": 1})
            self.assertEqual(report["preqc"]["groups"]["PD"]["stages"]["strict_linked_pre_qc_rows"], 1)
            self.assertEqual(report["all_visit_rank_preparation"]["sample_rows"], 4)
            with zipfile.ZipFile(output / one.RANK_BUNDLE_NAME) as archive:
                self.assertEqual(len(archive.read("rank_features.f32")), 4 * 4 * 4)
                self.assertFalse(json.loads(archive.read("MANIFEST.json"))["scores_computed"])
            with zipfile.ZipFile(output / one.CAPSULE_NAME) as archive:
                self.assertEqual(set(one.SOURCE_MEMBERS) | {"MANIFEST.json"}, set(archive.namelist()))


if __name__ == "__main__":
    unittest.main()
