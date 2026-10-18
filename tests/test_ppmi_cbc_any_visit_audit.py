import csv
import tempfile
import unittest
from pathlib import Path

from pipeline.audit_ppmi_cbc_any_visit import audit_table, read_baseline_collection


class AnyVisitCoverageTests(unittest.TestCase):
    def test_collection_date_does_not_require_bldcoll(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "collection.csv"
            path.write_text("PATNO,EVENT_ID,BLDCLLDT\n1,BL,01/2020\n", encoding="utf-8")
            dates, codes = read_baseline_collection(path, {("1", "BL")})
        self.assertEqual(dates, {("1", "BL"): (2020, 1)})
        self.assertFalse(codes)

    def test_screening_same_month_counts_without_results_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "metadata_only.csv"
            with path.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["PATNO", "EVENT_ID", "LCOLLDT"])
                writer.writerows([
                    [1, "SC", "01/2020"], [1, "SC", "01/2020"],
                    [2, "SC", "12/2019"], [2, "BL", "01/2020"],
                    [3, "SC", "01/2020"], [99, "SC", "01/2020"],
                ])
            result = audit_table(path, {"1": (2020, 1), "2": (2020, 1), "3": None})
        self.assertEqual(result["eligible_participants_with_any_row_any_visit"], 3)
        self.assertEqual(result["eligible_participants_with_any_row_same_BL_collection_month"], 2)
        self.assertEqual(result["same_month_SC_priority_participants_without_BL"], 1)
        self.assertEqual(result["same_month_BL_priority_participants"], 1)
        self.assertEqual(result["by_EVENT_ID"]["SC"]["CBC_table_rows"], 4)
        self.assertEqual(result["by_EVENT_ID"]["SC"]["same_BL_collection_month_unique_participants"], 1)


if __name__ == "__main__":
    unittest.main()
