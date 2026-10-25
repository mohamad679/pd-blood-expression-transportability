"""Tiny synthetic check that predict warnings retain context and outputs."""
import json
import tempfile
import unittest
import warnings
from pathlib import Path

import numpy as np

from pipeline.prediction_warning_capture import instrument_predict


class WarningEstimator:
    def predict(self, values):
        warnings.warn("synthetic clipping warning", RuntimeWarning)
        return np.asarray(values) * 2


class PredictionWarningCaptureTest(unittest.TestCase):
    def test_warning_is_logged_with_fit_context_and_prediction_is_unchanged(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "warnings.jsonl"
            context = dict(scope="GEO/outer/2", stage="inner", fold=3,
                           model="M5", **{"lambda": 0.01, "rho": 1.0})
            model = instrument_predict(WarningEstimator(), context, path)
            with warnings.catch_warnings(record=True) as visible:
                warnings.simplefilter("always")
                output = model.predict([1, 2])
            self.assertTrue(np.array_equal(output, [2, 4]))
            self.assertEqual(len(visible), 1)
            row = json.loads(path.read_text().strip())
            for key, value in context.items():
                self.assertEqual(row[key], value)
            self.assertEqual(row["warning_category"], "RuntimeWarning")
            self.assertEqual(row["message"], "synthetic clipping warning")

            updated = dict(context, stage="outer_holdout_prediction", fold=None)
            instrument_predict(model, updated, path)
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                model.predict([3])
            rows = [json.loads(line) for line in path.read_text().splitlines()]
            self.assertEqual(rows[-1]["stage"], "outer_holdout_prediction")
            self.assertIsNone(rows[-1]["fold"])


if __name__ == "__main__":
    unittest.main()
