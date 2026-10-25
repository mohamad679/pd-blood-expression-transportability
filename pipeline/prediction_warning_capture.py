"""Attach Python warnings emitted by model.predict to the active fit context."""
from __future__ import annotations

import json
import warnings
from functools import wraps
from pathlib import Path


def instrument_predict(estimator, context, log_path):
    """Wrap one estimator's predict call without changing returned predictions.

    Calling this again updates the context on an already-instrumented model,
    which lets the same fitted model be tagged for its later held-out call.
    """
    capture = getattr(estimator, "_prediction_warning_capture", None)
    if capture is not None:
        capture["context"] = dict(context)
        capture["log_path"] = Path(log_path)
        return estimator

    original = estimator.predict
    capture = {"context": dict(context), "log_path": Path(log_path)}

    @wraps(original)
    def predict_with_warning_capture(*args, **kwargs):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = original(*args, **kwargs)
        for warning in caught:
            row = dict(capture["context"])
            row.update(warning_category=warning.category.__name__,
                       message=str(warning.message), filename=warning.filename,
                       lineno=warning.lineno)
            destination = capture["log_path"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            with destination.open("a", encoding="utf-8") as stream:
                stream.write(json.dumps(row, sort_keys=True) + "\n")
                stream.flush()
            # Preserve the existing visible stderr behavior while retaining a
            # structured record that can be tied to the fit/fold/candidate.
            warnings.showwarning(warning.message, warning.category,
                                 warning.filename, warning.lineno)
        return result

    estimator.predict = predict_with_warning_capture
    estimator._prediction_warning_capture = capture
    return estimator
