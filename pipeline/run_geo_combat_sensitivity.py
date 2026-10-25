"""Run the preregistered partial-ComBat sensitivity on frozen GEO splits.

This exploratory runner evaluates three registered arms only:
inductive without disease, inductive preserving disease, and transductive
all-expression without disease. PPMI is never loaded.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
import time
import argparse
from contextlib import redirect_stdout
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from neuroCombat import neuroCombat, neuroCombatFromTraining

import geo_pipeline as gp


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "private_data/geo_combat_sensitivity_20261020"
EXPECTED_SPLIT_SHA256 = "2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2"
PINNED_NEUROCOMBAT = "0.2.12"
MIN_TRAIN_BATCH_N = 3
ARMS = (
    "inductive_no_disease",
    "inductive_preserve_disease",
    "transductive_no_disease",
)
# Upper bound for the frozen train_suite grid: the initial 5 lambda values,
# plus at most one new boundary lambda in each of two expansion rounds.
MAX_CHECKED_FITS_PER_SCOPE = 330


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    temp.replace(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sample_hash(samples: np.ndarray) -> str:
    return gp.array_hash(np.asarray(samples, dtype=str))


def with_expression(data: gp.Data, expression: np.ndarray) -> gp.Data:
    return gp.Data(expression, data.age, data.sex, data.y, data.source,
                   data.donor, data.sample, data.genes, data.marker_indices)


def load_batch_lookup() -> tuple[dict[str, str], dict[str, str]]:
    lookup: dict[str, str] = {}
    hashes: dict[str, str] = {}
    for cohort in ("GSE57475", "GSE99039"):
        path = ROOT / "prepared" / f"{cohort}_sample_metadata.json"
        rows = json.loads(path.read_text())
        hashes[path.name] = sha256(path)
        for row in rows:
            sample = str(row["gsm"])
            batch = row.get("batch")
            if not batch or sample in lookup:
                raise ValueError(f"Missing/duplicate batch metadata for {sample}")
            # Prefix cohort because PROBE chip labels and GENEPARK batches are
            # cohort-specific technical factors, even if their text matches.
            lookup[sample] = f"{cohort}:{batch}"
    return lookup, hashes


def get_batches(data: gp.Data, batch_lookup: dict[str, str]) -> np.ndarray:
    try:
        return np.asarray([batch_lookup[str(sample)] for sample in data.sample], dtype=str)
    except KeyError as exc:
        raise ValueError(f"No technical batch for sample {exc.args[0]}") from exc


def eligible_mask(batch: np.ndarray) -> tuple[np.ndarray, dict[str, int]]:
    counts = {str(name): int(n) for name, n in pd.Series(batch).value_counts().items()}
    eligible = {name for name, count in counts.items() if count >= MIN_TRAIN_BATCH_N}
    mask = np.asarray([str(name) in eligible for name in batch], dtype=bool)
    return mask, counts


def disease_design_is_full_rank(batch: np.ndarray, disease: np.ndarray) -> bool:
    mask, _ = eligible_mask(batch)
    batches, labels = batch[mask].astype(str), disease[mask].astype(int)
    levels = sorted(set(batches))
    if len(levels) < 2 or len(set(labels)) < 2:
        return False
    batch_columns = [(batches == level).astype(float) for level in levels[1:]]
    design = np.column_stack([np.ones(len(labels)), *batch_columns, labels.astype(float)])
    return int(np.linalg.matrix_rank(design)) == design.shape[1]


def precheck_train_only(train: gp.Data, inner: list, batch_lookup: dict[str, str],
                        preserve_disease: bool) -> list[dict]:
    partitions = [(f"inner/{i}", train.take(ti)) for i, (ti, _) in enumerate(inner)]
    partitions.append(("outer_refit", train))
    records = []
    for name, part in partitions:
        batch = get_batches(part, batch_lookup)
        keep, counts = eligible_mask(batch)
        eligible_batches = len(set(batch[keep]))
        issue = None
        if eligible_batches < 2:
            issue = "fewer_than_two_eligible_training_batches"
        elif preserve_disease and not disease_design_is_full_rank(batch, part.y):
            issue = "batch_disease_design_not_full_rank"
        records.append({
            "partition": name,
            "training_sample_count": int(len(part.y)),
            "eligible_batch_count": int(eligible_batches),
            "corrected_training_sample_count": int(keep.sum()),
            "uncorrected_training_sample_count": int((~keep).sum()),
            "batch_counts": counts,
            "disease_design_full_rank": (
                disease_design_is_full_rank(batch, part.y) if preserve_disease else None
            ),
            "not_evaluable_reason": issue,
        })
    return records


def fit_partial_combat(x: np.ndarray, batch: np.ndarray,
                       disease: np.ndarray | None, log_path: Path) -> dict:
    keep, counts = eligible_mask(batch)
    eligible = set(batch[keep].astype(str))
    if len(eligible) < 2:
        raise ValueError("At least two eligible batches are required for ComBat")
    covars = pd.DataFrame({"batch": batch[keep].astype(str)})
    categorical = []
    if disease is not None:
        covars["disease"] = disease[keep].astype(int)
        categorical = ["disease"]
    # neuroCombat's informational stdout is voluminous. Save it per arm while
    # leaving numerical RuntimeWarnings on stderr so they remain visible.
    with log_path.open("a", encoding="utf-8") as log, redirect_stdout(log):
        result = neuroCombat(
            dat=x[keep].T,
            covars=covars,
            batch_col="batch",
            categorical_cols=categorical,
            eb=True,
            parametric=True,
            mean_only=False,
            ref_batch=None,
        )
    return {"keep": keep, "counts": counts, "estimates": result["estimates"]}


def apply_partial_combat(x: np.ndarray, batch: np.ndarray, state: dict,
                         log_path: Path) -> tuple[np.ndarray, np.ndarray]:
    known = set(np.asarray(state["estimates"]["batches"], dtype=str))
    corrected = np.asarray([str(name) in known for name in batch], dtype=bool)
    result = x.copy()
    if corrected.any():
        with log_path.open("a", encoding="utf-8") as log, redirect_stdout(log):
            transformed = neuroCombatFromTraining(
                dat=x[corrected].T,
                batch=batch[corrected].astype(str),
                estimates=state["estimates"],
            )["data"]
        result[corrected] = transformed.T
    return result, corrected


class NotEvaluable(RuntimeError):
    pass


class Progress:
    def __init__(self, total_units: int, completed_units: int = 0):
        self.total_units = total_units
        self.completed_units = completed_units
        self.started = time.monotonic()
        self.last_printed_at = 0.0
        self.is_tty = sys.stdout.isatty()

    def update(self, label: str, local_fits: int = 0, force: bool = False) -> None:
        local_fraction = (min(local_fits, MAX_CHECKED_FITS_PER_SCOPE) /
                          MAX_CHECKED_FITS_PER_SCOPE) * 0.95
        fraction = min(1.0, (self.completed_units + local_fraction) / self.total_units)
        now = time.monotonic()
        # A time-based refresh keeps the bar alive even when one fit takes a
        # long time, without writing a terminal line for every fit.
        if not force and now - self.last_printed_at < 5.0:
            return
        elapsed = now - self.started
        eta_text = (f"~{elapsed / fraction * (1 - fraction) / 60:.1f}m"
                    if fraction > 0 else "calculating")
        width = 24
        filled = round(width * fraction)
        bar = "#" * filled + "-" * (width - filled)
        fits_text = f" | fits {local_fits}/{MAX_CHECKED_FITS_PER_SCOPE}" if local_fits else ""
        message = (f"Progress [{bar}] {fraction * 100:5.1f}% | "
                   f"{self.completed_units}/{self.total_units} arm-scopes{fits_text} | "
                   f"elapsed {elapsed / 60:.1f}m | ETA {eta_text} | {label}")
        if self.is_tty and not force:
            print("\r" + message[:220], end="", flush=True)
        else:
            if self.is_tty:
                print("\r" + message[:220], flush=True)
            else:
                print(message, flush=True)
        self.last_printed_at = now

    def complete_unit(self, label: str) -> None:
        self.completed_units += 1
        if self.is_tty:
            print("\r" + " " * 220 + "\r", end="", flush=True)
        self.update(label, force=True)


def scope_data(split: dict, all_data: gp.Data, lookup: dict[str, int]) -> tuple[gp.Data, gp.Data, list]:
    train = all_data.take(np.asarray([lookup[s] for s in split["train"]], dtype=int))
    test = all_data.take(np.asarray([lookup[s] for s in split["test"]], dtype=int))
    local = {str(sample): i for i, sample in enumerate(train.sample)}
    inner = [
        (np.asarray([local[s] for s in fold["train"]], dtype=int),
         np.asarray([local[s] for s in fold["validation"]], dtype=int))
        for fold in split["inner"]
    ]
    if len(inner) != 5:
        raise ValueError(f"Expected five frozen inner folds for {split['scope']}")
    return train, test, inner


def run_inductive(scope: str, arm: str, train: gp.Data, test: gp.Data,
                  inner: list, batch_lookup: dict[str, str], output: Path,
                  progress: Progress) -> dict:
    preserve = arm == "inductive_preserve_disease"
    preflight = precheck_train_only(train, inner, batch_lookup, preserve)
    blocked = [r for r in preflight if r["not_evaluable_reason"]]
    if blocked:
        raise NotEvaluable(json.dumps(blocked, ensure_ascii=False))

    original_prep = gp.Prep
    tracker: dict = {"fits": [], "applications": [], "fit_counter": 0}

    class CombatPrep:
        def __init__(self, fit_data, max_genes=1000):
            fit_number = tracker["fit_counter"]
            tracker["fit_counter"] += 1
            fit_label = (f"inner/{fit_number}" if fit_number < len(inner)
                         else "outer_refit")
            batch = get_batches(fit_data, batch_lookup)
            disease = fit_data.y if preserve else None
            log_path = output / "neurocombat.log"
            state = fit_partial_combat(fit_data.X, batch, disease, log_path)
            corrected_train, train_mask = apply_partial_combat(
                fit_data.X, batch, state, log_path
            )
            self.state = state
            self.fit_label = fit_label
            self.fit_hash = sample_hash(fit_data.sample)
            self.applied_hashes: set[str] = set()
            self.base = original_prep(with_expression(fit_data, corrected_train), max_genes)
            tracker["fits"].append({
                "partition": fit_label,
                "fit_sample_hash": self.fit_hash,
                "training_sample_count": int(len(fit_data.y)),
                "eligible_batch_count": int(len(state["estimates"]["batches"])),
                "corrected_training_sample_count": int(train_mask.sum()),
                "uncorrected_training_sample_count": int((~train_mask).sum()),
                "batch_counts": state["counts"],
                "selected_gene_count": int(len(self.base.selected)),
            })

        def __getattr__(self, name):
            return getattr(self.base, name)

        def transform(self, data, model):
            batch = get_batches(data, batch_lookup)
            corrected, mask = apply_partial_combat(
                data.X, batch, self.state, output / "neurocombat.log"
            )
            identity = sample_hash(data.sample)
            if identity not in self.applied_hashes:
                tracker["applications"].append({
                    "fit_partition": self.fit_label,
                    "fit_sample_hash": self.fit_hash,
                    "evaluation_sample_hash": identity,
                    "evaluation_sample_count": int(len(data.y)),
                    "corrected_evaluation_sample_count": int(mask.sum()),
                    "uncorrected_evaluation_sample_count": int((~mask).sum()),
                })
                self.applied_hashes.add(identity)
            return self.base.transform(with_expression(data, corrected), model)

        def save(self, path):
            self.base.save(path)

    fit_counter = {"n": 0}
    original_checked_fit = gp.checked_fit

    def monitored_checked_fit(*args, **kwargs):
        result = original_checked_fit(*args, **kwargs)
        fit_counter["n"] += 1
        progress.update(f"{arm} | {scope} | checked fits {fit_counter['n']}",
                        local_fits=fit_counter["n"])
        return result

    try:
        gp.Prep = CombatPrep
        gp.checked_fit = monitored_checked_fit
        prep, fits, decisions = gp.train_suite(train, inner, output, scope)
        predictions = gp.predict_suite(prep, fits, decisions, test) if len(test.y) else {}
    finally:
        gp.Prep = original_prep
        gp.checked_fit = original_checked_fit

    if predictions:
        np.savez_compressed(output / "heldout_predictions.npz",
                            samples=test.sample, y=test.y, **predictions)
    tracker["precheck"] = preflight
    tracker["fit_count"] = fit_counter["n"]
    tracker["expected_frozen_inner_fold_count"] = len(inner)
    tracker["PPMI_accessed"] = False
    tracker["heldout_diagnosis_used_for_correction"] = False
    tracker["heldout_expression_used_for_correction"] = False
    write_json(output / "combat_coverage.json", tracker)
    return {
        "status": "complete",
        "scope": scope,
        "arm": arm,
        "interpretation": "inductive train-only partial ComBat sensitivity",
        "seconds": None,
        "metrics": {name: gp.metrics(test.y, p) for name, p in predictions.items()},
        "coverage_file": "combat_coverage.json",
        "model_fit_count": fit_counter["n"],
        "PPMI_accessed": False,
    }


def run_transductive(scope: str, train: gp.Data, test: gp.Data,
                     split_data: gp.Data, batch_lookup: dict[str, str],
                     inner: list, output: Path, progress: Progress) -> dict:
    batch = get_batches(split_data, batch_lookup)
    eligible, counts = eligible_mask(batch)
    if len(set(batch[eligible])) < 2:
        raise NotEvaluable("fewer_than_two_eligible_batches_in_full_scope")
    # This arm intentionally uses all unlabeled expression in the scope,
    # including held-out samples; it is not an inductive estimate.
    log_path = output / "neurocombat.log"
    state = fit_partial_combat(split_data.X, batch, disease=None, log_path=log_path)
    corrected, correction_mask = apply_partial_combat(
        split_data.X, batch, state, log_path
    )
    progress.update(f"transductive ComBat correction complete | {scope}", force=True)
    index = {str(sample): i for i, sample in enumerate(split_data.sample)}
    corrected_train = with_expression(
        train, corrected[np.asarray([index[str(s)] for s in train.sample], dtype=int)]
    )
    corrected_test = with_expression(
        test, corrected[np.asarray([index[str(s)] for s in test.sample], dtype=int)]
    )
    original_checked_fit = gp.checked_fit
    fit_counter = {"n": 0}

    def monitored_checked_fit(*args, **kwargs):
        result = original_checked_fit(*args, **kwargs)
        fit_counter["n"] += 1
        progress.update(f"transductive_no_disease | {scope} | checked fits {fit_counter['n']}",
                        local_fits=fit_counter["n"])
        return result

    try:
        gp.checked_fit = monitored_checked_fit
        prep, fits, decisions = gp.train_suite(corrected_train, inner, output, scope)
        predictions = (gp.predict_suite(prep, fits, decisions, corrected_test)
                       if len(test.y) else {})
    finally:
        gp.checked_fit = original_checked_fit

    if predictions:
        np.savez_compressed(output / "heldout_predictions.npz",
                            samples=test.sample, y=test.y, **predictions)
    split_lookup = {str(s): i for i, s in enumerate(split_data.sample)}
    train_positions = np.asarray([split_lookup[str(s)] for s in train.sample], dtype=int)
    test_positions = np.asarray([split_lookup[str(s)] for s in test.sample], dtype=int)
    inner_coverage = []
    for fold_index, (fit_ix, validation_ix) in enumerate(inner):
        fit_ids = train.sample[fit_ix]
        validation_ids = train.sample[validation_ix]
        fit_positions = np.asarray([split_lookup[str(s)] for s in fit_ids], dtype=int)
        validation_positions = np.asarray([split_lookup[str(s)] for s in validation_ids], dtype=int)
        inner_coverage.append({
            "partition": f"inner/{fold_index}",
            "training_sample_hash": sample_hash(fit_ids),
            "validation_sample_hash": sample_hash(validation_ids),
            "training_sample_count": int(len(fit_positions)),
            "training_corrected_sample_count": int(correction_mask[fit_positions].sum()),
            "training_uncorrected_sample_count": int((~correction_mask[fit_positions]).sum()),
            "validation_sample_count": int(len(validation_positions)),
            "validation_corrected_sample_count": int(correction_mask[validation_positions].sum()),
            "validation_uncorrected_sample_count": int((~correction_mask[validation_positions]).sum()),
        })
    coverage = {
        "scope": scope,
        "arm": "transductive_no_disease",
        "full_scope_sample_hash": sample_hash(split_data.sample),
        "full_scope_sample_count": int(len(split_data.sample)),
        "full_scope_batch_counts": counts,
        "full_scope_corrected_sample_count": int(correction_mask.sum()),
        "full_scope_uncorrected_sample_count": int((~correction_mask).sum()),
        "training_corrected_sample_count": int(correction_mask[train_positions].sum()),
        "training_uncorrected_sample_count": int((~correction_mask[train_positions]).sum()),
        "evaluation_corrected_sample_count": int(correction_mask[test_positions].sum()),
        "evaluation_uncorrected_sample_count": int((~correction_mask[test_positions]).sum()),
        "inner_partition_coverage": inner_coverage,
        "heldout_diagnosis_used_for_correction": False,
        "heldout_expression_used_for_correction": bool(len(test.y)),
        "PPMI_accessed": False,
        "model_fit_count": fit_counter["n"],
    }
    write_json(output / "combat_coverage.json", coverage)
    return {
        "status": "complete",
        "scope": scope,
        "arm": "transductive_no_disease",
        "interpretation": "transductive all-expression partial ComBat; held-out expression informs correction",
        "metrics": {name: gp.metrics(test.y, p) for name, p in predictions.items()},
        "coverage_file": "combat_coverage.json",
        "model_fit_count": fit_counter["n"],
        "PPMI_accessed": False,
    }


def markdown_report(rows: list[dict], split_hash: str) -> str:
    lines = [
        "# GEO partial-ComBat sensitivity",
        "",
        "Exploratory GEO-only sensitivity using frozen folds. The three arms are",
        "inductive without disease, inductive preserving disease, and explicitly",
        "transductive all-expression without disease. Small batches remain uncorrected.",
        "ComBat input is the existing within-sample gene-rank matrix.",
        "",
        "| Scope | Arm | Status | M4 AUROC | M5 AUROC |",
        "|---|---|---|---:|---:|",
    ]
    for row in rows:
        metrics = row.get("metrics", {})
        m4 = metrics.get("M4", {}).get("AUROC")
        m5 = metrics.get("M5", {}).get("AUROC")
        m4_text = f"{m4:.3f}" if isinstance(m4, (int, float)) else "—"
        m5_text = f"{m5:.3f}" if isinstance(m5, (int, float)) else "—"
        lines.append(f"| {row['scope']} | {row['arm']} | {row['status']} | {m4_text} | {m5_text} |")
    lines.extend([
        "",
        "AUROC rows are descriptive sensitivity results and do not separate causal disease signal from batch. The transductive arm uses held-out expression for correction and is not external inductive performance. Partial correction retains small/unseen-batch samples at original ranks.",
        "",
        f"Frozen split manifest SHA-256: `{split_hash}`.",
    ])
    return "\n".join(lines) + "\n"


def recover_completed_scope(scope_dir: Path, arm: str, scope: str) -> dict | None:
    """Rebuild a result only when the runner's terminal completion artifacts exist."""
    not_evaluable_path = scope_dir / "NOT_EVALUABLE.json"
    if not_evaluable_path.is_file():
        row = json.loads(not_evaluable_path.read_text())
        if row.get("status") == "not_evaluable" and row.get("scope") == scope and row.get("arm") == arm:
            return row
        raise ValueError(f"Invalid NOT_EVALUABLE record in {scope_dir}")

    coverage_path = scope_dir / "combat_coverage.json"
    if not coverage_path.is_file():
        return None
    coverage = json.loads(coverage_path.read_text())
    if arm == "transductive_no_disease" and coverage.get("arm") != arm:
        raise ValueError(f"Coverage arm mismatch in {scope_dir}")
    if coverage.get("scope") not in (None, scope):
        raise ValueError(f"Coverage scope mismatch in {scope_dir}")
    predictions_path = scope_dir / "heldout_predictions.npz"
    metrics = {}
    if predictions_path.is_file():
        with np.load(predictions_path) as saved:
            y = saved["y"]
            metrics = {name: gp.metrics(y, saved[name])
                       for name in ("M1", "M2", "M3", "M4", "M5")
                       if name in saved.files}
    elif scope != "pooled/final":
        return None
    return {
        "status": "complete",
        "scope": scope,
        "arm": arm,
        "interpretation": (
            "transductive all-expression partial ComBat; held-out expression informs correction"
            if arm == "transductive_no_disease"
            else "inductive train-only partial ComBat sensitivity"
        ),
        "seconds": None,
        "metrics": metrics,
        "coverage_file": "combat_coverage.json",
        "model_fit_count": coverage.get("fit_count", coverage.get("model_fit_count")),
        "PPMI_accessed": False,
        "recovered_after_interruption": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=OUT,
        help="new empty output directory; use a new path to preserve an interrupted run",
    )
    parser.add_argument(
        "--resume-from", type=Path,
        help="reuse completed arm-scopes from an interrupted run, leaving its files unchanged",
    )
    args = parser.parse_args()
    output_root = args.output.expanduser().resolve()
    if version("neuroCombat") != PINNED_NEUROCOMBAT:
        raise ValueError(f"Expected neuroCombat=={PINNED_NEUROCOMBAT}")
    if output_root.exists() and any(output_root.iterdir()):
        raise ValueError(f"Refusing to reuse non-empty output directory: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)

    data, input_hashes = gp.load_geo()
    splits = gp.split_manifest(data, input_hashes)
    split_payload = json.dumps(splits, indent=2, ensure_ascii=False, allow_nan=False).encode()
    split_hash = hashlib.sha256(split_payload).hexdigest()
    if split_hash != EXPECTED_SPLIT_SHA256:
        raise ValueError(f"Frozen split hash mismatch: {split_hash}")
    batch_lookup, batch_hashes = load_batch_lookup()
    if set(batch_lookup) != set(map(str, data.sample)):
        raise ValueError("Batch metadata does not exactly align with GEO samples")

    scopes = sorted(splits["splits"], key=lambda item: item["scope"])
    sample_lookup = {str(sample): i for i, sample in enumerate(data.sample)}
    started = time.monotonic()
    summary_rows: list[dict] = []
    completed_keys: set[tuple[str, str]] = set()
    if args.resume_from:
        previous_root = args.resume_from.expanduser().resolve()
        previous_progress_path = previous_root / "PROGRESS.json"
        previous_started_path = previous_root / "STARTED.json"
        if not previous_progress_path.is_file() or not previous_started_path.is_file():
            raise ValueError(f"Cannot resume: missing STARTED.json or PROGRESS.json in {previous_root}")
        previous_started = json.loads(previous_started_path.read_text())
        previous_progress = json.loads(previous_progress_path.read_text())
        for field, expected in (("split_manifest_sha256", split_hash),
                                ("input_sha256", input_hashes),
                                ("batch_metadata_sha256", batch_hashes),
                                ("neuroCombat", PINNED_NEUROCOMBAT)):
            if previous_started.get(field) != expected:
                raise ValueError(f"Cannot resume: prior {field} does not match current inputs")
        prior_count = int(previous_progress.get("completed_arm_scopes", 0))
        results_path = previous_root / "RESULTS.jsonl"
        if results_path.is_file():
            summary_rows = [json.loads(line) for line in results_path.read_text().splitlines() if line]
        rows_by_key = {(str(row["arm"]), str(row["scope"])): row for row in summary_rows}
        # The original runner retained only the last row in PROGRESS.json.
        # Recover every completed scope from its final coverage/prediction
        # artifacts; an interrupted fit log alone is never treated as complete.
        for arm in ARMS:
            arm_root = previous_root / arm
            if not arm_root.is_dir():
                continue
            for old_scope_dir in arm_root.iterdir():
                if not old_scope_dir.is_dir():
                    continue
                scope = old_scope_dir.name.replace("__", "/")
                if scope not in {item["scope"] for item in scopes}:
                    raise ValueError(f"Unexpected prior scope directory: {old_scope_dir}")
                key = (arm, scope)
                recovered = recover_completed_scope(old_scope_dir, arm, scope)
                if recovered is not None:
                    rows_by_key.setdefault(key, recovered)
        if len(rows_by_key) < prior_count:
            raise ValueError("Cannot safely resume: fewer complete artifact sets than PROGRESS records")
        valid_keys = {(arm, item["scope"]) for item in scopes for arm in ARMS}
        if not set(rows_by_key).issubset(valid_keys):
            raise ValueError("Cannot safely resume: result history contains an unknown arm or scope")
        summary_rows = [rows_by_key[(arm, item["scope"])]
                        for item in scopes for arm in ARMS
                        if (arm, item["scope"]) in rows_by_key]
        for row in summary_rows:
            key = (str(row["arm"]), str(row["scope"]))
            old_scope_dir = previous_root / key[0] / key[1].replace("/", "__")
            if row.get("status") not in ("complete", "not_evaluable"):
                raise ValueError(f"Cannot resume scope with unexpected status: {row.get('status')}")
            artifact_row = recover_completed_scope(old_scope_dir, key[0], key[1])
            if artifact_row is None or artifact_row.get("status") != row.get("status"):
                raise ValueError(f"Cannot verify prior completion artifacts: {old_scope_dir}")
            new_scope_dir = output_root / key[0] / key[1].replace("/", "__")
            new_scope_dir.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(old_scope_dir, new_scope_dir)
            completed_keys.add(key)
        if summary_rows:
            (output_root / "RESULTS.jsonl").write_text(
                "".join(json.dumps(row, ensure_ascii=False, allow_nan=False) + "\n"
                        for row in summary_rows)
            )
    progress = Progress(len(scopes) * len(ARMS), completed_units=len(summary_rows))
    for row in summary_rows:
        print(f"Reusing completed {row['arm']} | {row['scope']}", flush=True)
    write_json(output_root / "STARTED.json", {
        "status": "running",
        "scope_count": len(scopes),
        "arm_count": len(ARMS),
        "arms": list(ARMS),
        "split_manifest_sha256": split_hash,
        "input_sha256": input_hashes,
        "batch_metadata_sha256": batch_hashes,
        "neuroCombat": PINNED_NEUROCOMBAT,
        "numpy_required": "1.26.4",
        "PPMI_accessed": False,
        "real_data_fits_started": True,
    })

    for scope_index, split in enumerate(scopes, start=1):
        scope = split["scope"]
        train, test, inner = scope_data(split, data, sample_lookup)
        designated_ids = list(split["train"]) + list(split["test"])
        if len(set(designated_ids)) != len(designated_ids):
            raise ValueError(f"Duplicate sample in frozen scope {scope}")
        full_scope = data.take(np.asarray([sample_lookup[s] for s in designated_ids], dtype=int))
        for arm in ARMS:
            if (arm, scope) in completed_keys:
                continue
            label = f"scope {scope_index}/{len(scopes)} {arm} {scope}"
            progress.update(label, force=True)
            arm_output = output_root / arm / scope.replace("/", "__")
            arm_output.mkdir(parents=True, exist_ok=False)
            scope_started = time.monotonic()
            try:
                if arm == "transductive_no_disease":
                    result = run_transductive(scope, train, test, full_scope,
                                              batch_lookup, inner, arm_output, progress)
                else:
                    result = run_inductive(scope, arm, train, test, inner,
                                           batch_lookup, arm_output, progress)
                result["seconds"] = time.monotonic() - scope_started
            except NotEvaluable as exc:
                result = {
                    "status": "not_evaluable",
                    "scope": scope,
                    "arm": arm,
                    "reason": str(exc),
                    "PPMI_accessed": False,
                    "metrics": {},
                }
                write_json(arm_output / "NOT_EVALUABLE.json", result)
            except Exception as exc:
                failure = {
                    "status": "failed",
                    "scope": scope,
                    "arm": arm,
                    "error_type": type(exc).__name__,
                    "message": str(exc),
                    "PPMI_accessed": False,
                }
                write_json(arm_output / "FAILURE.json", failure)
                write_json(output_root / "FAILURE.json", {
                    **failure,
                    "completed_arm_scopes": len(summary_rows),
                    "total_arm_scopes": len(scopes) * len(ARMS),
                })
                raise
            summary_rows.append(result)
            with (output_root / "RESULTS.jsonl").open("a", encoding="utf-8") as stream:
                stream.write(json.dumps(result, ensure_ascii=False, allow_nan=False) + "\n")
            write_json(output_root / "PROGRESS.json", {
                "completed_arm_scopes": len(summary_rows),
                "total_arm_scopes": len(scopes) * len(ARMS),
                "latest": result,
            })
            progress.complete_unit(f"completed {len(summary_rows)}/{len(scopes) * len(ARMS)} {arm} {scope}")

    summary = {
        "status": "complete",
        "analysis": "registered partial-ComBat sensitivity on frozen GEO folds",
        "arms": list(ARMS),
        "scope_count": len(scopes),
        "arm_scope_count": len(summary_rows),
        "elapsed_minutes": (time.monotonic() - started) / 60,
        "results": summary_rows,
        "split_manifest_sha256": split_hash,
        "input_sha256": input_hashes,
        "batch_metadata_sha256": batch_hashes,
        "neuroCombat": PINNED_NEUROCOMBAT,
        "minimum_training_batch_n": MIN_TRAIN_BATCH_N,
        "PPMI_accessed": False,
        "pooled_final_has_no_outer_test_metric": True,
        "interpretation": "Exploratory sensitivity; not causal separation of disease and batch. Transductive results use held-out expression for correction.",
    }
    write_json(output_root / "SUMMARY.json", summary)
    (output_root / "REPORT.md").write_text(markdown_report(summary_rows, split_hash))
    print(f"GEO_COMBAT_SENSITIVITY_COMPLETE: {output_root}", flush=True)


if __name__ == "__main__":
    main()
