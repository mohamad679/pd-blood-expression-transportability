"""Synthetic feasibility gate for the registered partial ComBat sensitivity."""
from __future__ import annotations

from importlib.metadata import version
import sys

import numpy as np
import pandas as pd
from neuroCombat import neuroCombat, neuroCombatFromTraining


PINNED_VERSION = "0.2.12"
MIN_TRAIN_BATCH_N = 3
PROGRESS_STEPS = 7


def show_progress(step: int, label: str) -> None:
    """Show synthetic-test progress without changing any test behavior."""
    width = 24
    filled = round(width * step / PROGRESS_STEPS)
    bar = "#" * filled + "-" * (width - filled)
    message = f"Progress [{bar}] {step}/{PROGRESS_STEPS} ({100 * step / PROGRESS_STEPS:5.1f}%) | {label}"
    if sys.stdout.isatty():
        print("\r" + message, end="", flush=True)
        if step == PROGRESS_STEPS:
            print(flush=True)
    else:
        print(message, flush=True)


def fit_partial_combat(x: np.ndarray, batch: np.ndarray, disease: np.ndarray | None):
    """Fit estimates using eligible training batches only; return original rows too."""
    counts = pd.Series(batch).value_counts().to_dict()
    eligible = {name for name, n in counts.items() if n >= MIN_TRAIN_BATCH_N}
    keep = np.asarray([name in eligible for name in batch])
    if len(eligible) < 2:
        raise ValueError("At least two eligible batches are required")
    covars = pd.DataFrame({"batch": batch[keep]})
    categorical = []
    if disease is not None:
        covars["disease"] = disease[keep]
        categorical = ["disease"]
    fitted = neuroCombat(
        dat=x[keep].T,
        covars=covars,
        batch_col="batch",
        categorical_cols=categorical,
        eb=True,
        parametric=True,
        mean_only=False,
        ref_batch=None,
    )
    return {"eligible": eligible, "keep": keep, "estimates": fitted["estimates"]}


def transform_partial_combat(x: np.ndarray, batch: np.ndarray, state: dict) -> tuple[np.ndarray, np.ndarray]:
    """Apply frozen estimates only to eligible known batches; leave others unchanged."""
    eligible_batches = set(state["estimates"]["batches"].astype(str))
    transform = np.asarray([str(name) in eligible_batches for name in batch])
    output = x.copy()
    if np.any(transform):
        adjusted = neuroCombatFromTraining(
            dat=x[transform].T,
            batch=np.asarray(batch[transform], dtype=str),
            estimates=state["estimates"],
        )["data"]
        output[transform] = adjusted.T
    return output, transform


def manual_frozen_transform(x: np.ndarray, batch: np.ndarray, estimates: dict) -> np.ndarray:
    """Independent direct implementation of the documented frozen-estimate formula."""
    levels = np.asarray(estimates["batches"], dtype=str)
    batch_index = np.asarray([int(np.flatnonzero(levels == str(name))[0]) for name in batch])
    pooled = np.asarray(estimates["var.pooled"])
    stand = np.asarray(estimates["stand.mean"])[:, 0]
    mod = np.asarray(estimates["mod.mean"])
    center = stand + mod.mean(axis=1)
    standardized = (x.T - center[:, None]) / np.sqrt(pooled)
    gamma = np.asarray(estimates["gamma.star"])[batch_index, :].T
    delta = np.asarray(estimates["delta.star"])[batch_index, :].T
    return (((standardized - gamma) / np.sqrt(delta)) * np.sqrt(pooled) + center[:, None]).T


def make_rank_matrix(seed: int, batches: list[str], labels: np.ndarray, genes: int = 80) -> np.ndarray:
    rng = np.random.default_rng(seed)
    gene_axis = np.linspace(-1.0, 1.0, genes)
    batch_effects = {name: rng.normal(0, 0.45, genes) for name in sorted(set(batches))}
    raw = np.column_stack([
        rng.normal(0, 1, genes) + batch_effects[name] + labels[i] * 0.15 * gene_axis
        for i, name in enumerate(batches)
    ])
    ranks = np.argsort(np.argsort(raw, axis=0), axis=0) / (genes - 1)
    return ranks.T


def main() -> None:
    show_progress(0, "starting synthetic checks")
    installed = version("neuroCombat")
    assert installed == PINNED_VERSION, f"Expected neuroCombat {PINNED_VERSION}, got {installed}"
    show_progress(1, "pinned neuroCombat version verified")

    train_batch = np.repeat(["A", "B", "C", "SMALL"], [8, 8, 8, 2])
    train_y = np.tile([0, 1], 13)
    train_x = make_rank_matrix(41021, train_batch.tolist(), train_y)
    held_batch = np.asarray(["A", "B", "C", "SMALL", "UNSEEN"])
    held_y = np.asarray([1, 0, 1, 0, 1], dtype=int)  # deliberately never passed to transform
    held_x = make_rank_matrix(41022, held_batch.tolist(), held_y)
    show_progress(2, "synthetic train and held-out batches created")

    # Batch count rule: the two-row batch is omitted only from ComBat estimation.
    state = fit_partial_combat(train_x, train_batch, disease=None)
    assert state["eligible"] == {"A", "B", "C"}
    assert not state["keep"][train_batch == "SMALL"].any()
    adjusted_train, train_mask = transform_partial_combat(train_x, train_batch, state)
    assert np.array_equal(adjusted_train[~state["keep"]], train_x[~state["keep"]])
    assert np.all(train_mask[state["keep"]]) and not np.any(train_mask[~state["keep"]])
    show_progress(3, "small training batch excluded and unchanged")

    # Held-out small/unseen batches fail closed to unchanged ranks; known large batches transform.
    adjusted_held, held_mask = transform_partial_combat(held_x, held_batch, state)
    assert held_mask.tolist() == [True, True, True, False, False]
    assert np.array_equal(adjusted_held[~held_mask], held_x[~held_mask])
    estimates = state["estimates"]
    manual = manual_frozen_transform(held_x[held_mask], held_batch[held_mask], estimates)
    np.testing.assert_allclose(adjusted_held[held_mask], manual, rtol=1e-10, atol=1e-10)
    assert np.isfinite(adjusted_held).all()
    show_progress(4, "frozen held-out transform and formula verified")

    # Fit/apply has no held-out label argument: changing labels alone cannot alter the transform.
    adjusted_held_again, held_mask_again = transform_partial_combat(held_x, held_batch, state)
    np.testing.assert_array_equal(adjusted_held_again, adjusted_held)
    np.testing.assert_array_equal(held_mask_again, held_mask)
    show_progress(5, "held-out labels confirmed absent from transformation")

    # Disease-preserving train fit is permitted, but held-out disease is never passed to apply.
    disease_state = fit_partial_combat(train_x, train_batch, disease=train_y)
    disease_adjusted, disease_mask = transform_partial_combat(held_x, held_batch, disease_state)
    assert disease_mask.tolist() == held_mask.tolist()
    assert np.isfinite(disease_adjusted).all()
    assert np.array_equal(disease_adjusted[~disease_mask], held_x[~disease_mask])
    show_progress(6, "disease-preserving fit verified without held-out labels")

    # An unseen batch is left unchanged by the partial wrapper; package API itself fails closed.
    try:
        neuroCombatFromTraining(
            dat=held_x[-1:, :].T,
            batch=np.asarray(["UNSEEN"], dtype=str),
            estimates=estimates,
        )
    except ValueError:
        unseen_api_fails_closed = True
    else:
        unseen_api_fails_closed = False
    assert unseen_api_fails_closed
    show_progress(7, "unseen-batch fail-closed behavior verified")

    print("PASS pinned neuroCombat version:", installed)
    print("PASS <3-sample training batch omitted from estimation and left unchanged")
    print("PASS known eligible batches transformed from frozen training estimates")
    print("PASS frozen transformation matches independent formula")
    print("PASS held-out diagnosis is not an input to transform")
    print("PASS small and unseen held-out batches remain uncorrected")
    print("PASS disease-preserving training fit applies without held-out diagnosis")
    print("SYNTHETIC COMBAT FEASIBILITY PASSED")


if __name__ == "__main__":
    main()
