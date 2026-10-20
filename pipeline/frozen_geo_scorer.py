"""Inference-only GEO scorer; no fitting and no PPMI command-line entry point.

The module accepts a fixed within-sample rank matrix, age and recorded sex.
It never reads diagnosis labels or estimates any statistic from test samples.
The separate post-lock executor must validate the public lock before use.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.special import expit


MODELS = ("M1", "M2", "M3", "M4", "M5")
PANELS = (
    "T cells", "B lineage", "Monocytic lineage", "Neutrophils",
    "Erythroid", "Platelet",
)


class FrozenGEOScorer:
    def __init__(self, final_dir: Path, marker_manifest: Path, sex_artifact: Path):
        self.final_dir = Path(final_dir)
        with np.load(self.final_dir / "preprocessing.npz", allow_pickle=False) as z:
            self.age_mean = float(z["age_mean"])
            self.age_sd = float(z["age_sd"])
            self.base_mean = np.array(z["base_mean"], dtype=np.float64)
            self.base_sd = np.array(z["base_sd"], dtype=np.float64)
            self.keep = np.array(z["keep"], dtype=np.int64)
            self.selected = np.array(z["selected"], dtype=np.int64)
            self.gene_mean = np.array(z["gene_mean"], dtype=np.float64)
            self.gene_sd = np.array(z["gene_sd"], dtype=np.float64)
            self.genes = np.array(z["genes"])
        if self.age_sd <= 0 or np.any(self.gene_sd <= 0):
            raise ValueError("Frozen GEO preprocessing has invalid scales")
        if len(self.genes) != 16364 or len(set(self.genes.tolist())) != len(self.genes):
            raise ValueError("Frozen GEO gene universe is invalid")
        if len(self.base_mean) != 10 or len(self.base_sd) != 10:
            raise ValueError("Frozen GEO baseline schema differs from protocol")
        if len(self.selected) != 1000 or np.any(self.selected < 0) or np.any(self.selected >= len(self.genes)):
            raise ValueError("Frozen GEO selected-gene block is invalid")

        manifest = json.loads(Path(marker_manifest).read_text(encoding="utf-8"))
        lookup = {gene: i for i, gene in enumerate(self.genes)}
        self.marker_indices = tuple(
            np.array([lookup[gene] for gene in manifest["cell_marker_panels"][panel]["gene_ids"]], dtype=int)
            for panel in PANELS
        )
        self.decisions = json.loads((self.final_dir / "decisions.json").read_text(encoding="utf-8"))
        self.models = {}
        for name in MODELS:
            with np.load(self.final_dir / f"{name}_model.npz", allow_pickle=False) as z:
                intercept = float(z["intercept"])
                coef = np.array(z["coef"], dtype=np.float64)
            if not np.isfinite(intercept) or not np.isfinite(coef).all():
                raise ValueError("Nonfinite frozen GEO model")
            self.models[name] = (intercept, coef)

        with np.load(sex_artifact, allow_pickle=False) as z:
            self.sex_gene_ids = np.array(z["gene_ids"])
            self.sex_selected = np.array(z["selected"], dtype=np.int64)
            self.sex_mean = np.array(z["mean"], dtype=np.float64)
            self.sex_sd = np.array(z["sd"], dtype=np.float64)
            self.sex_intercept = float(z["intercept"])
            self.sex_coef = np.array(z["coef"], dtype=np.float64)
        if (len(self.sex_selected) != 5 or np.any(self.sex_sd <= 0)
                or not np.array_equal(self.genes[self.sex_selected], self.sex_gene_ids)):
            raise ValueError("Frozen GEO sex-marker artifact disagrees with gene universe")

    def validate_input(self, ranks: np.ndarray, age: np.ndarray, sex: np.ndarray, genes: np.ndarray) -> None:
        if not np.array_equal(np.asarray(genes), self.genes):
            raise ValueError("Input gene IDs/order differ from frozen GEO")
        if ranks.ndim != 2 or ranks.shape[1] != len(self.genes):
            raise ValueError("Input rank-matrix dimensions are invalid")
        if age.shape != (len(ranks),) or sex.shape != (len(ranks),):
            raise ValueError("Input age/sex dimensions are invalid")
        if (not np.isfinite(ranks).all() or ranks.min() <= 0 or ranks.max() >= 1
                or np.max(np.abs(ranks.mean(axis=1) - .5)) > 1e-12):
            raise ValueError("Input within-sample ranks violate GEO invariants")
        if not np.isfinite(age).all() or np.any(age <= 0) or np.any(age >= 120):
            raise ValueError("Input age is invalid")
        if not np.isin(sex, [0, 1]).all():
            raise ValueError("Input recorded sex must be 0/1")

    def design(self, ranks: np.ndarray, age: np.ndarray, sex: np.ndarray, model: str) -> np.ndarray:
        age_z = (age - self.age_mean) / self.age_sd
        raw = np.column_stack((
            age_z, age_z ** 2, sex, age_z * sex,
            *(ranks[:, ix].mean(axis=1) for ix in self.marker_indices),
        ))
        baseline_z = (raw - self.base_mean) / np.where(self.base_sd > 0, self.base_sd, 1.)
        keep = [int(i) for i in self.keep if i < 4] if model in ("M1", "M3") else self.keep
        baseline = baseline_z[:, keep] if model != "M2" else np.empty((len(ranks), 0))
        if model in ("M1", "M4"):
            return np.ascontiguousarray(baseline)
        genes = (ranks[:, self.selected] - self.gene_mean) / self.gene_sd
        return np.ascontiguousarray(np.column_stack((baseline, genes)))

    def predict_models(self, ranks: np.ndarray, age: np.ndarray, sex: np.ndarray, genes: np.ndarray) -> dict[str, np.ndarray]:
        ranks = np.asarray(ranks, dtype=np.float64)
        age = np.asarray(age, dtype=np.float64)
        sex = np.asarray(sex, dtype=np.float64)
        self.validate_input(ranks, age, sex, genes)
        output = {}
        for name in MODELS:
            actual = self.decisions.get(name, {}).get("selected", {}).get("reuses_baseline", name)
            if actual != name:
                if actual not in output:
                    raise ValueError("Frozen no-gene branch refers to unavailable baseline")
                output[name] = output[actual].copy()
                continue
            design = self.design(ranks, age, sex, name)
            intercept, coef = self.models[name]
            if design.shape[1] != len(coef):
                raise ValueError("Frozen model coefficient and design dimensions differ")
            output[name] = expit(intercept + design @ coef)
            if not np.isfinite(output[name]).all():
                raise ValueError("Nonfinite frozen GEO prediction")
        return output

    def predict_sex_control(self, ranks: np.ndarray, genes: np.ndarray) -> np.ndarray:
        ranks = np.asarray(ranks, dtype=np.float64)
        if ranks.ndim != 2 or not np.array_equal(np.asarray(genes), self.genes):
            raise ValueError("Sex-control input dimensions or gene order differ")
        design = (ranks[:, self.sex_selected] - self.sex_mean) / self.sex_sd
        result = expit(self.sex_intercept + design @ self.sex_coef)
        if not np.isfinite(result).all():
            raise ValueError("Nonfinite sex-control prediction")
        return result
