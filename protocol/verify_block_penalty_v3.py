"""Independent synthetic verification; never loads GEO or PPMI participants."""

import os

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[variable] = "1"

import importlib.metadata
import json
import platform
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

from block_logistic_v3 import fit_block_logistic


def fit_reference(B, y, offset=None):
    """Unpenalized binomial likelihood optimized independently by SciPy."""
    design = np.column_stack((np.ones(len(B)), B))
    offset = np.zeros(len(y)) if offset is None else offset

    def objective(beta):
        eta = design @ beta + offset
        return np.mean(np.logaddexp(0, eta) - y * eta)

    def gradient(beta):
        return design.T @ (expit(design @ beta + offset) - y) / len(y)

    result = minimize(objective, np.zeros(design.shape[1]), jac=gradient,
                      method="BFGS", options={"gtol": 1e-10, "maxiter": 5000})
    if np.max(np.abs(gradient(result.x))) > 1e-7:
        raise AssertionError("Independent reference likelihood did not converge")
    return result.x


def main():
    rng = np.random.default_rng(20261015)
    n, k, p = 1200, 4, 24
    B = rng.normal(size=(n, k))
    G = rng.normal(size=(n, p))
    G[:, 0] += 0.5 * B[:, 0]
    X = np.column_stack((B, G))
    X = (X - X.mean(axis=0)) / X.std(axis=0)
    B, G = X[:, :k], X[:, k:]
    eta = -0.2 + B @ np.array([0.9, -0.7, 0.5, 0.25]) + 0.6 * G[:, 0]
    y = rng.binomial(1, expit(eta)).astype(float)
    ref = fit_reference(B, y)
    checks = []
    caught = []
    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always")
        # With no gene block, alpha must not change the unpenalized fit.
        for ratio in (0.0, 0.5, 1.0):
            for alpha in (0.0, 0.01, 1.0, 100.0):
                model = fit_block_logistic(B, y, k, alpha, ratio)
                error = np.max(np.abs(np.r_[model.intercept_, model.coef_] - ref))
                assert error < 2e-5, ("baseline invariance", alpha, ratio, error)
                checks.append({"test": "baseline_only_matches_unpenalized_reference",
                               "alpha": alpha, "l1_ratio": ratio, "max_coefficient_error": float(error)})
        # Baseline coefficients may move when fitted gene effects move. Their
        # loss gradient, not their numerical invariance, verifies no penalty.
        for ratio in (0.0, 0.5, 1.0):
            for alpha in (0.01, 0.2, 1.0):
                model = fit_block_logistic(X, y, k, alpha, ratio)
                pred = model.predict(X)
                grad = X.T @ (pred - y) / n
                base_residual = max(abs(float(np.mean(pred-y))), np.max(np.abs(grad[:k])))
                beta = model.coef_[k:]
                active = np.abs(beta) > 1e-8
                residual = grad[k:] + alpha * (1-ratio) * beta
                active_kkt = np.max(np.abs(residual[active] + alpha * ratio * np.sign(beta[active]))) if active.any() else 0.0
                inactive_kkt = max(0.0, float(np.max(np.abs(residual[~active])) - alpha * ratio)) if (~active).any() else 0.0
                conditional_ref = fit_reference(B, y, offset=G @ beta)
                conditional_error = np.max(np.abs(np.r_[model.intercept_, model.coef_[:k]] - conditional_ref))
                assert base_residual < 2e-6, ("baseline gradient", base_residual)
                assert max(active_kkt, inactive_kkt) < 2e-6, ("gene KKT", active_kkt, inactive_kkt)
                assert conditional_error < 2e-5, ("conditional reference", conditional_error)
                checks.append({"test": "partial_penalty_KKT_and_independent_conditional_reference",
                               "alpha": alpha, "l1_ratio": ratio,
                               "baseline_loss_gradient_max": float(base_residual),
                               "gene_KKT_violation_max": float(max(active_kkt, inactive_kkt)),
                               "conditional_baseline_coefficient_error": float(conditional_error),
                               "nonzero_genes": int(active.sum())})
        # Strong L1 penalty sets all genes to zero and recovers M4 exactly.
        strong = fit_block_logistic(X, y, k, 1.0, 1.0)
        strong_error = np.max(np.abs(np.r_[strong.intercept_, strong.coef_[:k]] - ref))
        assert np.max(np.abs(strong.coef_[k:])) < 1e-9
        assert strong_error < 2e-5
        # Negative control: deliberately penalize every column. The same
        # unpenalized-baseline gradient check must detect this wrong mask.
        wrong = fit_block_logistic(X, y, 0, 0.2, 0.5)
        wrong_grad = np.max(np.abs(B.T @ (wrong.predict(X)-y) / n))
        assert wrong_grad > 1e-3, "Verification failed to detect a deliberately wrong mask"
        caught = [str(w.message) for w in recorded]
    assert not caught, caught
    packages = {dist.metadata["Name"]: dist.version for dist in importlib.metadata.distributions()}
    payload = {"status": "passed", "data": "synthetic_only", "seed": 20261015,
               "n": n, "baseline_columns": k, "gene_columns": p,
               "python": platform.python_version(), "packages": dict(sorted(packages.items())),
               "checks": checks, "strong_L1_baseline_error": float(strong_error),
               "deliberately_wrong_mask_baseline_gradient": float(wrong_grad),
               "warnings": caught, "GEO_or_PPMI_models_fitted": False,
               "note": "Unpenalized means zero penalty derivative, not coefficient invariance when gene effects vary."}
    path = Path(__file__).with_name("block_penalty_v3_verification.json")
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    requirements = [f"{name}=={version}" for name, version in sorted(packages.items())
                    if name.lower() not in {"pip", "setuptools", "wheel"}]
    Path(__file__).with_name("requirements_glum_v3.lock.txt").write_text("\n".join(requirements)+"\n")
    print(json.dumps({"status": "passed", "checks": len(checks)+2,
                      "glum": importlib.metadata.version("glum"),
                      "strong_L1_baseline_error": float(strong_error),
                      "wrong_mask_detected": True, "real_models_fitted": False}))


if __name__ == "__main__":
    main()
