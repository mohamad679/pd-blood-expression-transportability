# Numerical amendment before locking — 16 October 2026

> **Provenance and status.** This English text was prepared after the first real PROBE fold began at 2026-10-15 22:07:32 UTC. Its source is the timestamped Persian v3.1 source (https://doi.org/10.5281/zenodo.23405371): `KKT_NUMERICAL_AMENDMENT_20261016_FA.md` (SHA-256 `ffad1329e0493ddd2f83536a77f3726901eabe910909ad23020e71dcae0ebe5f`). It introduces no method change. If wording conflicts, the original Persian preregistration and its published amendment govern. Statements of registration, execution, and publication status in the translated body describe its original snapshot, not the reader's current date. See [TRANSLATION_SOURCE_MAP.json](../TRANSLATION_SOURCE_MAP.json).

This amendment is not public registration and arises solely from synthetic-data outputs. Objective, features, weights, lambda/rho, KKT=1e-6 threshold and selection criterion remain fixed. No GEO result or individual PPMI score informed the change.

## Documented rationale

All 48 initial scale/KKT comparisons with 40 and 1000 genes passed. However, in the strongly platform-correlated 1000-gene fixture, pure ridge with irls-cd and lambda=1 took about 31.19 seconds; lambda=0.1,rho=0 had produced no output after more than 244.3 seconds. The first run was deliberately interrupted; interruption is not proof of mathematical nonconvergence, and the candidate was neither removed nor assigned an unreliable AUC for use. The four completed fits and first code/plan/report versions are preserved in kkt_history_first_pass_20261016.

A diagnostic solve of the same problem with independent SciPy took 0.0933 seconds, giving KKT_max=2.3892e-8 and objective=0.3551254255474776. The L2-specific glum irls-ls solve took about 1.026 seconds, giving KKT_max=6.4817e-14 and objective=0.35512542554740933; the objective difference was about 6.8e-14, with no warning. This comparison was performed after the slowdown was observed and was not part of the initial plan; it diagnoses and justifies the numerical amendment rather than silently replacing the stress result.

## Bounded change specified before the rerun

- Use `solver=irls-ls` only for `rho==0 AND lambda>0`; it solves the same GLM with an exclusively L2 penalty. Other candidates and baseline-only remain irls-cd.
- All initial settings and the single cold retry, weights, P1/P2, normalized loss, KKT criterion and SciPy agreement thresholds remain fixed.
- If the new solver fails KKT, apply the same fixed retry; remaining failure is invalidity, not permission to remove the candidate.
- Repeat all 48 scale comparisons with the same seeds and the full stress grid from the beginning; new filenames use solver_v2. The previous scale comparison and partial stress outputs are not overwritten.
- The lambda/rho grid is neither contracted nor expanded in this amendment. If the new stress run has unresolved failures, report numbers before continuing. After passing, the same-seed 23/66/9 and null reruns are required as planned.

Public Zenodo registration, final v3.1 locking and PPMI scoring have still not occurred. This amendment only selects a suitable numerical method for the same mathematical objective and does not guarantee prediction results in advance.

Solver capability sources: [official glum algorithm documentation](https://glum.readthedocs.io/en/latest/background.html) and [API/scale_predictors settings](https://glum.readthedocs.io/en/latest/glm.html). The installed version throughout these runs is glum 3.4.1; equation verification relies on independent SciPy as well as documentation.
