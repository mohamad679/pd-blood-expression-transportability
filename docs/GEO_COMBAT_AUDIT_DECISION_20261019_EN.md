# GEO ComBat audit decision (2026-10-19)

## Scope and estimand

ComBat remains a bounded, exploratory development-data sensitivity, not part of the primary pipeline. It assesses how conclusions change under alternative handling of the recorded GEO technical batches; it does not establish causal separation of batch and disease. Use only the two GEO development cohorts and the already frozen outer and inner partitions. PPMI is excluded.

The actual batch variable is cohort-specific: the supplied batch field in GSE99039 (GENEPARK/TRO), and the title-derived chip proxy in GSE57475 (PROBE). These variables are not assumed to be equivalent technical factors across cohorts. Technical metadata never becomes a predictor in M1–M5.

## Prespecified contrasts

1. **Training-only, without disease covariate (inductive):** fit ComBat parameters using expression and batch labels from the current training partition only, with no disease covariate. Apply those saved parameters to its held-out partition. For nested CV, repeat from scratch inside each inner training partition and apply to its validation partition; refit on the full outer-training partition before the outer test. No held-out expression or diagnosis is used to estimate correction parameters.
2. **Training-only, preserving disease (inductive):** same partitioning and application rule, with the binary disease indicator included as a biological covariate only for samples in the current training partition. It is never supplied for validation or test samples. This is a distinct sensitivity because batch and disease imbalance can make their effects weakly identifiable.
3. **All-expression, without disease covariate (transductive, registered diagnostic):** estimate unsupervised batch correction from the expression and batch metadata of all samples in the relevant development evaluation scope, without disease labels, then evaluate with the frozen folds. This deliberately uses held-out expression to estimate correction and is therefore transductive; it must be labeled and must not be described as inductive external performance.

There is **no all-expression ComBat arm with disease as a covariate**. Doing so would require supplying held-out disease labels to ComBat, which is prohibited. An unknown/test label must never be imputed, inferred, or substituted to make the covariate matrix complete. Thus the registered comparison is not a full factorial design: the label-preserving comparison is train-only, while the all-expression contrast is unsupervised and explicitly transductive.

Within every train-only arm, all ComBat parameters and every downstream data-dependent step are re-estimated at the relevant fit boundary. Inner validation/test labels are used only for the already specified evaluation metrics. Keep original fixed folds, model definitions, tuning grid, convergence criteria, and reporting rules unchanged. A previously unseen held-out batch cannot receive parameters learned from held-out samples in an inductive arm; if the pinned implementation cannot transform it from training estimates, report that partition as not evaluable for that arm rather than estimating from test data.

## Operational rules added before execution (2026-10-20)

### Input and locked implementation

- Input to ComBat is the existing **within-sample gene-rank matrix**, not raw or normalized intensity. The corrected values therefore remain adjusted rank-scale values; ComBat's usual continuous-expression interpretation does not automatically carry over. This is a limitation to state with the results.
- The diagnostic implementation is pinned to `neuroCombat==0.2.12` in the dedicated [ComBat requirements lock](../requirements_glum_v3_combat.lock.txt). The primary-analysis lock is left unchanged. Use the package defaults explicitly (`eb=True`, `parametric=True`, `mean_only=False`, `ref_batch=None`) unless the synthetic gate identifies a software incompatibility; any method change requires a new decision before real fits.
- Apply the correction before the existing training-only source-stratified MAD selection. Recompute MAD selection and all model fitting/tuning from each training partition. Do not choose among ComBat arms by which AUROC is largest.

### Small-batch and design-rank rules

- Evaluate batch counts separately in every actual ComBat fitting partition, for each cohort, arm, outer fold, and inner fold. Only batches with at least **3 training samples** contribute to ComBat parameter estimation. A batch below this threshold is omitted from the ComBat fit; its training samples remain in the model's training data on their original within-sample ranks, without batch correction. Do not merge batches or remove samples from the analysis population.
- In train-only arms, apply frozen training estimates only to held-out samples from batches represented by at least 3 samples in that training fit. Samples from a batch below threshold or absent in training remain on their original within-sample ranks, without correction. Never use held-out expression or diagnosis to estimate a batch parameter. Record corrected and uncorrected sample counts separately for training and evaluation data in every partition.
- In the transductive arm, apply the same minimum of 3 to the full unlabeled expression set used to estimate correction. Samples in smaller batches remain uncorrected. Label this arm transductive in all tables and prose.
- For the disease-preserving arm, also check the training design matrix for rank deficiency/complete batch-disease confounding using only batches included in the ComBat fit. If the disease coefficient is not identifiable in that partition, report the arm as non-evaluable there; do not change covariates or batches after seeing outcomes.
- Apply these rules to both cohorts. This is a **partial ComBat** sensitivity: large batches are corrected and small/unseen batches are retained uncorrected. Report correction coverage and this mixed-processing limitation; do not describe it as complete batch removal.

### Metadata-only eligibility decision (2026-10-20)

The frozen split manifest was reconstructed from sample metadata only and its SHA-256 matched the registered value `2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2`. No expression values were read. Under the initial rule that any training batch below 3 made the whole partition non-evaluable:

| Cohort | Outer train/test partitions invalid | Inner train/validation fits invalid |
|---|---:|---:|
| GSE57475 (PROBE) | 0/5 | 8/25 |
| GSE99039 (GENEPARK) | 5/5 | 24/25 |

The initial rule would remove every GENEPARK outer fold and nearly every GENEPARK inner fit. Metadata also showed 16/16 PROBE chip batches contain at least 3 samples (142/142 samples), while 10/11 GENEPARK batches do (238/239 samples); the only full-cohort batch below threshold is the singleton `08/27/09_TRO`. Each outer training fold in both cohorts still has at least two batches meeting the threshold. Before looking at any ComBat result, the whole-partition exclusion rule is therefore amended to the partial-correction rule above: omit only small batches from ComBat estimation, retain their samples uncorrected, and report correction coverage. The complete metadata-only tally is in [GEO_COMBAT_BATCH_FEASIBILITY_20261020_EN.md](GEO_COMBAT_BATCH_FEASIBILITY_20261020_EN.md).

### Software feasibility gate (must pass before real fits)

The upstream `neuroCombatFromTraining` API is described by its maintainers as **under development** and accepts only the new data, batch labels, and training estimates; the documented interface does not accept new-sample biological covariates. Therefore, pinning a version is not by itself evidence that every registered inductive arm is implementable. Before any real fit, add and pass synthetic tests showing that (a) train-only estimates are applied to held-out samples without using held-out expression to re-estimate ComBat parameters, (b) held-out diagnosis is never read by the transform, (c) known batches transform correctly and unseen batches fail closed, and (d) the disease-preserving arm retains the specified training covariate contribution without needing held-out labels. If the pinned API cannot satisfy these conditions, stop before real-data fits and report the limitation for a new decision. Do not silently substitute an all-expression/transductive fit for an inductive arm.

The same synthetic gate must verify deterministic handling of the rank matrix, the minimum-count checks, and a known-batch train/apply example against an independently computed reference. No ComBat fit or GEO result has been run by adding these rules.

## Interpretation cautions

When disease composition differs by batch, omitting disease from the adjustment may attribute genuine disease-associated expression to batch and remove part of it. Including disease as a covariate can instead preserve or accentuate estimated group differences under confounding and unstable identification. Nygaard et al. describe bias and inflated confidence in differential-expression inference when group and batch are imbalanced; this supports caution about both choices, not a claim that a particular direction or magnitude of bias must occur in this project. Report the arms as sensitivity analyses and avoid selecting one based on the more favorable AUROC.

Reference: Nygaard, Rødland & Hovig (2016), “Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses,” *Biostatistics* 17(1), 29–39. https://doi.org/10.1093/biostatistics/kxv027

Software references: [neuroCombat 0.2.12 on PyPI](https://pypi.org/project/neuroCombat/0.2.12/) and [official neuroCombat API/source](https://github.com/Jfortin1/neuroCombat). The upstream source marks `neuroCombatFromTraining` as under development and documents that test batches must occur in training; this is why the synthetic feasibility gate above is mandatory.

## Implementation clarification recorded after the run (2026-10-20)

The operational rule above said ComBat precedes downstream data-dependent steps, but did not name the panel feature path explicitly. In the implemented sensitivity, M4 panel scores are recomputed from the ComBat-adjusted expression/rank matrix for samples in eligible batches; samples in small or unseen batches retain the original rank matrix. M5 receives the same adjusted panel features plus its separately selected gene block. This paragraph clarifies observed code behavior **after outcomes were available**; it is not a new prespecified analysis rule. The original M1–M5 pipeline and its stored results were not changed. See the separate [execution audit](GEO_COMBAT_EXECUTION_AUDIT_20261020_EN.md) for the timing and numerical-warning record.
