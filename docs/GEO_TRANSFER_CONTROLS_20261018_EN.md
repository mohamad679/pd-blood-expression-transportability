# GEO sex and age transfer controls

**Run date:** 2026-10-18. **Population:** healthy controls in the two development
sources only. **PPMI access or scoring:** none.

The registered development sex gate passed in both directions. The age control
did not demonstrate the registered development capacity in either direction.
The latter limits interpretation; it does not authorize a new transformation,
feature panel, tuning grid, or endpoint. These controls do not establish
microarray-to-RNA-seq transfer or resolve PPMI sample linkage.

The [machine-readable aggregate](../results/geo_controls_20261018/aggregate_summary.json)
contains exact metrics, candidate scores, bootstrap intervals and hashes. The
[execution plan](GEO_CONTROL_EXECUTION_PLAN_20261018_EN.md) and code were committed
as `36db979` before these real control fits. The plan's original bytes and hash
are preserved. This private Git implementation record is not described as a new
public preregistration. The controls were specified in the existing public v3.1
protocol; operational details of this implementation are disclosed in the plan.

## Participants and separation

| Training source | Training controls | Female / male | External test source | Test controls | Female / male |
|---|---:|---:|---|---:|---:|
| GSE57475 / PROBE | 49 | 23 / 26 | GSE99039 / GENEPARK-TRO | 119 | 85 / 34 |
| GSE99039 / GENEPARK-TRO | 119 | 85 / 34 | GSE57475 / PROBE | 49 | 23 / 26 |

Only `case=0` records were used. No PD case entered preprocessing, tuning,
training or evaluation for these controls. Each operational donor has exactly
one included record; no donor spans the training and test sources or an inner
training/validation pair. PROBE donor identity retains the original limitation:
GSM is an operational donor identifier where a separate biological subject ID
was unavailable. The checks cannot prove absence of undocumented donor overlap.

The four prepared input hashes match the frozen input audit. The core disease
pipeline remains byte-for-byte unchanged, SHA-256
`6f7f6cf42608002b185a63044387c9b9f5e624fc0b9ff5ae156030e3655e1a2c`.
The private auxiliary split manifest hash is
`2e97f70f530824069cf36e9c8b4f1bbeef3d3d8513ffb0a6bc13b67f52682e44`.
Actual memberships, individual predictions, preprocessing and fitted arrays
remain in ignored `private_data/geo_controls_20261018/`.

## Sex: development gate passed

| Direction | AUROC | 95% bootstrap CI | BA at 0.5, descriptive | Development AUROC >= 0.95 |
|---|---:|---|---:|---|
| PROBE → GENEPARK-TRO | 1.000000 | [1.000000, 1.000000] | 0.500000 | Pass |
| GENEPARK-TRO → PROBE | 0.993311 | [0.974875, 1.000000] | 0.500000 | Pass |

The model is glum 3.4.1 binomial/logit ridge regression, trained only on the
five registered Y markers: KDM5D, DDX3Y, RPS4Y1, UTY, EIF1AY. XIST remains
absent; no marker was added or dropped based on these results. The intercept is
unpenalized, all five marker slopes are penalized, and auxiliary weights are
uniform. Feature means/SDs come only from the respective training controls.

Five sex-stratified donor folds selected lambda by mean validation AUROC, with
the registered larger-lambda tie preference. In each source the initial grid
was extended upward twice, to 10 and then 100. All candidate means were tied
within tolerance: 0.991667 in PROBE and 1.000000 in GENEPARK-TRO. Lambda=100 was
therefore selected in both. Both selections remain at the grid edge after the
two permitted extensions; no further search was performed. Large shrinkage
and differences in source sex proportions make fixed-threshold classification
uninformative here. The observed BA=0.5 is retained, not repaired by threshold
selection. The gate was explicitly based on AUROC, not BA, before this run.

Across 72 checked sex fits, every fit passed the fixed KKT <= 1e-6 rule; the
largest residual was 2.6747e-8. There were no captured fit warnings, prediction
warnings, or retries. Transfer AUROC using linear scores exactly matched that
using glum probabilities in both directions. These auxiliary fits did not
reproduce the disease-tuning sigmoid issue.

The 5,000 sex-stratified donor bootstraps condition on the developed models and
observed test records. A [1,1] empirical interval when all observed pairs are
correctly ranked does **not** establish perfect population performance or
include training uncertainty. The protocol's development gate uses the two
point estimates. The separate PPMI criterion, AUROC >= 0.95 and lower CI bound
>= 0.90 in PPMI controls, has **not** been evaluated.

## Age: development capacity not demonstrated

| Direction | Spearman | 95% bootstrap CI | MAE, years | MAE of training-median predictor | MAE difference, years |
|---|---:|---|---:|---:|---:|
| PROBE → GENEPARK-TRO | -0.058118 | [-0.238535, 0.121006] | 7.038111 | 7.025210 | +0.012901 |
| GENEPARK-TRO → PROBE | 0.258599 | [-0.051083, 0.522189] | 17.482476 | 7.693878 | +9.788598 |

For the age-MAE difference against the same training-median comparator, paired
bootstrap CIs were [-0.542604, 0.569721] and [6.680570, 12.797813] years,
respectively. The absolute MAE CIs were [5.940462, 8.189174] and
[15.235214, 19.735341] years. No Spearman bootstrap replicate was undefined.

The control used 1,000 genes selected independently in each training partition
from the fixed 15,720-gene autosomal universe by positive MAD, with gene-ID
tie breaking. Means and SDs were learned inside training partitions only.
The estimator was sklearn Ridge with an unpenalized intercept and the SVD
solver; alpha penalizes slopes under a summed squared-error objective. Five
donor folds selected alpha from the registered 0.01, 0.1, 1, 10, 100 grid by
mean validation MAE. The selected alpha was 100 in both sources, with inner
MAE 8.416929 and 6.575356 years. It is the upper grid endpoint; the age protocol
did not authorize an adaptive alpha extension, so none was performed.

Training-median ages were 59 years for PROBE and 61 years for GENEPARK-TRO.
Observed control age ranges were 47–87 and 35–82 years, respectively. The
registered capacity rule required Spearman >= 0.4 **and** MAE better than the
constant comparator in **both** directions. Neither direction satisfied it.
The results show poor transfer of this fixed age estimator. They do not isolate
platform, age distribution, sample size or biology as its cause. They neither
establish nor refute transferable PD gene information.

## Consequences and unchanged boundaries

1. The development prerequisite for the sex gate is established. PPMI execution
   is still conditional on final linkage/eligibility, remaining diagnostics,
   the second public lock, and the separate PPMI sex gate.
2. This age model has insufficient development support. A later poor PPMI age
   result cannot be presented as evidence of a PPMI-specific processing defect.
   Age is not a stopping gate and cannot be used to retune the holdout method.
3. Y-marker success is a limited positive control; it does not certify general
   autosomal transfer, cell-mixture adjustment, biological donor identity or PD
   discrimination. CBC support is still outstanding.
4. No phenotype was overwritten, no discordant sample was excluded, and no
   disease model, feature transformation, threshold or endpoint was changed.
5. No new Zenodo version is created for this auxiliary step. Include these
   outcomes in the next scheduled public milestone and final pre-PPMI lock.

## Verification and reproduction on macOS

The [synthetic verification](../results/geo_controls_20261018/synthetic_verification.json)
passed eight checks before the real controls: held-out mutation invariance,
training-only standardization, donor separation, sex-target selection rather
than PD labels, recoverable synthetic sex/age signals, KKT acceptance, detection
of sigmoid-induced mixed-class tail ties, and paired age bootstrap behavior.
The study run finished in 66.12 seconds using the pinned Python 3.11 environment
with eight OpenMP threads. The parallel retrospective diagnostic is a separate
workload; the timing is not a general compute estimate.

From the repository root, with the frozen prepared inputs present, use **fresh**
private output and aggregate paths; do not overwrite the completed run:

```bash
OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 \
UV_CACHE_DIR=/private/tmp/pdgeo_uv_cache \
uv run --offline \
  --default-index https://pypi.tuna.tsinghua.edu.cn/simple \
  --with-requirements requirements_glum_v3.lock.txt \
  python pipeline/geo_controls.py \
  --output private_data/geo_controls_NEW_RUN \
  --aggregate /private/tmp/geo_controls_NEW_RUN.json
```

The offline command assumes dependencies are cached. Individual-level outputs
remain local. Only this report, code, the execution plan, synthetic verification
and aggregate results belong in Git. The completed control code hash is
`3a4c769b3445172452df41ce1556b46a6a0d87fd8b47fef52d58f23fd661f850`.
