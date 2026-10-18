# GEO sex/age controls and sigmoid diagnostic: execution plan

Written on 2026-10-18, after the core GEO results and before running these
auxiliary controls. This implements the already registered controls in
`protocol/PROTOCOL_v3_1_EN.md`; it is not a new preregistration, a replacement
endpoint, or an authorization to score PPMI. No PPMI input is read.

## Fixed auxiliary implementation details

- Use the four frozen prepared GEO inputs, verified against the input audit.
  Use healthy controls only, once per operational donor. A duplicated donor
  causes failure rather than choosing a visit after examining outcomes.
- Run both directions, training on one GEO and evaluating on the other. All
  scaling, positive-MAD selection and tuning use the training source only.
- Use five shuffled donor folds. For sex, stratify by recorded sex; for age,
  use shuffled donor KFold without outcome binning. Use the existing scope-seed
  function, scopes `controls/{source}/sex/inner` and
  `controls/{source}/age/inner`. Save actual membership privately and its hash
  publicly. No PD case participates in either fit or validation.
- Auxiliary fits use uniform sample weights. The source-by-disease balancing
  rule for M1–M5 does not apply to a sex or age target in healthy controls.
- Sex: the five fixed Y markers, training mean/SD scaling, glum 3.4.1 ridge
  binomial logit, unpenalized intercept and all five slopes penalized. Use the
  registered lambda grid, mean inner-fold AUROC, tie tolerance 1e-6 favoring
  larger lambda, and at most two one-decade edge extensions as in the core
  implementation. Use the unchanged KKT backend and retry/failure policy.
  Capture prediction warnings as well as fit warnings. Primary probabilities
  use the same glum prediction pathway as inner tuning in the locked pipeline;
  compare their AUROC with the linear score as a numerical diagnostic only.
- Age: autosomal genes only, selected by positive training MAD with gene-ID tie
  breaking, at most 1,000; training mean/SD scaling. Use sklearn Ridge with
  an unpenalized intercept, solver `svd`, alpha 0.01, 0.1, 1, 10, 100.
  Its objective is sum of squared residuals plus alpha times squared slopes,
  not a sample-size-normalized loss. Select by mean inner-fold MAE; ties
  within 1e-6 favor larger alpha. Report the constant training-median-age
  comparator. These implementation choices are recorded before control scores.
- Use 5,000 percentile donor bootstraps, sex-stratified for sex AUROC and
  ordinary within-control donor resampling for age Spearman/MAE. Paired
  resamples also report the age MAE difference against the median comparator.
  Bootstrap scope strings are `controls/{train}/sex/bootstrap` and
  `controls/{train}/age/bootstrap`. Record undefined age correlations rather
  than substituting zero. Models are fixed during bootstrap.
- Development sex gate: AUROC >= 0.95 in **both** directions. BA at 0.5 is
  descriptive. Age capacity: Spearman >= 0.4 and MAE below the training-median
  comparator in both directions; failure only limits claims. Sex failure
  means the registered conditional confirmatory execution does not proceed,
  unless an independently demonstrated technical bug justifies an amendment
  and repetition of the same fixed control. No rescue transformation is chosen.

## Sigmoid warning diagnostic

Installed glum 3.4.1 `LogitLink.inverse` uses scipy expit and clips probabilities
to `[50*eps, 1-50*eps]` when extremes occur. This may create ties and hence may
affect inner AUROC. Passed fit KKT checks do not establish absence of this effect.

Keep the locked core pipeline and all original run files unchanged. Reconstruct
the existing manifest's training/validation partitions and replay the recorded
inner fits with the same preprocessing, penalties, weights and numerical backend.
For every replay, compare glum probabilities, expit probabilities and linear-score
AUROCs; record clipped-tail case/control pair counts, warning context, and
preprocessing hashes. Reconstruct the original candidate means and choices from
the recorded grid, and report whether unclipped or linear-score rankings would
change any choice. These alternatives are diagnostics, not replacement results.
Any material scoring correction requires a timestamped public amendment before
an affected scientific rerun. Report failed reproduction explicitly.

Synthetic checks precede real controls: preprocessing cannot learn from test
data; sex labels must not be confused with PD status; a known transferable sex
and age signal should be recoverable; bootstrap must preserve pairing; and the
sigmoid diagnostic must detect a deliberately constructed mixed-class tail tie.

Aggregate results go in Git; expression, identifiers, memberships, per-fit logs,
predictions and model arrays remain in ignored `private_data/`. Record hashes of
code, inputs, config and the locked core pipeline. Neither control results nor
diagnostics permit changing M1–M5 or selecting a new endpoint.
