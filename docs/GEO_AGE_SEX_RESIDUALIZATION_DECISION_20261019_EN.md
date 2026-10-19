# GEO age/sex residualization implementation decision — 19 October 2026

The v3.1 configuration registers `train_OLS_only` but does not specify the
covariate basis or whether gene selection precedes residualization. Before any
full sensitivity run, use the following fixed implementation:

1. In each inner-training partition, fit an **unweighted** per-gene OLS using
   an intercept and the same covariate basis as M1: training-standardized age,
   age squared, recorded sex, and age-by-sex. Reuse that partition's age mean
   and SD and OLS coefficients unchanged to residualize its validation data.
   Standardize age using that training partition's own mean and SD **before**
   constructing age squared and age-by-sex. Repeat from scratch in the outer
   training set for the final refit; never estimate nuisance coefficients
   from the outer holdout.
2. For pooled training (including the pooled final fit), fit one OLS across
   the entire training partition, without a source indicator or source-specific
   coefficients. The same pooled training coefficients and age standardization
   are applied to a held-out cohort; separate source coefficients would not be
   defined for an unseen source such as PPMI.
3. Residualize the added gene block **before** applying the existing
   source-stratified positive-MAD selection, gene-ID tie break, top-1,000 cap,
   and training-only mean/SD scaling. This makes selection target variation
   remaining after the registered age/sex adjustment.
4. Leave the existing M1/M4 baseline construction and marker panel inputs
   unchanged; the sensitivity concerns only the added gene block.

Interpretation is bounded in advance: because age differs between PD cases and
controls, residualization may remove disease-related variation correlated with
age as well as age-associated shortcut variation. Therefore a performance drop
after residualization, by itself, is not evidence that the removed signal was a
shortcut.

The one-scope pilot checks that the fit code is invoked separately for every
inner-training partition, OLS residuals satisfy training normal equations,
baseline transforms remain exactly unchanged, and no fit hash uses an outer
holdout sample. This is an implementation choice resolving details left open
by the preregistration, not a new sensitivity hypothesis.
