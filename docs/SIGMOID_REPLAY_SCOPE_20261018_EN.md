# Retrospective sigmoid replay scope clarification

The initial all-scope diagnostic replay was stopped to bound this technical
investigation to the observed warning. The completed sex/age control outputs
are unaffected. `--reported-contexts-only` replays **all** candidates and five
folds for the six model/scope contexts in which the original execution reported
sigmoid warnings (GSE99039 outer 0/M5, outer 2/M3, transfer/M5 and pooled
final/M2/M3/M5), plus their M1/M4 baseline candidates: 490 inner fits. The recorded
fit time for these is about 28.5 minutes, versus roughly 100 minutes for all
scopes. This is a retrospective diagnostic scope, not a change in scientific
models, selection or endpoints. Results cannot establish reproduction of every
unreplayed scope. The first driver attempt failed before fitting because it
passed a duplicate context keyword; that driver bug was corrected. Initial and
interrupted diagnostic directories are retained, and no scientific result is
replaced by either attempt.

