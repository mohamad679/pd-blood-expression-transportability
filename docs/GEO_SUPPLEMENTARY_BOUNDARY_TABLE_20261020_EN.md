# Supplementary Table S1. Persistent lambda-boundary selections in GEO

A boundary is listed only when the saved decision's `persistent_boundary` flag is true after two permitted grid expansions. These are descriptive tuning decisions, not performance estimates. All integrated final runs are included; timing pilots and interrupted copies are excluded. GSE57475 is PROBE, and GSE99039 is GENEPARK. A transfer scope names its training cohort.

| Analysis | Scopes | M2 flagged / examined | M3 flagged / examined | M5 flagged / examined | All flagged / examined |
|---|---:|---:|---:|---:|---:|
| Original | 13 | 1/13 | 0/13 | 1/13 | 2/39 |
| X/Y exclusion | 13 | 3/13 | 0/13 | 0/13 | 3/39 |
| Age/sex residualization | 13 | 1/13 | 1/13 | 0/13 | 2/39 |
| Panel-gene exclusion | 13 | 0/0 | 0/0 | 2/13 | 2/13 |
| Age/sex matching | 6 | 1/6 | 0/6 | 1/6 | 2/18 |
| ComBat | 39 | 7/39 | 5/39 | 3/39 | 15/117 |

| Analysis | Scope or ComBat arm / scope | Model | Lambda | Rho | Grid expansions |
|---|---|---:|---:|---:|---:|
| Original | GSE57475/outer/2 | M2 | 1e-06 | 0.5 | 2 |
| Original | GSE57475/transfer | M5 | 1e-06 | 1.0 | 2 |
| X/Y exclusion | GSE57475/outer/0 | M2 | 1e-06 | 1.0 | 2 |
| X/Y exclusion | GSE57475/outer/2 | M2 | 1e-06 | 0.5 | 2 |
| X/Y exclusion | GSE57475/outer/4 | M2 | 1e-06 | 0.0 | 2 |
| Age/sex residualization | GSE57475/transfer | M2 | 1e-06 | 0.5 | 2 |
| Age/sex residualization | GSE57475/transfer | M3 | 1e-06 | 0.5 | 2 |
| Panel-gene exclusion | GSE57475/outer/0 | M5 | 1e-06 | 0.0 | 2 |
| Panel-gene exclusion | GSE57475/transfer | M5 | 1e-06 | 0.5 | 2 |
| Age/sex matching | GSE99039/outer/1 | M2 | 1e-06 | 0.0 | 2 |
| Age/sex matching | GSE99039/outer/1 | M5 | 1e-06 | 0.0 | 2 |
| ComBat | inductive_no_disease/GSE57475/outer/0 | M2 | 1e-06 | 0.0 | 2 |
| ComBat | inductive_no_disease/GSE57475/outer/2 | M2 | 1e-06 | 0.5 | 2 |
| ComBat | inductive_no_disease/GSE57475/transfer | M2 | 1e-06 | 0.5 | 2 |
| ComBat | inductive_preserve_disease/GSE57475/outer/0 | M2 | 1e-06 | 0.0 | 2 |
| ComBat | inductive_preserve_disease/GSE57475/outer/1 | M2 | 1e-06 | 0.0 | 2 |
| ComBat | inductive_preserve_disease/GSE57475/transfer | M2 | 1e-06 | 0.5 | 2 |
| ComBat | inductive_preserve_disease/GSE57475/transfer | M3 | 1e-06 | 0.5 | 2 |
| ComBat | inductive_preserve_disease/GSE99039/outer/2 | M5 | 1e-06 | 0.0 | 2 |
| ComBat | transductive_no_disease/GSE57475/outer/1 | M3 | 1e-06 | 0.5 | 2 |
| ComBat | transductive_no_disease/GSE57475/outer/3 | M3 | 1e-06 | 1.0 | 2 |
| ComBat | transductive_no_disease/GSE57475/outer/3 | M5 | 1e-06 | 0.0 | 2 |
| ComBat | transductive_no_disease/GSE57475/outer/4 | M3 | 1e-06 | 0.5 | 2 |
| ComBat | transductive_no_disease/GSE57475/transfer | M3 | 1e-06 | 0.0 | 2 |
| ComBat | transductive_no_disease/GSE57475/transfer | M5 | 1e-06 | 0.5 | 2 |
| ComBat | transductive_no_disease/GSE99039/outer/1 | M2 | 1e+02 | 0.0 | 2 |

The separate GEO sex-marker capacity control selected lambda=100 at the upper edge in both directional training sources after two extensions. It used a distinct auxiliary model and grid and is excluded from the M2/M3/M5 counts above.

Source: saved final-run decision artifacts; extraction method in `docs/GEO_BOUNDARY_INVENTORY_METHOD_20261020_EN.md`. No PPMI data were accessed.
