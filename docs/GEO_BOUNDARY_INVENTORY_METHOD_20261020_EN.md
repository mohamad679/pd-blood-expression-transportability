# Method fixed before the GEO lambda-boundary inventory

This inventory reads the completed, integrated GEO run directories only: the original run, X/Y exclusion, age/sex residualization, panel-gene exclusion, age/sex matching, and three ComBat arms. Timing pilots and interrupted copies are excluded. It performs no fit, score, or model selection.

For each saved M2, M3, or M5 decision, count a persistent boundary only when its recorded `persistent_boundary` field is `true`. Report the selected lambda, rho, and number of grid expansions. Count every analyzed model decision in the denominator, including pooled/final training decisions, but distinguish those from held-out evaluation scopes. Include the panel-gene sensitivity's M5-only decisions. Summarize each analysis and list every flagged decision in Supplementary Table S1. An auxiliary sex-control edge selection is described separately because it uses a different tuning grid and is not an M1–M5 disease-model decision.

The table is a descriptive inventory of tuning behavior, not an additional performance analysis or an independent test. Boundary selections do not by themselves establish numerical failure or overfitting. No PPMI data are accessed.
