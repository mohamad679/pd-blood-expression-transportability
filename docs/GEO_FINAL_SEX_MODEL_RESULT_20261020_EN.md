# Final GEO-control sex model: fit and audit

The prespecified five-Y-marker auxiliary model was fitted on 168 frozen GEO healthy controls (49 PROBE/GSE57475 and 119 GENEPARK/GSE99039; recorded sex counts 108 and 60). PPMI data were not loaded. This final pooled model is reserved for the PPMI *control* transfer gate after the public second lock and does not change M1–M5 or the completed directional GEO control results.

Five source-and-sex-stratified inner folds were frozen before tuning. Among 36 successful fits, lambdas 1, 10, and 100 tied at mean inner AUROC 0.9984848485. The registered tie rule chose the larger lambda. After the maximum two one-decade extensions, the selected lambda of 100 remained on the upper grid boundary; no further search was performed. The number is an inner GEO control tuning result, not external discrimination.

Independent post-run checks confirmed the artifact SHA-256 against the actual `sex_model.npz`, 36/36 passed KKT records, maximum candidate KKT 2.2643e-08, final KKT 2.9965e-11, zero fit or prediction warnings, five disjoint training/validation folds, ordered five-marker arrays, and the registered tie decision. The private model, fold membership, and fit logs remain under ignored `private_data/` on the researcher's Mac. The committed aggregate manifest contains no participant identifiers, coefficients, or predictions.

| Frozen item | SHA-256 |
|---|---|
| Final model artifact | `aa92a1516328477b1c4ba1a61e7bcec43263cc38574f57ec99c5f26465348200` |
| Inner split manifest | `6d3bde988887349c089d60e66899833deb3364b8a0a5848fe8f65312f9bd1e19` |
| Fit script | `365f8b8dbf1f523611831687c6556611291241c6780fc93a1f64e76962229f0d` |
| Pre-fit decision | `eb672d073aad4815dfbe98170afbd69adca0343d6a70fbbf3af88ee83ed01ba3` |

Sources: [pre-fit decision](GEO_FINAL_SEX_MODEL_DECISION_20261020_EN.md), [aggregate result](../results/geo_final_sex_model_20261020.json), and private run manifest `private_data/geo_final_sex_model_20261020/SUMMARY.json`. The second public lock and all PPMI expression-derived evaluations remain pending.
