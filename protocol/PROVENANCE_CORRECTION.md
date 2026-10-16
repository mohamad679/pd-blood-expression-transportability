# Provenance correction to the first English edition

This public version corrects the English presentation published as
https://doi.org/10.5281/zenodo.23406064. It does not alter the original
preregistration (https://doi.org/10.5281/zenodo.23405371), the pre-execution
implementation amendment (https://doi.org/10.5281/zenodo.23405802), any model,
threshold, seed, grid, sample, code, or outcome. One real PROBE fold started
2026-10-15 22:07:32 UTC, before either English edition. Its results were viewed.
The correction is public only from this new Zenodo version's actual timestamp.

The first English edition incorrectly named newly created English files as the
documents saved or hashed before synthetic output. The actual source files were:

| Original source | SHA-256 | Evidence and English counterpart |
|---|---|---|
| `SYNTHETIC_V3_1_PLAN_FA.md` | `45d15dc4576b92301016250bdf5eb369eb2a3b0a56499e4ef1eb24e6544fab1a` | Original plan, translated as `SYNTHETIC_V3_1_PLAN_EN.md` |
| `KKT_STEP_ZERO_v3_1_FA.md` | `4c59a9d6c95b223edfcffcaeaad91cad0e4a7139f42068d180e6628064db5cfb` | Named in the original `KKT_STEP_ZERO_PREEXEC_SHA256.json`; translated later as `KKT_STEP_ZERO_v3_1_EN.md` |
| `KKT_NUMERICAL_AMENDMENT_20261016_FA.md` | `ffad1329e0493ddd2f83536a77f3726901eabe910909ad23020e71dcae0ebe5f` | Original solver amendment, translated later as `KKT_NUMERICAL_AMENDMENT_20261016_EN.md` |

The original Persian v3.1 text and public v3.1.1 amendment govern if a
translation differs. Status statements inside translations refer to their
original 16 October 2026 freeze, not the current publication state. The full
English pilot report is included in this ZIP. The methods and the actual
time ordering are unchanged; the presentation and source attribution are
corrected. The [difference table](TRANSLATION_DIFFERENCES.md) records these
changes. No PPMI disease, sex, age, or cell-marker score has been calculated.

The 1,000-development-permutation scope in the original protocol is **not**
reduced by this correction. Any change to its model set or count needs a
separate, timestamped methodological amendment before the first affected
permutation run. It remains pending.
