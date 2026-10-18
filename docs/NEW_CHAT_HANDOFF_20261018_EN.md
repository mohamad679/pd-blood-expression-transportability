# New-chat continuation brief — 18 October 2026

This is a durable handoff from the long project conversation, updated after
the cloud run completed on 19 October 2026. Read it as context
for continuing the existing research, not as authorization to repeat completed
analyses. Later verified evidence supersedes this dated snapshot. Refer to the
published protocol/amendments for exact scientific rules rather than substituting
this summary for them.

## 1. Current state: cloud run complete, infrastructure removed

**The fixed B=200 exploratory GEO permutation control has completed.** Its
`SUMMARY.json` reports `status: complete`, `B: 200`, eight workers and
`PPMI_accessed: false`. All fixed replicate IDs 0–199 have completion markers
in the returned archive. The six aggregate results and their limits are in
[GEO_PERMUTATION_RESULTS_20261019_EN.md](GEO_PERMUTATION_RESULTS_20261019_EN.md).
The final summary is not the confirmatory external PPMI test.

The 159 MiB returned archive is available in ignored local storage at
`private_data/server_return/pdgeo-results-20261019.tar.gz` and also in
`~/Downloads/`. It passed archive inspection; SHA-256 matches the original
server file:

```text
50d9b536f5f0f46a176605859134e23f143bd3994840d885074dc79ca5a7fead
```

The user deleted the Hetzner VM after verification. Console screenshots showed
no servers and no Primary IPs. **There is no live cloud controller to inspect,
resume or restart.** Historical operational instructions in later sections
explain how this stage was done; do not treat them as current steps. The next
scientific work is bounded remaining diagnostics, PPMI identity/eligibility,
final GEO models/auxiliary control and the second public pre-PPMI lock. No PPMI
score has been computed.

## 2. Workspace, access and communication

- Local project folder: `/Users/mohsenshamsijazeb/Documents/ChatGPT/Potential paper`.
- Repository: `/Users/mohsenshamsijazeb/Documents/ChatGPT/Potential paper/pd-blood-expression-transportability`.
- Private remote: `https://github.com/mohamad679/pd-blood-expression-transportability.git`.
- Branch at handoff: `main`. Recent commits before this document: `12b6811`
  (read-only progress monitor), `54f62ed` (verified public amendment),
  `9ee6dc0` (encrypted transfer fallback), `727738e` (cloud preparation).
- New chats in the same local project can access attached files, subject to
  permissions; they do not automatically inherit this conversation transcript.
  The complete final server checkpoint archive is now held locally in ignored
  private storage.
- Speak Persian to the user. Keep operational answers short, prompt and concrete.
  Give one next step, specify **Mac Terminal** or **server after SSH**, and wait
  for the user's result before the next dependent step. Avoid long command menus.
- Repository documentation and reports are English. Commit/push meaningful
  completed changes, respecting immutable scientific bindings and data boundaries.
- The user now works on a Mac. Use Bash/zsh, not PowerShell. Do not send them
  back to Windows merely to review the completed GEO cloud job.
- Follow applicable `AGENTS.md`. Prefer codebase-memory-mcp for code discovery;
  if indexed functions are unavailable, use the documented fallback. Never
  assume a new chat has remote credentials or an authenticated SSH session.

## 3. Historical server operations — VM deleted after result return

The server details and commands below record the completed stage. The VM and
its Primary IP no longer exist; do not attempt SSH or recreate infrastructure
for an ordinary results review.

Hetzner VM `pdgeo-permutation`, IPv4 `188.245.16.195`, root directory
`/root/pdgeo`. CPX62 in Nuremberg: 16 shared AMD vCPUs (EPYC-Genoa), 32 GB RAM,
640 GB disk; Ubuntu 26.04.1 LTS x86_64. Exact Python 3.11.9 and the 20 locked
packages are installed, including glum 3.4.1. Do not change the environment.

**On the Mac**, the practical working SSH route has been:

```bash
ssh -o PubkeyAuthentication=no root@188.245.16.195
```

The user enters their root password privately. The assistant does not know it.
A dedicated encrypted key exists at `~/.ssh/pdgeo_hetzner`; its public key was
installed correctly, but the user had passphrase difficulties and password SSH
was the reliable route. Do not confuse the public key with its passphrase.
The firewall `pdgeo-ssh` permits TCP 22 from the user's public IPv4 /32 and
ICMP; outbound is allowed. If the user's IP changes, inspect that restriction
before changing anything broadly.

**On the server**, these are read-only inspection commands; give the one needed:

```bash
python3 /root/pdgeo-monitor.py --once
tail -n 30 /root/pdgeo/permutation.log
pgrep -af pipeline.permutation_executor
tmux list-sessions
```

The existing full-run session was launched as follows. This is a historical
record, **not an instruction to rerun it**:

```bash
tmux new-session -d -s pdgeo-run -c /root/pdgeo 'bash pipeline/server_stage.sh run 8 > permutation.log 2>&1'
```

Older sessions included `pdgeo`, `pdgeo-preflight`, and `pdgeo-pilot`.
An old pilot PID 3883 is historical; do not use it as the current process ID.
`duplicate session` is not permission to kill a session and relaunch a controller.
`Another controller is active` means inspect the existing controller first.

Read-only monitor source is `pipeline/monitor_permutations.py`; deployed copy is
`/root/pdgeo-monitor.py` (local temporary copy `/private/tmp/pdgeo-monitor.py`).
It does not fit models or write checkpoints. For the recurrent arrow-key echo
noise (`^[OA`, `^[OB`, `^[[B`), the display wrapper is:

```bash
bash -c 'trap "stty echo" EXIT; stty -echo; clear; python3 /root/pdgeo-monitor.py --interval 30'
```

Ctrl-C stops this **observer only**, not the detached computation. `stty sane`
can restore terminal settings if needed. Do not send Ctrl-C to an attached
foreground computation. The monitor counts completion markers but is not a
replacement for the executor's full artifact/hash validation.

During execution the user could turn off the Mac because tmux continued on
the VM. The original screenshot quoted about EUR 0.254/hour for the VM plus
IPv4; the monthly figure was a cap, not a required full-month purchase. The
VM and Primary IP were deleted after the result returned. The actual final
invoice has not been inspected.

The private output/checkpoint archive was verified before the user deleted the
VM. The subsequent Primary IP inventory was empty.

### Current output locations

- Run root: `/root/pdgeo/private_data/geo_permutations_cloud_v1/`.
- `RUN_BINDING.json`, `PLAN.json`, `IDENTITY_VERIFIED.json`, `RUNTIME_PILOT.json`.
- `rep_000` through `rep_199`, with atomic `COMPLETE.json` and attempt artifacts.
- Attempts retain `STARTED.json`, `progress.json`, `FAILURE.json`, `RESULTS.json`,
  warning records, fit logs, preprocessing/model artifacts and predictions.
- `SUMMARY.json` is produced only after all 200 valid completions. The monitor's
  COMPLETE state also requires the completed summary with B=200.
- `/root/pdgeo/permutation.log` and `/root/pdgeo/runtime_pilot.log`.
- `/root/pdgeo/results/cloud_preflight/preflight_passed.json` and
  `ZENODO_PUBLICATION_RECEIPT.json`.

The summary, logs and complete checkpoint tree were returned in the verified
private archive named in Section 1. Do not publish private GEO fit artifacts.
Inspect the actual output schema before postprocessing.

### Resolved infrastructure problems — do not repeat them

- Web/VNC console pasting corrupted shifted characters (`&&`, `%`, `>`), making
  SSH key commands ineffective. Use native Mac Terminal SSH for shell commands.
- The Mac `.zshrc` warning about a missing Google Cloud SDK path is unrelated.
- DNS initially failed; systemd-resolved was absent. `/etc/resolv.conf` was
  backed up as `/root/resolv.conf.before-pdgeo` and replaced with nameservers
  1.1.1.1 and 8.8.8.8. GitHub then resolved successfully. No need to repeat this
  unless a new DNS failure is demonstrated.
- Bulk SCP from Mac stalled at 2–6 KB/s. The GEO-only 37 MiB allowlisted bundle
  was instead encrypted and downloaded by the server over HTTPS in seconds.
  No PPMI data were included. It was checksum-verified and extracted already.
- Repo privacy caused unauthenticated GitHub release URLs to return 404. A
  short-lived signed asset URL in the transferred `fetch_geo_bundle.sh` solved
  it, without sending the GitHub account token to the server. That URL can
  expire; the original download command is not a durable public URL.
- The decryption pass is private, ignored local/server storage. Never print it,
  commit it, or ask the user to paste it into a chat. Setup is complete now.
- Two overlapping preflight controllers caused a lock error; the original
  controller completed. Its passing evidence was returned and inspected.

## 4. Verified public gates and provenance

The combined amendment is **published and effective**, version
`3.1.4-permutation-cbc`, DOI [10.5281/zenodo.23449296](https://doi.org/10.5281/zenodo.23449296).
Publication timestamp: `2026-10-18T13:46:35.803828+00:00`. All four public files
were downloaded unauthenticated and matched local bytes. Do not publish it
again or treat the historical draft title as the current execution gate.

```text
binding_sha256:
4be7bd1bc5e92da5d094fc6f59dbd45df033b254b3dc31043b4af8ed1ccc0310
preflight_report_sha256:
2b21d17727adcbe71022a411e4db870b46d616d34d8839da562b0b5b7f80f0ac
cloud config sha256:
740a1355c8c94df47a508ca336a09479a5d324524e374fdb1a8043425cfe985e
```

Linux preflight passed 33 checks and 650 synthetic fit attempts (max KKT
2.6839175878284882e-8). Eight workers/two threads were observed. Conservative
memory budget: 5,996,040,192 bytes against 32,848,162,816 bytes physical RAM;
zero `pswpout` increase. All 20 observed M4/M5 fold comparisons reproduced
candidate/selection decisions, AUROC error zero, max probability error
7.649436639667329e-14. The 10.302-second serial / 2.681-second parallel
synthetic benchmark is not a full-run ETA. The old Mac swap gate failure is
historical and was not waived; Linux passed its own frozen gate.

Evidence: `results/cloud_preflight/preflight_passed.json`,
`results/cloud_preflight/ZENODO_PUBLICATION_RECEIPT.json`, and
`results/zenodo_permutation_cbc_v3_1_4/` (public record, verification and files).
The returned archive is also in the user's Downloads as `PREFLIGHT_RETURN.tar.gz`.

Zenodo lineage, author Mohammad Javad Asgari, ORCID `0009-0002-2313-571X`
(do not invent an affiliation):

| Record | DOI suffix |
|---|---|
| Concept | 23405370 |
| Parent v3.1 | 23405371 |
| Implementation v3.1.1 | 23405802 |
| English translation v3.1.1-en | 23406064 |
| English provenance-corrected v3.1.2-en | 23411084 |
| Completed GEO milestone v3.1.3-geo-cv | 23445851 |
| Joint permutation/CBC v3.1.4 | 23449296 |

All suffixes are under `10.5281/zenodo.`. English protocol translations were
created after the first PROBE fit; original Persian registration and published
amendments govern discrepancies. Never claim every English file predates all
fits or that translated bytes have source-language hashes. The second public
pre-PPMI scoring lock is **still outstanding**.

Frozen configuration/certificate files retain historical status booleans and
older repetition counts where appropriate to their original snapshot. Do not
edit bound bytes to make their status text look current. Use the latest public
amendment, receipts and living status instead.

## 5. Scientific objective and locked design

This is a methodology/transportability study of whole-blood PD transcriptome
classifiers and age/sex/cell/platform shortcuts, including added gene value over
a strong baseline. It is beyond feasibility but is not a finished external
validation or paper. Q1/Q2 quality is an aim, not a promise; positive effects,
clinical utility and causal shortcut fractions are not assumed.

Development cohorts: GSE57475/PROBE, n=142 (93 PD, 49 controls), and the
GSE99039 GENEPARK-TR subset, n=239 (120 idiopathic PD, 119 controls; 199 other
records excluded). Total 381, 16,364 shared genes (15,720 autosomal, 629 X,
15 Y). GENEPARK is one TRO site with metadata/batch limitations. PROBE uses
GSM as operational donor identity with incomplete provenance; do not claim
biological independence has been proved.

PPMI Project 133, IR3 B38 baseline Long RNA-seq, is external only. The roughly
42 GB original archive was on Windows; do not assume it is on this Mac. Full
sample-ID segment parsing corrected an earlier first-dot parser: 4,871 count
files linked one-to-one to release identities; 4,752 filename/header agreed,
119 were header-only. Provisional 311 PD / 161 healthy controls is **not a final
locked population**. Identity/alias, eligibility, disease concordance, gene
coverage and private crosswalk lock remain. Feasibility metadata/labels were
seen; do not describe a completely untouched metadata holdout.

At handoff **no PPMI disease, sex, age or cell-marker model scores were computed**.
PPMI person-level data stay local, outside Git and the cloud VM. Sex concordance
cannot detect same-sex swaps; longitudinal cosine is diagnostic, not automatic
reassignment/exclusion. Identity failures cannot be rescued by outcome tuning.

Preprocessing: valid array probes map to unique Entrez/GENCODE v29 protein-coding
canonical genes and aggregate by median. Within-sample midrank percentile is
`(rank - 0.5)/G`; RNA counts use exon-union length adjustment before ranking.
Select at most 1,000 positive-MAD genes using training data only; pooled
selection uses the minimum within-training-source MAD. Fit training-only
means/SD; no PPMI centering/tuning. Ties/zeros/platform noise remain limitations.

All M1–M5 use the same weighted binomial block-penalty logistic family:

- M1: age, age squared, sex, age-by-sex; unpenalized baseline.
- M2: genes only. M3: M1 plus genes.
- M4: M1 plus six fixed marker-panel means. M5: exact M4 baseline plus genes.
- Baseline penalty factors are zero and gene factors one; no artificial
  large-column scaling to approximate an unpenalized block.
- Panels: T 14, B 8, monocytic 7, neutrophils 15, erythroid 11, platelet 118
  genes. Means of percentile expression are proxies, not measured fractions,
  full MCP-counter or xCell. Erythroid/globin may reflect platform processing.

Weights `n/(2*K*n[source,class])` give equal source-by-class mass and sum to n,
across models and inner/outer/transfer/final scopes including K=1. Preprocessing
and external evaluation are unweighted. Kish n: pooled 334.046 versus 381;
PROBE 128.366; GENEPARK 238.996. Within-validation-source macro-AUROC avoids
pooled source/prevalence shortcuts (source-only pooled AUROC was 0.5724765).

Nested five outer by five inner donor-separated, class-stratified folds; two
development sources, not a three-source LOCO design. Seed 20261014. Lambda grid
`[1e-4, 1e-3, 1e-2, .1, 1]`, rho `[0, .5, 1]`; ties within 1e-6 favor
baseline, larger lambda, then larger rho. Only two bounded decade expansions
are allowed for selected nonzero-gene edge candidates. No post-result grid rescue.

Numerical validity: weighted normalized mean-NLL KKT <=1e-6. Original solver
max iterations 1,000 / inner 100,000, gradient tolerance 1e-7, backend predictor
scaling disabled; one cold retry 5,000 / 1,000,000 / 1e-8 for the same candidate.
An invalid candidate invalidates its dependent analysis; never silently drop
candidates/folds/replicates. A warning alone is not a numerical failure.

Primary endpoint remains external PPMI AUROC(M5)-AUROC(M4), one confirmatory
comparison, 5,000 paired class-stratified donor percentile-bootstrap samples,
fixed trained models. Delta*=0.05 is a research margin, not clinically validated.
Exclude gains >=0.05 only if the upper 95% CI is below 0.05; near-zero noisy
results may be inconclusive. If M5 selects no genes and equals M4, report that
branch; a degenerate difference CI is not general proof of no gene information.
AUROC is primary; prevalence-aware AUPRC, Brier, log loss and calibration are
descriptive. Case-control probabilities are not population risk.

## 6. Completed scientific work and findings

All 13 core GEO scopes completed: ten outer-fold analyses, two directional
transfers and one pooled final fit. Do not rerun them simply because a new chat
starts. Mean outer AUROC M4/M5: PROBE 0.56146/0.59111; GENEPARK
0.66295/0.68608. Transfer M5-minus-M4: -0.178 PROBE-to-TRO, -0.009 reverse.
These are development results, not external PPMI performance.

Pooled final M5 lambda .001, rho 1: **L1/lasso**, 210 nonzero genes of 1,000;
M3 lambda .01/rho 1 also L1; M2 lambda .1/rho 0 is L2. An already released short
report wrongly called final M5 ridge. `REPORTING_CORRECTION_20261018_EN.md`
corrects interpretation without rewriting frozen report/model bytes.

All 3,146 post-pilot fit records passed KKT (max about 4.1093e-8). All 3,415
fit-record membership, train-only preprocessing and Kish checks across the
13 original scopes matched. Six sigmoid warning contexts were reproduced in
a 490-fit audit; explicit clipping changed neither AUROC nor candidate selection.
No scoring-method amendment/rerun was justified by that audit.

GEO controls are complete. Sex uses only five available Y markers (KDM5D,
DDX3Y, RPS4Y1, UTY, EIF1AY); XIST absent, do not add it. Control-only ridge
transfer AUROC 1.000 and 0.993 passed the two development directions. External
PPMI sex gate is AUROC >=.95 and CI lower >=.90 in controls, still unexecuted.
Threshold-.5 balanced accuracy is descriptive and can be .5 despite high AUROC.
Final GEO-trained auxiliary sex fit still needs inclusion in the second lock.

Autosomal control-only ridge age model uses alpha .01/.1/1/10/100, MAE tuning,
larger alpha for ties. Spearman -0.058 and .259 and no MAE gain over the
training median did not establish age transfer capacity. Both selected alpha
100 at the top boundary. Do not expand the grid after seeing this; distinguish
poor transport from possible fixed-grid limitations. Age is non-stopping.

## 7. Active permutation control — exact scope

Published amendment fixes **200 IDs 0–199, M4/M5 only, ten within-GEO outer
folds each**. It excludes directional transfers and pooled final fitting.
Pilot 0 counts within 200, not in addition. The scientific config on the server
is `protocol/permutation_config_cloud_v1.json`, not the old Mac two-worker file.

Shuffle one coherent donor-label map within cohort by fixed outer-test-fold
block, permuting both training and evaluation labels, preserving block class
counts. Rebuild inner stratification from current permuted training labels.
Recompute all label-dependent weights, tuning, allowed expansions and fits;
label-free caches are valid only with matching training-membership/code hashes.
Seed uses the first four big-endian SHA256 bytes of
`20261014:permutation/{rep_id}/{cohort}`, with donors sorted before shuffling.
Worker scheduling/restart must not change labels or statistics.

Per cohort, unweighted mean of five outer AUROCs for M4, M5 and their paired
delta. Six exploratory p-values use `(1 + count[Tb >= Tobs - 1e-12])/201`.
Minimum p is 1/201 (~.00498); Monte Carlo SE near .05 is ~.015. This is a global
association null; delta under global shuffling is not a conditional test of
added genes given baseline. Chance-like results do not certify absence of leakage.

No outcome-driven stopping, seed replacement, shortened B/denominator, skipped
terminal failures, or lock deletion. Genuine interruption can resume validated
completed replicates with the same seeds and retained attempts. Inspect the
executor and original error before any recovery; don't treat a terminal
numerical failure as an interruption. Report invalid/incomplete analyses honestly.

## 8. CBC correction and published selection rule

User supplied current and archived Blood Chemistry & Hematology and archived
Whole Blood Sample Collection exports. Among the nine fixed numeric markers,
56,104 shared exact measurement keys had identical values and units, zero
differences/ambiguous duplicates. This is not whole-table identity.

Earlier exact-visit matching incorrectly added a restriction and yielded 0/161
baseline controls. That broad non-evaluability conclusion was withdrawn. All
161 have SC CBC rows; 122 have known BL collection month; **31** have same-month
CBC rows (all at SC) in both exports. 31 is an availability upper bound, not
marker-specific numeric n. No association/expression score was calculated.

The published joint amendment uses current export as primary, archive for
concordance only. Same participant, any visit, CBC collection month/year equal
to BL whole-blood collection month/year; exclude adjacent months. Prefer BL,
then SC, then other visits per marker. Collapse only identified exact duplicates
with agreeing values/units. At the highest available priority tier, unresolved
multiple collection occasions/conflicts exclude the participant-marker; never
fall back to a lower tier or break ties by CSV order or result magnitude.
Monthly dates cannot establish same/nearest/earliest day and omit some genuine
within-30-day pairs across adjacent months. `BLDCOLL` codes remain uninterpreted
and are neither an eligibility filter nor an execution blocker.

Per-marker usable n>=30, Spearman>=.4 and CI lower>0 provide limited supportive
evidence. n<30 is non-evaluable; do not change matching rules to inflate n.
CBC does not certify all autosomal/PD transfer, total lymphocytes are not T/B
specific, RBC is not reticulocytes, and globin processing remains a limitation.
CBC and age are non-stopping diagnostics; neither changes fixed M4/M5 panels.

## 9. Remaining work after the completed run

1. Review the returned full B=200 summary, numerical validity, warnings and
   checkpoint evidence using the private local archive. Aggregate results are
   documented; deeper analysis can continue without a VM. The complete
   archive has already been transferred and verified, and the user deleted
   the VM and Primary IP.
2. Complete remaining bounded diagnostics/sensitivities from the governing
   protocol. The original 13-scope membership/preprocessing audit is already
   complete: do not repeat it without a new reason. Remaining work includes
   bounded leakage/provenance checks, bootstrap coverage simulation, registered
   X/Y removal, train-only age/sex residualization, age/sex matching within
   three years without replacement, exclusion of panel genes from the M5-added
   gene block, technical metadata diagnostics and the development ComBat audit
   as specified. Inspect exact protocol before implementing; do not invent
   replacement controls or infer completion from this list.
3. Finalize PPMI identity/crosswalk, eligibility, gene coverage, QC exclusions
   by diagnosis, treatment/source metadata and publication conditions. Private
   person-level material remains local. CBC numeric completeness/selection
   metadata can be audited under the published rule while waiting; any PPMI
   expression-derived scoring requires its proper prior lock.
4. Fit/freeze the required final GEO auxiliary sex model; verify saved final
   M1–M5 models/preprocessing/splits. Publish the **second pre-PPMI lock** of
   code/environment, trained models/preprocessing and a non-identifying hash/
   reference for the private population/crosswalk **before any PPMI score**.
5. Perform prespecified external controls. Identity/sex/transformation failures
   stop confirmatory PD scoring unless a proven technical implementation bug
   is transparently corrected/amended before repetition. Otherwise report
   non-evaluability. Age, CBC and longitudinal diagnostics do not authorize
   target tuning. Then execute fixed PPMI disease scoring only when gates pass.
6. Compute registered external metrics/paired intervals and sensitivities with
   no target refitting. Interpret null/inconclusive findings, limitations and
   feasibility honestly; write the manuscript and final release. Finishing
   permutations alone does **not** finish this project.

## 10. Reading map and evidence precedence

Read `docs/STATUS.md` and this handoff first. For detailed science/history:

- `docs/COMPLETE_PROJECT_REPORT_20261018_EN.md` — retrospective full report.
- `docs/PROJECT_CONTINUATION_REPORT_20261018_EN.md` and
  `docs/PROJECT_HISTORY_AND_ROADMAP.md` — continuation/history.
- `docs/GEO_CV_RESULTS_20261018_EN.md`,
  `results/geo_development_20261018/aggregate_summary.json`.
- `docs/GEO_TRANSFER_CONTROLS_20261018_EN.md`,
  `results/geo_controls_20261018/aggregate_summary.json`.
- `docs/SIGMOID_AUDIT_20261018_EN.md`,
  `docs/REPORTING_CORRECTION_20261018_EN.md`.
- `results/permutation_preflight_20261018/geo_membership_audit.json`.
- `docs/PPMI_CBC_ANY_VISIT_AUDIT_20261018_EN.md` and
  `docs/PPMI_COLLECTION_LINKAGE_AUDIT_20261018_EN.md`.
- `protocol/PROTOCOL_v3_1_EN.md`, `protocol/PREREG_v3_1_FULL_EN.md`,
  `protocol/PROVENANCE_CORRECTION.md`, `protocol/TRANSLATION_SOURCE_MAP.json`,
  `protocol/TRANSLATION_DIFFERENCES.md`, and published amendment files in
  `results/zenodo_permutation_cbc_v3_1_4/`.
- `docs/HETZNER_GEO_RUNBOOK_EN.md` is the original operational guide; its
  initial not-yet-executed statements are historical. Its setup steps are
  complete and must not be repeated as a default next action.
- `protocol/PERMUTATION_AMENDMENT_DRAFT_20261018_EN.md`,
  `docs/PERMUTATION_EXECUTOR_PREFLIGHT_20261018_EN.md`, and
  `docs/CLOUD_HANDOFF_20261018_EN.md` are historical preparation records,
  superseded for current execution status by the published amendment and
  returned Linux certificate. Preserve their original bytes/provenance.

Frozen core pipeline SHA256:
`6f7f6cf42608002b185a63044387c9b9f5e624fc0b9ff5ae156030e3655e1a2c`.
Split manifest SHA256:
`2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2`.
The certificate lists every bound source/config hash. Check those before any
scientific change; read-only monitoring/documentation changes are not permission
to alter the engine. Four ignored `prepared/` GEO inputs are frozen;
`data/manifests/input_audit.json` records their hashes. A verified 167-file local
backup of original post-pilot GEO outputs is in ignored
`private_data/geo_runs_20261018/`. Server permutation outputs are a different
directory and still require return/backup.

## First reply expected from the next assistant

After reading this handoff and the aggregate result report, briefly confirm in
Persian that the B=200 run completed, the private archive is verified locally,
and the Hetzner VM/Primary IP were removed. PPMI scoring and its second public
lock are still pending. Continue with one concrete next research step based on
the published protocol; do not attempt server SSH, request secrets or recite
the whole history.
