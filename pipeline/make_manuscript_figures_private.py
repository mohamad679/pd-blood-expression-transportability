#!/usr/bin/env python3
"""Make aggregate manuscript figures; keep PPMI-derived panels owner-only.

Reads saved predictions locally without printing individual rows or identifiers.
No fitting, threshold selection, or recalibration is performed.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="pd-figures-mpl-"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve

from paired_bootstrap import paired_bootstrap


ROOT = Path(__file__).resolve().parents[1]
COLORS = {"M4": "#007E87", "M5": "#C84B3A"}
TRANSFER = {"PROBE → GENEPARK": "GSE57475", "GENEPARK → PROBE": "GSE99039"}


def saved(path: Path, m5_key: str = "M5") -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with np.load(path, allow_pickle=False) as z:
        y = np.asarray(z["y"], dtype=int)
        p4 = np.asarray(z["M4"], dtype=float)
        p5 = np.asarray(z[m5_key], dtype=float)
        samples = np.asarray(z["samples"]) if "samples" in z.files else np.arange(len(y))
    if (set(np.unique(y)) != {0, 1} or len(samples) != len(set(map(str, samples)))
            or not all(np.isfinite(p).all() and np.all((p >= 0) & (p <= 1)) for p in (p4, p5))):
        raise ValueError("Invalid saved prediction block")
    return y, p4, p5, samples


def save_figure(fig, out: Path, stem: str) -> None:
    fig.savefig(out / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(out / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def roc_figure(ppmi_file: Path, aggregate: dict, out: Path) -> None:
    panels = []
    for label, scope in TRANSFER.items():
        panels.append((label, *saved(ROOT / f"private_data/geo_runs_20261018/{scope}_transfer/heldout_predictions.npz")[:3]))
    y, p4, p5, _ = saved(ppmi_file)
    if len(y) != 472 or any(abs(roc_auc_score(y, pred) - aggregate["disease"]["metrics"][name]["AUROC"]) > 1e-12
                            for name, pred in (("M4", p4), ("M5", p5))):
        raise ValueError("PPMI saved predictions do not reproduce primary aggregate")
    panels.append(("GEO pooled → PPMI", y, p4, p5))
    fig, axes = plt.subplots(1, 3, figsize=(11.8, 3.6), sharex=True, sharey=True)
    for ax, (title, truth, a, b) in zip(axes, panels, strict=True):
        for name, pred in (("M4", a), ("M5", b)):
            fpr, tpr, _ = roc_curve(truth, pred)
            ax.plot(fpr, tpr, lw=2.1, color=COLORS[name], label=f"{name}  {roc_auc_score(truth, pred):.3f}")
        ax.plot([0, 1], [0, 1], color="#888888", lw=.8, ls="--")
        ax.set(title=title, xlabel="False-positive rate", xlim=(0, 1), ylim=(0, 1))
        ax.legend(frameon=False, loc="lower right", fontsize=8)
        ax.grid(alpha=.15)
    axes[0].set_ylabel("True-positive rate")
    fig.suptitle("Fixed-model discrimination across cohort transfers", fontsize=12)
    fig.tight_layout()
    save_figure(fig, out, "figure_1_transfer_ROC_private")


def forest_figure(out: Path) -> None:
    arms = [
        ("Original", Path("private_data/geo_runs_20261018"), "original"),
        ("X/Y excluded", Path("private_data/geo_xy_sensitivity_full_20261019"), "standard"),
        ("Age/sex residualized", Path("private_data/geo_age_sex_resid_full_20261019"), "standard"),
        ("Panel genes excluded", Path("private_data/geo_panel_gene_exclusion_full_20261019"), "standard"),
        ("ComBat: no disease", Path("private_data/geo_combat_sensitivity_20261020_retry1/inductive_no_disease"), "standard"),
        ("ComBat: disease in train", Path("private_data/geo_combat_sensitivity_20261020_retry1/inductive_preserve_disease"), "standard"),
        ("ComBat: transductive", Path("private_data/geo_combat_sensitivity_20261020_retry1/transductive_no_disease"), "standard"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), sharey=True)
    for ax, (direction, scope) in zip(axes, TRANSFER.items(), strict=True):
        for k, (label, base, style) in enumerate(arms):
            folder = f"{scope}_transfer" if style == "original" else f"{scope}__transfer"
            key = "M5_panel_gene_excluded" if label == "Panel genes excluded" else "M5"
            y, p4, p5, ids = saved(ROOT / base / folder / "heldout_predictions.npz", key)
            estimate, _ = paired_bootstrap(y, p4, p5, ids, repetitions=5000, seed=20261014)
            delta, lo, hi = estimate["delta_AUROC"], *estimate["CI95"]
            ax.errorbar(delta, k, xerr=np.array([[delta - lo], [hi - delta]]), fmt="o",
                        color="#334F77" if k < 4 else "#A15A29", capsize=2.5, ms=4)
        ax.axvline(0, color="#777777", lw=.9)
        ax.set(title=direction, xlabel="M5 − M4 AUROC (descriptive 95% CI)", xlim=(-.32, .17))
        ax.set_yticks(range(len(arms)), [a[0] for a in arms])
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=.15)
    fig.suptitle("GEO transfer sensitivity: paired fixed-score bootstrap", fontsize=12)
    fig.text(.5, -.015, "ComBat transductive arm uses held-out expression; all CIs condition on fixed fits.",
             ha="center", fontsize=8)
    fig.tight_layout()
    save_figure(fig, out, "figure_2_GEO_transfer_sensitivity_forest")


def reliability_figure(secondary: dict, out: Path) -> None:
    fig, (ax, counts) = plt.subplots(2, 1, figsize=(5.9, 6.7), sharex=True,
                                    gridspec_kw={"height_ratios": [3, 1]})
    ax.plot([0, 1], [0, 1], color="#888888", lw=.9, ls="--", label="Perfect calibration")
    centers = np.arange(10) / 10 + .05
    for name in ("M4", "M5"):
        bins = secondary["models"][name]["calibration_bins"]
        xx = [b["mean_predicted_PD"] for b in bins if b["n"]]
        yy = [b["observed_PD_fraction"] for b in bins if b["n"]]
        ax.plot(xx, yy, marker="o", ms=4, lw=1.7, label=name, color=COLORS[name])
        counts.step(centers, [b["n"] for b in bins], where="mid", lw=1.7, color=COLORS[name], label=name)
    ax.set(ylabel="Observed PD fraction", xlim=(0, 1), ylim=(0, 1),
           title="PPMI calibration of frozen GEO models")
    ax.grid(alpha=.15)
    ax.legend(frameon=False, fontsize=8)
    counts.set(xlabel="Mean predicted PD probability / bin midpoint", ylabel="n per bin")
    counts.grid(axis="y", alpha=.15)
    counts.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    save_figure(fig, out, "figure_3_PPMI_reliability_private")


def batch_table(out: Path) -> None:
    with (ROOT / "results/geo_technical_metadata_audit_20261019/class_by_batch.csv").open(newline="") as handle:
        data = [r for r in csv.DictReader(handle) if r["source"] == "GSE99039"]
    if len(data) != 11 or sum(int(r["PD"]) + int(r["control"]) for r in data) != 239:
        raise ValueError("GENEPARK technical batch table differs from audit")
    data.sort(key=lambda r: r["batch_or_chip_proxy"])
    rows = [[r["batch_or_chip_proxy"], r["control"], r["PD"], r["n"], f'{float(r["PD_percent"]):.1f}%'] for r in data]
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    ax.axis("off")
    tbl = ax.table(cellText=rows, colLabels=["Recorded batch", "Control", "PD", "Total", "PD fraction"],
                   loc="center", cellLoc="center", colColours=["#DCE8EB"] * 5,
                   colWidths=[.38, .13, .11, .12, .18])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    tbl.scale(1, 1.42)
    for i, row in enumerate(data, start=1):
        if int(row["control"]) <= 2 or int(row["PD"]) <= 2:
            for j in range(5):
                tbl[(i, j)].set_facecolor("#FBE6DD")
    ax.set_title("GENEPARK recorded batch × diagnosis", fontsize=12, pad=14)
    fig.tight_layout()
    save_figure(fig, out, "table_2_GENEPARK_batch")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ppmi-predictions", required=True, type=Path)
    parser.add_argument("--ppmi-primary-aggregate", required=True, type=Path)
    parser.add_argument("--ppmi-secondary-aggregate", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    a = parser.parse_args()
    private_root = (Path.home() / "PPMI_private").resolve()
    out = a.output_dir.resolve()
    if private_root not in out.parents or out.exists():
        raise ValueError("Use a fresh output directory inside ~/PPMI_private")
    primary = json.loads(a.ppmi_primary_aggregate.read_text(encoding="utf-8"))
    secondary = json.loads(a.ppmi_secondary_aggregate.read_text(encoding="utf-8"))
    if primary["status"] != "completed_after_sex_gate" or secondary["status"] != "complete":
        raise ValueError("PPMI aggregate completion status is missing")
    if primary["ordered_population_sha256"] != secondary["ordered_population_sha256"]:
        raise ValueError("PPMI aggregates refer to different populations")
    os.umask(0o077)
    out.mkdir(mode=0o700)
    roc_figure(a.ppmi_predictions, primary, out)
    forest_figure(out)
    reliability_figure(secondary, out)
    batch_table(out)
    for path in out.iterdir():
        path.chmod(0o600)
    print(json.dumps({"status": "complete", "figure_count": 4, "output_dir": str(out)}))


if __name__ == "__main__":
    main()
