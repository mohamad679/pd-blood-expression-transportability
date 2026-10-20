"""Inventory persistent lambda boundaries from completed GEO decision artifacts."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNS = (
    ("Original", "private_data/geo_runs_20261018", 13),
    ("X/Y exclusion", "private_data/geo_xy_sensitivity_full_20261019", 13),
    ("Age/sex residualization", "private_data/geo_age_sex_resid_full_20261019", 13),
    ("Panel-gene exclusion", "private_data/geo_panel_gene_exclusion_full_20261019", 13),
    ("Age/sex matching", "private_data/geo_age_sex_matching_full_20261019", 6),
    ("ComBat", "private_data/geo_combat_sensitivity_20261020_retry1", 39),
)
MODELS = ("M2", "M3", "M5")


def nice_scope(path: Path) -> str:
    return str(path).replace("__outer__", "/outer/").replace("_outer_", "/outer/").replace(
        "__transfer", "/transfer"
    ).replace("_transfer", "/transfer").replace("pooled__final", "pooled/final").replace(
        "pooled_final", "pooled/final"
    )


def main() -> None:
    rows = []
    summary = []
    for analysis, relative, expected_scope_count in RUNS:
        base = ROOT / relative
        files = sorted(base.rglob("decisions.json")) + sorted(base.rglob("M5_decision.json"))
        if len(files) != expected_scope_count:
            raise ValueError(f"{analysis}: expected {expected_scope_count} scopes, found {len(files)}")
        examined = Counter()
        flagged = Counter()
        for file in files:
            content = json.loads(file.read_text())
            decisions = {"M5": content} if file.name == "M5_decision.json" else content
            for model in MODELS:
                decision = decisions.get(model)
                if decision is None:
                    continue
                if "selected" not in decision or "persistent_boundary" not in decision:
                    raise ValueError(f"Incomplete decision: {file} {model}")
                examined[model] += 1
                if not decision["persistent_boundary"]:
                    continue
                selected = decision["selected"]
                if selected.get("baseline") or decision.get("expansion_rounds") != 2:
                    raise ValueError(f"Unexpected boundary state: {file} {model}")
                flagged[model] += 1
                scope = nice_scope(file.parent.relative_to(base))
                rows.append(
                    {
                        "analysis": analysis,
                        "scope": scope,
                        "model": model,
                        "lambda": selected["lambda"],
                        "rho": selected["rho"],
                        "expansion_rounds": decision["expansion_rounds"],
                    }
                )
        summary.append(
            {
                "analysis": analysis,
                "scope_count": len(files),
                "decisions_examined": dict(examined),
                "persistent_boundary_count": dict(flagged),
            }
        )
    payload = {
        "status": "complete",
        "purpose": "descriptive inventory of persistent lambda-boundary flags",
        "PPMI_accessed": False,
        "summary": summary,
        "selections": rows,
    }
    output_json = ROOT / "results/geo_boundary_selection_inventory_20261020.json"
    output_json.write_text(json.dumps(payload, indent=2) + "\n")

    lines = [
        "# Supplementary Table S1. Persistent lambda-boundary selections in GEO",
        "",
        "A boundary is listed only when the saved decision's `persistent_boundary` flag is true after two permitted grid expansions. These are descriptive tuning decisions, not performance estimates. All integrated final runs are included; timing pilots and interrupted copies are excluded. GSE57475 is PROBE, and GSE99039 is GENEPARK. A transfer scope names its training cohort.",
        "",
        "| Analysis | Scopes | M2 flagged / examined | M3 flagged / examined | M5 flagged / examined | All flagged / examined |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for item in summary:
        n = item["decisions_examined"]
        f = item["persistent_boundary_count"]
        total_n = sum(n.values())
        total_f = sum(f.values())
        lines.append(
            f"| {item['analysis']} | {item['scope_count']} | "
            + " | ".join(f"{f.get(model, 0)}/{n.get(model, 0)}" for model in MODELS)
            + f" | {total_f}/{total_n} |"
        )
    lines += [
        "",
        "| Analysis | Scope or ComBat arm / scope | Model | Lambda | Rho | Grid expansions |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['analysis']} | {row['scope']} | {row['model']} | "
            f"{row['lambda']:.0e} | {row['rho']:.1f} | {row['expansion_rounds']} |"
        )
    lines += [
        "",
        "The separate GEO sex-marker capacity control selected lambda=100 at the upper edge in both directional training sources after two extensions. It used a distinct auxiliary model and grid and is excluded from the M2/M3/M5 counts above.",
        "",
        "Source: saved final-run decision artifacts; extraction method in `docs/GEO_BOUNDARY_INVENTORY_METHOD_20261020_EN.md`. No PPMI data were accessed.",
        "",
    ]
    output_md = ROOT / "docs/GEO_SUPPLEMENTARY_BOUNDARY_TABLE_20261020_EN.md"
    output_md.write_text("\n".join(lines))
    print(f"{len(rows)} persistent boundaries across {sum(sum(s['decisions_examined'].values()) for s in summary)} model decisions")
    print(output_json)
    print(output_md)


if __name__ == "__main__":
    main()
