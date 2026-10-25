"""Structural version comparison without copying private manuscript prose."""

import csv
import difflib
import hashlib
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRIVATE = Path.home() / "PPMI_private"
ORIGINAL = PRIVATE / "manuscript_20261021/MANUSCRIPT_PRIVATE.md"
REVISED = PRIVATE / "manuscript_20261021_revision1/MANUSCRIPT_PRIVATE_REVISED.md"
OUT = ROOT / "revision/outputs/manuscript_version_inventory_P1.csv"


def describe(path):
    text = path.read_text()
    lines = text.splitlines()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "title": next((line[2:] for line in lines if line.startswith("# ")), "<absent>"),
            "word_count": len(text.split()),
            "h2_sections": len([line for line in lines if line.startswith("## ")]),
            "h3_sections": len([line for line in lines if line.startswith("### ")]),
            "table_captions": len(re.findall(r"\*\*Table\s+\d+", text)),
            "figure_legends": len(re.findall(r"\*\*Figure\s+\d+", text))}


def main():
    a, b = describe(ORIGINAL), describe(REVISED)
    a_lines, b_lines = ORIGINAL.read_text().splitlines(), REVISED.read_text().splitlines()
    matcher = difflib.SequenceMatcher(None, a_lines, b_lines, autojunk=False)
    changed = sum(max(i2-i1, j2-j1) for op, i1, i2, j1, j2 in matcher.get_opcodes() if op != "equal")
    a["approx_changed_lines_vs_original"] = 0
    b["approx_changed_lines_vs_original"] = changed
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(a))
        writer.writeheader()
        writer.writerows([a,b])
    print("Manuscript structural inventory complete; current working candidate: revision1")


if __name__ == "__main__":
    main()
