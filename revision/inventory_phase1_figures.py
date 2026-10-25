"""Inventory file identities only; never read or export figure pixels to stdout."""

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRIVATE = Path.home() / "PPMI_private"
OUT = ROOT / "revision/outputs/figure_inventory_P1.csv"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    figures = sorted((PRIVATE / "manuscript_figures_20261021_retry1").glob("*"))
    supplements = sorted((PRIVATE / "manuscript_20261021_revision1").glob("ONLINE_RESOURCE*"))
    rows = [{"kind": "figure" if path in figures else "online_resource", "path": str(path),
             "bytes": path.stat().st_size, "sha256": digest(path)}
            for path in figures + supplements if path.is_file()]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Inventory complete: {len(figures)} figure files and {len(supplements)} online resources")


if __name__ == "__main__":
    main()
