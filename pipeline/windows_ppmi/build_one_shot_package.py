#!/usr/bin/env python3
"""Build the deterministic, code-only Windows one-shot tool ZIP."""

import hashlib
import json
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
OUTPUT = REPO / "tools" / "PPMI_ONE_SHOT_TOOLS_20261020.zip"
SOURCES = {
    "prepare_ppmi_remaining_inputs.py": HERE / "prepare_ppmi_remaining_inputs.py",
    "audit_ppmi_preqc_qc.py": REPO / "pipeline" / "audit_ppmi_preqc_qc.py",
    "audit_ir3_all_once.py": HERE / "audit_ir3_all_once.py",
    "prepare_ppmi_windows_bundle.py": HERE / "prepare_ppmi_windows_bundle.py",
    "recover_ir3_links_from_cache.py": HERE / "recover_ir3_links_from_cache.py",
    "README.txt": HERE / "README_PPMI_ONE_SHOT_EXPORT.txt",
}


def write_member(archive: zipfile.ZipFile, name: str, blob: bytes) -> None:
    info = zipfile.ZipInfo(name, date_time=(2026, 10, 20, 0, 0, 0))
    info.compress_type = zipfile.ZIP_DEFLATED
    info.external_attr = 0o644 << 16
    archive.writestr(info, blob, compresslevel=9)


def main() -> None:
    contents = {name: path.read_bytes() for name, path in SOURCES.items()}
    manifest = {
        "schema": "PPMI_one_shot_code_package_v1",
        "contains_ppmi_data": False,
        "files_sha256": {name: hashlib.sha256(blob).hexdigest()
                         for name, blob in sorted(contents.items())},
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(OUTPUT, "w") as archive:
        for name, blob in sorted(contents.items()):
            write_member(archive, name, blob)
        write_member(archive, "MANIFEST.json",
                     (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode())
    with zipfile.ZipFile(OUTPUT) as archive:
        assert archive.testzip() is None
    print(str(OUTPUT))
    print("sha256 " + hashlib.sha256(OUTPUT.read_bytes()).hexdigest())
    print("bytes " + str(OUTPUT.stat().st_size))


if __name__ == "__main__":
    main()
