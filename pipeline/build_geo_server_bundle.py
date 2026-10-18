"""Build a GEO-only transfer archive from an explicit allowlist, never the workspace."""
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    names = ['requirements_glum_v3.lock.txt', 'docs/HETZNER_GEO_RUNBOOK_EN.md',
             'protocol/config_v3_1.json', 'protocol/v3_marker_manifest.json',
             'protocol/permutation_config_cloud_v1.json',
             'protocol/numerical_kkt_v3_1.py', 'protocol/block_logistic_v3_1.py',
             'data/manifests/input_audit.json']
    names += ['pipeline/' + name for name in (
        'geo_pipeline.py', 'permutation_executor.py', 'verify_permutation_executor.py',
        'finalize_permutation_preflight.py', 'server_setup.sh', 'server_stage.sh',
        'test_server_portability.py')]
    expected = json.loads((ROOT/'data/manifests/input_audit.json').read_text())['input_sha256']
    for name, digest in expected.items():
        path = ROOT/'prepared'/name
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen GEO input mismatch: ' + name)
        names.append('prepared/' + name)
    # Identity needs only original ten-fold predictions and tuning decisions.
    # No PPMI paths, tables, models, caches, drafts or credentials are included.
    for accession in ('GSE57475', 'GSE99039'):
        for fold in range(5):
            for file in ('heldout_predictions.npz', 'decisions.json'):
                names.append(f'private_data/geo_runs_20261018/{accession}_outer_{fold}/{file}')
    output = ROOT/'private_data/server_transfer/GEO_SERVER_BUNDLE.tar.gz'
    output.parent.mkdir(parents=True, exist_ok=True)
    for name in names:
        path = ROOT/name
        if path.is_symlink() or not path.is_file():
            raise ValueError('Missing or symlinked allowlisted file: ' + name)
    manifest = ''.join(f'{hashlib.sha256((ROOT/name).read_bytes()).hexdigest()}  {name}\n'
                       for name in names)
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    provenance = dict(source_commit=revision, files=len(names),
                      all_files_individually_hashed=True, PPMI_files_included=False,
                      note='Source commit plus TRANSFER_SHA256.txt identifies the exact package.')
    with tempfile.TemporaryDirectory() as work:
        temporary = Path(work)
        (temporary/'TRANSFER_SHA256.txt').write_text(manifest)
        (temporary/'TRANSFER_PROVENANCE.json').write_text(json.dumps(provenance, indent=2)+'\n')
        with tarfile.open(output, 'w:gz') as archive:
            for name in names:
                archive.add(ROOT/name, arcname='pdgeo/'+name, recursive=False)
            for name in ('TRANSFER_SHA256.txt', 'TRANSFER_PROVENANCE.json'):
                archive.add(temporary/name, arcname='pdgeo/'+name)
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    output.with_suffix(output.suffix+'.sha256').write_text(f'{digest}  {output.name}\n')
    print(json.dumps(dict(path=str(output), bytes=output.stat().st_size,
                         sha256=digest, files=len(names), PPMI_files_included=False), indent=2))


if __name__ == '__main__':
    main()
