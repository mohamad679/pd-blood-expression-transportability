"""Produce a public, aggregate preflight certificate only after all gates pass."""
import argparse
import json
from pathlib import Path

try:
    from . import permutation_executor as pe
except ImportError:
    import permutation_executor as pe


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--synthetic', type=Path, required=True)
    ap.add_argument('--run-root', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    synthetic = json.loads(args.synthetic.read_text())
    identity_path = args.run_root / 'IDENTITY_VERIFIED.json'
    identity = json.loads(identity_path.read_text())
    binding = json.loads((args.run_root / 'RUN_BINDING.json').read_text())
    if (synthetic.get('status') != 'synthetic_passed_identity_pending'
            or synthetic.get('memory_benchmark_completed') is not True
            or synthetic['code_sha256'] != pe.code_hashes()
            or synthetic['verifier_sha256'] != pe.sha_file(pe.ROOT / 'pipeline/verify_permutation_executor.py')
            or synthetic['config_sha256'] != pe.sha_file(pe.CONFIG)
            or identity.get('status') != 'passed'
            or identity['binding_sha256'] != pe.object_hash(binding)
            or binding['code_sha256'] != pe.code_hashes()
            or binding['config_sha256'] != pe.sha_file(pe.CONFIG)
            or binding['package_versions'] != pe.locked_environment()
            or len(identity['comparisons']) != 20):
        raise ValueError('Preflight evidence does not satisfy all release gates')
    pe.atomic_json(args.output, dict(status='passed', code_sha256=pe.code_hashes(),
        config_sha256=pe.sha_file(pe.CONFIG), verifier_sha256=synthetic['verifier_sha256'],
        synthetic_report_sha256=pe.sha_file(args.synthetic),
        observed_identity_report_sha256=pe.sha_file(identity_path),
        binding_sha256=pe.object_hash(binding), synthetic_checks=synthetic['checks'],
        synthetic_fit_attempts=synthetic['synthetic_fit_attempts'], max_KKT=synthetic['max_KKT'],
        resource_benchmark=synthetic['resource_benchmark'],
        observed_identity_comparisons=identity['comparisons'],
        real_randomized_fits_performed=False, PPMI_accessed=False,
        public_amendment_still_required=True))
    print('All preflight gates passed; publish the bound amendment before real permutations.')


if __name__ == '__main__':
    main()
