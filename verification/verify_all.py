#!/usr/bin/env python3
"""Regenerate implemented checks and bind claim status to actual execution.

Exit 0: all declared computational claims covered. Exit 1: failure or bad
source/binding. Exit 2: implemented checks passed; proof obligations remain.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from smale_verify.runner import run_exact_checks
from verify_v19_changes import source_audit


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def save(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path,
        help='source root containing the files listed in claim_status.json; default: package parent')
    parser.add_argument('--output-dir', type=Path,
        help='execution records directory; default: verification/certificates')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    source_root = (args.source_root or here.parent).resolve()
    output = (args.output_dir or here/'v19-records').resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((here/'claim_status.json').read_text(encoding='utf-8'))
    code_files = sorted([*here.glob('*.py'), *here.glob('smale_verify/*.py'),
        here/'smale_verify/stages.json', here/'program_index.json', here/'claim_status.json',*here.glob('vendor/*.whl')])
    code_before = {str(path.relative_to(here)).replace('\\','/'):digest(path) for path in code_files}
    failures = []
    if set(manifest.get('claim_classes',{})) != set(manifest['claims']):
        failures.append('claim-class map does not cover every declared group')
    source_records = {}
    for name, expected in manifest['source_files'].items():
        path = source_root/name
        actual = digest(path) if path.is_file() else None
        source_records[name] = dict(expected=expected, actual=actual,
                                    matches=actual is not None and actual == expected.upper())
        if not source_records[name]['matches']:
            failures.append('source missing or changed: '+name)
    print('Source revision: '+('PASS' if not failures else 'FAIL'), flush=True)
    binding = {}
    if not failures:
        try:
            binding = source_audit(source_root)
            if binding['outdated_active_code_labels']:
                failures.append('program index contains removed source labels')
            if any(value == 'UNMAPPED' for value in binding['documented_replacements'].values()):
                failures.append('an appendix program has neither an implementation nor a documented replacement')
        except Exception as exc:
            failures.append('source/program binding: '+str(exc))
    if failures:
        save(output/'verification_run.json',dict(status='FAIL',exit_code=1,
             source_revision=source_records,source_program_binding=binding,failures=failures,
             scope='Preflight failure; computational checks were not run.'))
        for failure in failures:
            print('FAIL '+failure,file=sys.stderr)
        return 1
    checks = {}
    print('Regenerating exact algebra and modular jet records...', flush=True)
    try:
        algebra = run_exact_checks()
        for record in algebra['checks']:
            name = record['check']
            if name in checks:
                raise ValueError('duplicate check identifier: '+name)
            checks[name] = record
            if record['status'] != 'PASS':
                failures.append('check not PASS: '+name)
        save(output/'exact_algebra_record.json', algebra)
        program_directory = output/'programs'
        program_directory.mkdir(parents=True, exist_ok=True)
        for program in sorted({record['program'] for record in algebra['checks']}):
            save(program_directory/(program+'.json'), dict(
                python=algebra['python'], arithmetic=algebra['arithmetic'],
                program=program, checks=[record for record in algebra['checks']
                    if record['program']==program], scope=algebra['scope']))
        print(f"Exact checks: {len(checks)} PASS", flush=True)
    except Exception as exc:
        failures.append('exact algebra: '+str(exc))
        print('Exact checks: FAIL '+str(exc), flush=True)
    print('Regenerating D4 rational Sturm census...', flush=True)
    d4_path = output/'d4_lemma_8_10.json'
    process = subprocess.run([sys.executable, str(here/'ak_diagram4_sign_census.py'),
        '--output', str(d4_path), '--summary'], capture_output=True, text=True,
        encoding='utf-8', errors='replace', check=False)
    if process.returncode:
        failures.append('D4 census: '+process.stderr.strip())
    else:
        try:
            d4 = json.loads(d4_path.read_text(encoding='utf-8'))
            expected_counts = dict(vectors=64, with_unique_positive_pair=16,
                                   without_positive_pair=48, second_obstruction_checks=8)
            if any(d4.get(k) != v for k, v in expected_counts.items()):
                raise ValueError('unexpected D4 census counts')
            checks['D4.8.10.sign_census'] = dict(check='D4.8.10.sign_census',
                status='PASS', program='ak_diagram4_sign_census', function='main',
                certificate='d4_lemma_8_10.json', **expected_counts)
            print('D4 census: PASS (64 signs, 16 candidates, 48 exclusions)', flush=True)
        except Exception as exc:
            failures.append('D4 record: '+str(exc))
    claims = {}
    for name, description in manifest['claims'].items():
        required = manifest.get('required_checks', {}).get(name, [])
        if description.startswith('NOT APPLICABLE:'):
            status = 'NOT_APPLICABLE'
        elif description.startswith('PASS:'):
            missing = [key for key in required if checks.get(key, {}).get('status') != 'PASS']
            if not required or missing:
                status = 'FAIL'
                failures.append('claim binding failed: '+name+' '+str(missing or ['no required checks']))
            else:
                status = 'PASS'
        else:
            status = 'INCOMPLETE'
        claims[name] = dict(status=status, description=description, required_checks=required)
    incomplete = [name for name, value in claims.items() if value['status'] == 'INCOMPLETE']
    summary = {status:sum(value['status'] == status for value in claims.values())
               for status in ('PASS','NOT_APPLICABLE','INCOMPLETE','FAIL')}
    code_after = {str(path.relative_to(here)).replace('\\','/'):digest(path) for path in code_files}
    if code_after != code_before:
        failures.append('verification code or manifest changed during execution')
    for name, expected in manifest['source_files'].items():
        if not (source_root/name).is_file() or digest(source_root/name) != expected.upper():
            failures.append('source changed during execution: '+name)
    exit_code = 1 if failures else 2 if incomplete else 0
    record = dict(status='FAIL' if failures else 'INCOMPLETE' if incomplete else 'PASS',
        exit_code=exit_code, python=sys.version.split()[0], source_root=str(source_root),
        arithmetic='exact rational polynomials, integer Sylvester determinants, rational Bernstein covers, algebraic-number jets, Sturm chains and finite fields; bundled SymPy 1.14.0; no floating proof steps',
        source_revision=source_records,
        source_program_binding=binding,
        code_sha256=code_before,
        code_and_source_unchanged_during_execution=code_after == code_before and all(
            (source_root/name).is_file() and digest(source_root/name) == expected.upper()
            for name, expected in manifest['source_files'].items()),
        regenerated_certificate_sha256={path.name:digest(path)
            for path in (output/'exact_algebra_record.json',d4_path) if path.is_file()},
        mathematical_certificate_sha256={path.name:digest(path) for path in (here/'certificates').glob('*.json')
            if path.name in ('d2_mass_resultants.json','d2_positive_mass_exclusions.json','d2_L4_same_sign.json','d2_K5_mass_ray.json','d2_mass_rays.json','d2_mass_ray_root_inputs.json','d2_ray_response.json','d2_ray_response_inputs.json','d6_leading_jets.json','d5_branch_index.json','d5_elementary_four.json','d2_ordinary_seed_inputs.json','d2_ordinary_seed.json','d2_fold_seed_inputs.json','d2_fold_seed.json','d2_equal_mass_orbits.json','d2_S5_exact_response.json','d2_general_fold_block.json','d2_regular_X.json','d2_Euler_global.json','d2_S6_S8_global.json','d2_K1_branches.json','d2_K4_arcs.json','d2_K5_all_labels.json','d5_source_audit.json','d5_census_inputs.json','d5_census.json','d5_dominant_inputs.json','d5_dominant.json','d5_inner_witnesses.json','d5_complete_four_products.json','d8_ray_diagrams.json','d8_pattern_boundaries.json')},
        executed_checks=sorted(checks), claim_summary=summary, claims=claims,
        failures=failures,
        scope='Computational coverage only, including documented alternative certificates. Analytic arguments and cited external results are not formally verified. This record does not certify an unconditional main theorem.')
    save(output/'verification_run.json',record)
    print('Claim coverage: '+json.dumps(summary))
    for name in incomplete:
        print('INCOMPLETE '+name+': '+claims[name]['description'])
    for failure in failures:
        print('FAIL '+failure, file=sys.stderr)
    print('Overall: '+record['status']+'; execution record: '+str(output/'verification_run.json'))
    return exit_code


if __name__ == '__main__':
    raise SystemExit(main())
