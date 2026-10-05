"""Run paper-named exact checks with explicit dependency inputs."""
import argparse
import importlib
import json
from pathlib import Path
import sys
from .exact import normalized_system

ROOT = Path(__file__).resolve().parents[1]

def run_exact_checks(programs=None):
    stages = json.loads((ROOT / 'smale_verify' / 'stages.json').read_text(encoding='utf-8'))
    context = normalized_system()
    selected = set(range(len(stages))) if programs is None else {i for i, stage in enumerate(stages) if stage['program'] in programs}
    if programs is not None:
        missing = set(programs) - {s['program'] for s in stages}
        if missing: raise ValueError('No implemented exact program: ' + ', '.join(sorted(missing)))
    pending = list(selected)
    while pending:
        i = pending.pop()
        for name in stages[i]['deps']:
            if name in context: continue
            producer = next((j for j in range(i-1,-1,-1) if name in stages[j]['exports']), None)
            if producer is None: raise RuntimeError(f"No producer for {stages[i]['program']}.{stages[i]['name']}: {name}")
            if producer not in selected:
                selected.add(producer)
                pending.append(producer)
    checks = []
    for i in sorted(selected):
        stage = stages[i]
        function = getattr(importlib.import_module(stage['program']), stage['name'])
        records, outputs = function(**{name: context[name] for name in stage['deps']})
        context.update(outputs)
        for record in records:
            checks.append(dict(record, program=stage['program'], function=stage['name'], prerequisite=programs is not None and stage['program'] not in programs))
    return dict(python=sys.version.split()[0], arithmetic='exact rational sparse Laurent polynomials, finite fields, bounded integer resultants, rational Bernstein covers and algebraic-number jets; bundled SymPy 1.14.0', checks=checks, scope='implemented computational checks only; see claim_status.json for unverified claims')

def algebra_cli(argv=None, default_program=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--program', action='append', help='paper program name, including ak_')
    parser.add_argument('--output', type=Path, help='save a reproducible execution record as JSON')
    parser.add_argument('--list', action='store_true', help='list all paper program names and coverage')
    args = parser.parse_args(argv)
    try:
        if args.list:
            print((ROOT / 'program_index.json').read_text(encoding='utf-8'))
            return 0
        programs = args.program or ([default_program] if default_program else None)
        result = run_exact_checks(programs)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
        print(json.dumps(result, indent=2))
        return 0
    except (AssertionError, ArithmeticError, ValueError, RuntimeError) as exc:
        print(f'FAIL: {exc}', file=sys.stderr)
        return 1

def program_cli(program):
    return algebra_cli(default_program=program)
