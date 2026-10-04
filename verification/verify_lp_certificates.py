#!/usr/bin/env python3
"""Independent exact verifier for rational LP witnesses (D5).

Input JSON stores A x <= b, with every rational encoded as an integer or a
string such as "-7/12".  Equalities must be represented by two inequalities.
Strict inequalities are represented by adding a common slack variable s and
including 0 < s <= 1 among the constraints.

An infeasibility certificate y has y_i >= 0, y^T A = 0, y^T b < 0.
A feasibility witness is an exact rational vector x with A x <= b.
This verifies a certificate; it does not construct the D5 regime model.
"""
from fractions import Fraction
import argparse
import json
from pathlib import Path
import sys


def rational(value, where):
    if isinstance(value, bool) or isinstance(value, float):
        raise ValueError(f"{where}: use integers or rational strings, never floats")
    try:
        return Fraction(value)
    except (ValueError, TypeError, ZeroDivisionError) as exc:
        raise ValueError(f"{where}: invalid rational {value!r}") from exc


def load_system(obj):
    n = obj["variables"]
    if not isinstance(n, int) or n < 1:
        raise ValueError("variables must be a positive integer")
    raw_a, raw_b = obj["A"], obj["b"]
    if len(raw_a) != len(raw_b):
        raise ValueError("A and b row counts differ")
    A = []
    b = []
    for i, (row, rhs) in enumerate(zip(raw_a, raw_b)):
        if len(row) != n:
            raise ValueError(f"A row {i} has wrong width")
        A.append([rational(x, f"A[{i}]") for x in row])
        b.append(rational(rhs, f"b[{i}]"))
    return n, A, b


def check_witness(n, A, b, x):
    if len(x) != n:
        raise ValueError("witness has wrong length")
    residuals = [sum(a * v for a, v in zip(row, x)) - rhs for row, rhs in zip(A, b)]
    if any(r > 0 for r in residuals):
        bad = [(i, str(r)) for i, r in enumerate(residuals) if r > 0]
        raise ValueError(f"feasibility witness violates rows {bad}")
    return residuals


def check_farkas(A, b, y):
    if len(y) != len(A):
        raise ValueError("Farkas vector has wrong length")
    if any(v < 0 for v in y):
        raise ValueError("Farkas multipliers must be nonnegative")
    n = len(A[0])
    combined = [sum(y[i] * A[i][j] for i in range(len(A))) for j in range(n)]
    rhs = sum(v * r for v, r in zip(y, b))
    if any(combined):
        raise ValueError(f"y^T A is not zero: {[str(v) for v in combined]}")
    if rhs >= 0:
        raise ValueError(f"Farkas contradiction needs y^T b < 0; got {rhs}")
    return rhs


def check_case(case):
    n, A, b = load_system(case)
    mode = case.get("kind")
    if mode == "feasible":
        x = [rational(v, f"x[{i}]") for i, v in enumerate(case["x"])]
        residuals = check_witness(n, A, b, x)
        return {"case": case.get("id", "unnamed"), "status": "PASS",
                "kind": mode, "variables": n, "constraints": len(A),
                "witness": [str(v) for v in x],
                "maximum_slack": str(max([-r for r in residuals], default=Fraction(0)))}
    if mode == "infeasible":
        y = [rational(v, f"y[{i}]") for i, v in enumerate(case["y"])]
        rhs = check_farkas(A, b, y)
        return {"case": case.get("id", "unnamed"), "status": "PASS",
                "kind": mode, "constraints": len(A), "certificate_rhs": str(rhs)}
    raise ValueError("kind must be feasible or infeasible")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path, nargs="?", help="JSON certificate file")
    parser.add_argument("--stdin", action="store_true", help="read one JSON certificate from standard input")
    args = parser.parse_args()
    if bool(args.certificate) == args.stdin:
        parser.error("provide exactly one certificate path or --stdin")
    obj = json.load(sys.stdin) if args.stdin else json.loads(args.certificate.read_text(encoding="utf-8"))
    cases = obj if isinstance(obj, list) else obj.get("cases", [obj])
    results = [check_case(case) for case in cases]
    print(json.dumps({"checker": "D5 exact rational LP certificates",
                      "certificate_file": "<stdin>" if args.stdin else args.certificate.name,
                      "case_count": len(results), "results": results,
                      "warning": "PASS certifies only supplied LP rows and witnesses, not that the rows exhaust the manuscript's D5 branches."},
                     indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        raise SystemExit(1)
