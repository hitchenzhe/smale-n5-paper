#!/usr/bin/env python3
"""Check that the LaTeX source tree matches verification/claim_status.json."""
import argparse
import hashlib
import json
from pathlib import Path
import sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root",type=Path,required=True,help="folder containing the source files listed in claim_status.json")
    args=parser.parse_args()
    manifest=Path(__file__).with_name("claim_status.json")
    expected=json.loads(manifest.read_text(encoding="utf-8"))["source_files"]
    errors=[]
    for relative,want in expected.items():
        path=args.source_root/Path(relative)
        if not path.is_file():
            errors.append(f"MISSING {relative}")
            continue
        got=hashlib.sha256(path.read_bytes()).hexdigest().upper()
        if got!=want.upper():errors.append(f"CHANGED {relative}: expected {want}, got {got}")
        else:print(f"OK {relative} {got}")
    if errors:
        print("\n".join(errors),file=sys.stderr)
        return 1
    print(f"PASS: all {len(expected)} source files match the recorded revision.")
    return 0


if __name__=="__main__":raise SystemExit(main())


