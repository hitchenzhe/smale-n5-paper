#!/usr/bin/env python3
"""Compatibility entry point for the paper-named exact algebra programs."""
from smale_verify.runner import algebra_cli
if __name__ == '__main__':
    raise SystemExit(algebra_cli())
