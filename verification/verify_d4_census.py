"""Compatibility entry for ak_diagram4_sign_census."""
import sys
from ak_diagram4_sign_census import main
if __name__=='__main__':
    try:main()
    except (ArithmeticError,AssertionError,ValueError,ZeroDivisionError) as exc:
        print(f'FAIL: {exc}',file=sys.stderr)
        raise SystemExit(1)
