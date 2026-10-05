"""Load the bundled, unmodified pure-Python exact symbolic dependencies."""
from pathlib import Path
import hashlib
import sys


def sympy():
    vendor=Path(__file__).resolve().parents[1]/'vendor'
    expected={'sympy-1.14.0-py3-none-any.whl':'e091cc3e99d2141a0ba2847328f5479b05d94a6635cb96148ccb3f34671bd8f5',
              'mpmath-1.3.0-py3-none-any.whl':'a0b2b9fe80bbcd81a6647ff13108738cfb482d481d826cc0e02f5b35e5c88d2c'}
    for name in ('sympy-1.14.0-py3-none-any.whl','mpmath-1.3.0-py3-none-any.whl'):
        path=vendor/name
        if not path.is_file():raise RuntimeError('Missing bundled exact symbolic dependency: '+name)
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected[name]:raise RuntimeError('Bundled symbolic dependency hash mismatch: '+name)
        if str(path) not in sys.path:sys.path.insert(0,str(path))
    import sympy as library
    if library.__version__!='1.14.0':raise RuntimeError('Expected SymPy 1.14.0')
    return library
