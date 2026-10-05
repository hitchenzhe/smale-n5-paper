"""Reconstruction of ak_diagram6_primitive_cubic from Table 5, paper page 56.

Only the checks reported by this module are implemented; see program_index.json.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from smale_verify.exact import (P, V, zero, equal, edge, equations, substitute,
    diff_poly, divide_poly, groebner, det, quotient_norm)

def verify_mixed_norm():
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    k=V('k'); AA=V('A'); BB=V('B')
    pol=AA*k**3-4*AA*k+BB*k**3-BB*k**2-4*BB*k+BB
    equal(quotient_norm(pol),9*(8*AA**2+16*AA*BB+9*BB**2)**2,'D6 primitive mixed norm')
    record('D6.9.9.mixed_norm_only',1,{'scope':'norm only; jet derivation and ranks remain'})
    return checks, {}

if __name__ == "__main__":
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram6_primitive_cubic'))
