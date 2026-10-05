"""Reconstruction of ak_diagram6_primitive_opposite from Table 5, paper page 56.

Only the checks reported by this module are implemented; see program_index.json.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from smale_verify.exact import (P, V, zero, equal, edge, equations, substitute,
    diff_poly, divide_poly, groebner, det, quotient_norm)

def verify_opposite_norms():
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    k=V('k'); BB=V('B')
    equal(quotient_norm(2*BB+k**3),16*(BB**2+2)**2,'D6 primitive opposite denominator norm')
    equal(quotient_norm(2*BB),16*BB**4,'D6 primitive opposite numerator norm')
    record('D6.9.9.opposite_norms_only',2,{'scope':'norms only; jet derivation and ranks remain'})
    return checks, {}

if __name__ == "__main__":
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram6_primitive_opposite'))
