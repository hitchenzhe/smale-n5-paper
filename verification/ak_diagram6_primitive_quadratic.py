"""Reconstruction of ak_diagram6_primitive_quadratic from Table 5, paper page 56.

Only the checks reported by this module are implemented; see program_index.json.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from smale_verify.exact import (P, V, zero, equal, edge, equations, substitute,
    diff_poly, divide_poly, groebner, det, quotient_norm)

def verify_same_shape_norm():
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    k=V('k')
    equal(quotient_norm(1-2*k**3),1089,'D6 primitive same norm')
    record('D6.9.9.same_shape_norm_only',1,{'scope':'norm only; jet derivation and ranks remain'})
    return checks, {}

if __name__ == "__main__":
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram6_primitive_quadratic'))
