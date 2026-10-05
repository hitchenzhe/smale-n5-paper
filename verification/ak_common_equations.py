"""Reconstruction of ak_common_equations from Table 5, paper page 56.

Only the checks reported by this module are implemented; see program_index.json.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from smale_verify.exact import (P, V, zero, equal, edge, equations, substitute,
    diff_poly, divide_poly, groebner, det, quotient_norm)

def verify_common_equations(D,Ew,Ez,m,ties,w,z):
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    equal(sum(m[j]*Ez[j] for j in m),sum(m[j]*z[j] for j in m),'mass center z')
    equal(sum(m[j]*Ew[j] for j in m),sum(m[j]*w[j] for j in m),'mass center w')
    zero(sum(m[j]*(w[j]*Ez[j]-z[j]*Ew[j]) for j in m),'total torque')
    I=sum(m[j]*z[j]*w[j] for j in m)
    U=sum(m[i]*m[j]*D[(i,j)] for i,j in D)
    vir=sum(m[j]*w[j]*Ez[j] for j in m)+sum(m[i]*m[j]*D[(i,j)]*ties[(i,j)] for i,j in D)
    equal(I-U,vir,'I-U ideal certificate')
    record('common.mass_center_torque_virial',4)
    return checks, {}

if __name__ == "__main__":
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_common_equations'))
