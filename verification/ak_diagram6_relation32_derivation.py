"""Reconstruction of ak_diagram6_relation32_derivation from Table 5, paper page 56.

Only the checks reported by this module are implemented; see program_index.json.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from smale_verify.exact import (P, V, zero, equal, edge, equations, substitute,
    diff_poly, divide_poly, groebner, det, quotient_norm)

def verify_torque_identity(D,W,Z,cross,m,rw,rz,t):
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    Sigma=sum(m[i]*m[j]*D[(i,j)]**3*Z[i]*W[j] for i,j in cross)
    lhs=sum(m[i]*(W[i]*rz[i]-Z[i]*rw[i]) for i in (1,2))+m[5]*W[5]*rz[5]
    Q=m[5]*(W[5]*Z[5]*(1-sum(m[k]*D[edge(k,5)]**3 for k in (1,2,3,4)))+
            sum(m[k]*D[edge(k,5)]**3*W[k]*Z[5] for k in (1,2))+
            sum(m[k]*D[edge(k,5)]**3*W[5]*Z[k] for k in (3,4)))+t**3*sum(m[i]*m[j]*D[(i,j)]**3*W[i]*Z[j] for i,j in cross)
    equal(lhs,-t*Sigma+t**2*Q,'D6 9.4 explicit Q')
    record('D6.9.4.torque_explicit_Q',1,{'Q_terms':len(Q.terms),'Q_t_degrees':Q.exponents('t')})
    return checks, {}

if __name__ == "__main__":
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram6_relation32_derivation'))
