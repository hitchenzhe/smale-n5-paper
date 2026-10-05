"""Reconstruction of ak_diagram1_virial_exclusion from Table 5, paper page 56.

Only the checks reported by this module are implemented; see program_index.json.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from smale_verify.exact import (P, V, zero, equal, edge, equations, substitute,
    diff_poly, divide_poly, groebner, det, quotient_norm)

def verify_triangle_virial(m):
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    tz={k:V('tz'+str(k)) for k in (3,4,5)}
    tw={k:V('tw'+str(k)) for k in (3,4,5)}
    tD={edge(i,j):V('tD%d%d'%(i,j)) for i,j in combinations((3,4,5),2)}
    tphi={k:sum(m[l]*tD[edge(k,l)]**3*(tz[k]-tz[l]) for l in (3,4,5) if k!=l) for k in (3,4,5)}
    tpsi={k:sum(m[l]*tD[edge(k,l)]**3*(tw[k]-tw[l]) for l in (3,4,5) if k!=l) for k in (3,4,5)}
    tchi={e:tD[e]**2*(tz[e[1]]-tz[e[0]])*(tw[e[1]]-tw[e[0]])-1 for e in tD}
    tU=sum(m[i]*m[j]*tD[(i,j)] for i,j in tD)
    equal(sum(m[k]*tz[k]*tpsi[k] for k in (3,4,5)),
          tU+sum(m[i]*m[j]*tD[(i,j)]*tchi[(i,j)] for i,j in tD),'D1.5.11 virial identity')
    for name in sorted({'tz3','tz4','tz5','tw3','tw4','tw5','tD34','tD35','tD45'}):
        lhs=2*diff_poly(tU,name)+sum(m[k]*(tpsi[k]*diff_poly(tz[k],name)+tphi[k]*diff_poly(tw[k],name)) for k in (3,4,5))
        rhs=sum(m[i]*m[j]*(tD[(i,j)]*diff_poly(tchi[(i,j)],name)-2*tchi[(i,j)]*diff_poly(tD[(i,j)],name)) for i,j in tD)
        equal(lhs,rhs,'D1.5.11 differential coefficient '+name)
    record('D1.5.11.triangle_virial_and_differential',10)
    return checks, {}

if __name__ == "__main__":
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram1_virial_exclusion'))
