"""Reconstruction of ak_diagram4_valuation_coverage from Table 5, paper page 56.

Only the checks reported by this module are implemented; see program_index.json.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from smale_verify.exact import (P, V, zero, equal, edge, equations, substitute,
    diff_poly, divide_poly, groebner, det, quotient_norm)

def verify_torque_expansion(D,m,pair,tri):
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    # D4 torque expansion in Lemma 8.4.  Represent t^(e/2) by an independent
    # Laurent variable r, so every affine exponent in e becomes an integer pair.
    t4,r4=V('D4_t'),V('D4_r_half_e')
    A4={i:V('D4_A'+str(i)) for i in pair};B4={i:V('D4_B'+str(i)) for i in pair}
    ZH4,WH4,u4,v4,omega4,qmass4=map(V,('D4_Zhat','D4_What','D4_u','D4_v','D4_omega','D4_q'))
    d4i={i:V('D4_d'+str(i)) for i in pair}
    zeta4={3:m[5]*u4,4:m[5]*v4,5:-m[3]*u4-m[4]*v4}
    eta4={3:r4**2*omega4,4:m[5]*qmass4,
           5:-m[3]/m[5]*r4**2*omega4-m[4]*qmass4}
    zd4={i:t4**-4*A4[i] for i in pair}
    wd4={i:t4**4*B4[i] for i in pair}
    for k in tri:
        zd4[k]=t4**-1*r4**-3*ZH4+t4**2*zeta4[k]
        wd4[k]=t4**5*r4**-1*WH4+t4**2*eta4[k]
    dd4={}
    dd4[(1,2)]=D[(1,2)]
    for i in pair:
        dd4[(i,3)]=t4*r4**-1*d4i[i]
        for k in (4,5):dd4[(i,k)]=t4*D[(i,k)]
    for e in ((3,4),(3,5),(4,5)):dd4[e]=t4**-2*D[e]
    TP4=sum((m[i]*m[k]*dd4[(i,k)]**3*(zd4[i]*wd4[k]-wd4[i]*zd4[k])
              for i in pair for k in tri),P())
    groups4={}
    for mon,coef in TP4.terms.items():
        ex=dict(mon);key=(ex.get('D4_t',0),ex.get('D4_r_half_e',0))
        rest=tuple((name,power) for name,power in mon if name not in {'D4_t','D4_r_half_e'})
        groups4[key]=groups4.get(key,P())+P({rest:coef})
    groups4={key:poly for key,poly in groups4.items() if poly.terms}
    lead4=(1,-1)
    lead4poly=m[3]*sum((m[i]*d4i[i]*(d4i[i]**2*A4[i]*omega4) for i in pair),P())
    equal(groups4.get(lead4,P()),lead4poly,'D4.8.4 leading torque coefficient')
    g4_pairs=[(0,1),(0,3),(3,0),(3,-3),(5,-5),(5,-2),(8,-2),(8,1)]
    expected4={lead4}|{(lead4[0]+a,lead4[1]+b) for a,b in g4_pairs}
    if set(groups4)!=expected4:
        raise AssertionError(('D4 torque exponents',set(groups4),expected4))
    record('D4.8.4.torque_expansion',9,
           {'leading_term':'m3 sum_i m_i d_i(d_i^2 A_i omega)',
            'remainder_exponents_in_(t,r)':[list(x) for x in sorted(expected4-{lead4})],
            'positive_for_0<e<2':True})
    return checks, {}

if __name__ == "__main__":
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram4_valuation_coverage'))
