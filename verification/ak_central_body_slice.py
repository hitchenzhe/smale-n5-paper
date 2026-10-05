"""Reconstruction of ak_central_body_slice from Table 5, paper page 56.

Only the checks reported by this module are implemented; see program_index.json.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from smale_verify.exact import (P, V, zero, equal, edge, equations, substitute,
    diff_poly, divide_poly, groebner, det, quotient_norm)

def verify_torque_expansion(A,B,D,m,pair,q,solve_square,tri,wh,zet,zh):
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    # D1 Lemma 5.4: expand the reduced pair torque on the exceptional o_3=0 chart.
    tv5=('L','M','R','S','F')
    def T5(v): return P({tuple((x,n) for x,n in zip(tv5,v) if n):F(1)})
    theta,omega=V('vartheta'),V('omega')
    eta5={3:T5((0,0,0,0,2))*omega,4:m[5]*theta,
          5:-m[3]/m[5]*T5((0,0,0,0,2))*omega-m[4]*theta}
    z5={i:T5((0,2,2,0,0))*zh+T5((-2,0,0,0,0))*A[i] for i in pair}
    w5={i:T5((0,0,0,2,0))*wh+T5((2,0,0,0,0))*B[i] for i in pair}
    for k in tri:
        z5[k]=-q*T5((0,2,2,0,0))*zh+T5((0,2,0,0,0))*zet[k]
        w5[k]=-q*T5((0,0,0,2,0))*wh+eta5[k]
    d5torque={}
    for i in pair:
        d5torque[(i,3)]=T5((1,0,0,0,-1))*V('d'+str(i))
        for k in (4,5):d5torque[(i,k)]=T5((1,0,0,0,0))*D[(i,k)]
    for e in [(3,4),(3,5),(4,5)]:d5torque[e]=T5((0,-1,0,0,0))*D[e]
    torque=sum(m[i]*m[k]*d5torque[(i,k)]**3*(z5[i]*w5[k]-w5[i]*z5[k]) for i in pair for k in tri)
    torque_groups={}
    for mon,coef in torque.terms.items():
        powers=dict(mon);exp5=tuple(powers.get(x,0) for x in tv5)
        rest=tuple((x,n) for x,n in mon if x not in tv5)
        torque_groups[exp5]=torque_groups.get(exp5,P())+P({rest:coef})
    torque_groups={e:c for e,c in torque_groups.items() if c.terms}
    lead_exp=(1,0,0,0,-1)
    lead_poly=m[3]*sum(m[i]*V('d'+str(i))*(V('d'+str(i))**2*A[i]*omega) for i in pair)
    equal(torque_groups[lead_exp],lead_poly,'D1.5.4 leading torque coefficient')
    g5=[(0,0,0,0,1),(0,0,0,0,3),(0,0,0,2,-2),(0,0,0,2,1),
        (4,2,0,0,-2),(4,2,0,0,1),(2,2,2,0,0),(2,2,0,2,-2),
        (4,2,2,0,-2),(2,2,2,0,1),(2,2,0,2,1),(4,2,2,0,1),(2,2,2,0,3)]
    observed={tuple(a-b for a,b in zip(exponent,lead_exp)) for exponent in torque_groups if exponent!=lead_exp}
    if not observed<=set(g5):raise AssertionError(('D1 torque: unexpected remainder exponents',observed-set(g5)))
    if len(observed)!=13:raise AssertionError(('D1 torque: expected all 13 exponent classes',observed))
    cone5=[(1,0,0,0,0),(0,1,0,0,0),(-1,0,0,1,0),(1,1,1,0,0),(0,0,0,0,1),(1,0,0,0,-1)]
    cone5_certs={}
    for e in g5:
        target=tuple(F(x,2) for x in e);witness=None
        for n in range(1,6):
            for cols in combinations(range(6),n):
                for rr in combinations(range(5),n):
                    coeff=solve_square([[cone5[k][r] for k in cols] for r in rr],[target[r] for r in rr])
                    if coeff is not None and all(c>=0 for c in coeff) and all(sum(coeff[h]*cone5[k][r] for h,k in enumerate(cols))==target[r] for r in range(5)):
                        witness=[F(0)]*6
                        for k,c in zip(cols,coeff):witness[k]=c
                        break
                if witness is not None:break
            if witness is not None:break
        if witness is None or not any(witness):raise AssertionError(('D1 torque cone certificate failed',e))
        cone5_certs[','.join(map(str,e))]=list(map(str,witness))
    record('D1.5.4.torque_expansion_and_cone',14,{'observed_remainder_exponents':len(observed),'cone_certificates':cone5_certs})
    return checks, {}

if __name__ == "__main__":
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_central_body_slice'))
