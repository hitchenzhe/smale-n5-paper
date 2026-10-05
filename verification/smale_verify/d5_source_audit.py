"""Exact coefficients and nonvanishing facts behind the D5 LP rows."""
from itertools import product
from pathlib import Path
import json
from .symbolic import sympy

def verify():
    s=sympy();b,u,theta=s.symbols('b u theta');r=1/(b+1);m=[b*b,b*b/r**2,1/r**2];offsets=[0,r,-r/b]
    centre=sum(mi*si for mi,si in zip(m,offsets))/sum(m)
    centred=[s.factor(si-centre)for si in offsets]
    expected=[-(b-1)*(b+1)/(b*(b*b+b+1)),(2*b+1)/(b*(b+1)*(b*b+b+1)),-b*(b+2)/((b+1)*(b*b+b+1))]
    if any(s.cancel(v-w)for v,w in zip(centred,expected)):raise ArithmeticError('Centred real offsets')
    # Phi=Psi=0: use u=y4/y3 and eliminate y5/y3 by Phi.
    v=-(m[0]+m[1]*u)/m[2]
    F=s.fraction(s.cancel(m[0]+m[1]*u**3+m[2]*v**3))[0]
    H=s.fraction(s.cancel(b/u**2+1/v**2-b-1))[0]
    resultant=s.factor(s.resultant(F,H,u));constant,factors=s.factor_list(resultant);records=[]
    for polynomial,power in factors:
        p=s.Poly(polynomial,b);positive=0 if polynomial==b else p.count_roots(0,s.oo)
        rec=dict(polynomial=str(polynomial),power=int(power),positive_roots=int(positive))
        if positive:
            if p.degree()!=4 or positive!=1:raise ArithmeticError('Unexpected Phi/Psi positive-mass factor')
            root=next(z for z in s.polys.rootoftools.CRootOf.all_roots(p)if z.is_real and z>0)
            field=s.QQ.algebraic_field(root);fp=s.Poly(F.subs(b,root),u,domain=field);hp=s.Poly(H.subs(b,root),u,domain=field)
            common=s.gcd(fp,hp)
            if common.degree()!=1 or not s.rem(s.Poly(u*u-1,u,domain=field),common).is_zero:raise ArithmeticError('Empty regime saturation failed')
            rec.update(common_gcd=str(common.as_expr()),excluded='u^2=1 contradicts 1/y4^2-1/y3^2=-r !=0')
        records.append(rec)
    # The source coefficient identities can be checked before imposing any
    # particular mass or shape: substitute s_k=theta-y_k^-2 exactly.
    yy=s.symbols('y3:6');mm=s.symbols('m3:6');ss=[theta-y**-2 for y in yy]
    Phi=sum(mk*y for mk,y in zip(mm,yy));Psi=sum(mk*y**3 for mk,y in zip(mm,yy));S1=sum(mk*sk*y for mk,sk,y in zip(mm,ss,yy))
    identities={
        'F':sum(mk*y**3*(theta-sk)for mk,y,sk in zip(mm,yy,ss))-Phi,
        'V1':sum(mk*sk*y**3 for mk,y,sk in zip(mm,yy,ss))-(theta*Psi-Phi),
        'V2':sum(mk*sk*(sk-theta)*y**3 for mk,y,sk in zip(mm,yy,ss))+S1,
        'B':sum(mk*y*(1+sk*y**2/2)for mk,y,sk in zip(mm,yy,ss))-(Phi+theta*Psi)/2}
    if any(s.expand(v)for v in identities.values()):raise ArithmeticError('Near-source leading coefficient identity')
    m1,m2,MT=s.symbols('m1 m2 MT');R=s.Matrix([[m1+m2,1,-1],[m2,1+m1/MT,m2/MT],[-m1,m1/MT,1+m2/MT]])
    if (R*s.Matrix([1,-m2,m1])).applyfunc(s.simplify)!=s.zeros(3,1):raise ArithmeticError('R common kernel')
    minors=[]
    from itertools import combinations
    for rows in combinations(range(3),2):
        for cols in combinations(range(3),2):
            value=s.factor(R.extract(rows,cols).det())
            num=s.fraction(value)[0];coeffs=s.Poly(num,m1,m2,MT).coeffs()
            if not coeffs or not(all(c>0 for c in coeffs)or all(c<0 for c in coeffs)):raise ArithmeticError('R coordinate-plane rank')
            minors.append(dict(rows=list(rows),columns=list(cols),determinant=str(value)))
    # Exact positive determinant closes both the inner RRR and the central
    # deep case once their two leading-coefficient relations are derived.
    pair=s.Matrix([[m1,m2],[-m2**2,m1**2]])
    if s.expand(pair.det()-m1**3-m2**3):raise ArithmeticError('Outer leading coefficient contradiction')
    sign_cases=[]
    for signs in product((1,-1),repeat=6):
        parities=[signs[k]*signs[k+3]for k in range(3)]
        tied=[list(p)for p in combinations(range(3),2)if parities[p[0]]==parities[p[1]]]
        if not tied:raise ArithmeticError('Inner cube pigeonhole failed')
        sign_cases.append(dict(signs=list(signs),equal_parity_pair=tied[0]))
    payload=dict(check='D5.near_coefficients_empty_regimes_R_planes_and_deep_exclusions',status='PASS',
        centred_offsets=list(map(str,centred)),near_coefficients=dict(F='Phi',G='Psi',V1='-(theta*Psi-Phi)',V2='-sum m*s*y',B='(Phi+theta*Psi)/2'),
        empty_regime_resultant=str(resultant),empty_regime_factors=records,R_minors=minors,inner_sign_cases=sign_cases,outer_pair_determinant='m1^3+m2^3 >0',
        nonvanishing_proof='The centred s4>0 and s5<0 for every b>0; s3 changes sign only at b=1. Off the slice there is exactly one offset on one side of zero. For y^-2=-s, its term alone supplies one real/imaginary component of Phi, Psi and sum m*s*y, so all three are nonzero. In a dominant-coincidence regime the special offset is nonzero and its edge has strictly lower force order than all other edges. On the central slice the two generic offset contributions have perpendicular phases and cannot cancel for a single source.',
        deep_proof='In the forbidden central deep case A1=A2=g<alpha1,alpha2 and beta1,beta2<-g, tied special forcings imply beta1=beta2=beta. Since q-beta=-3(g+beta)/2>0, the pair-centre W identity forces m1*w1hat+m2*w2hat=0. The special differences share leading coefficient -xi3hat. The cancelled weighted special forcings have equal squares, giving m1^2*w2hat=m2^2*w1hat. The verified positive determinant makes both leading W coefficients zero, a contradiction.',
        inner_proof='At inner RRR, the pair-centre W identity likewise gives opposite weighted W coefficients. Uniform forcing gives one f at all three triangle bodies. If f!=0, squaring its three signed sums yields only two parity-dependent values, forcing two distinct offset cubes to coincide. The exhaustive 64 sign check supplies a repeated parity pair. Thus f=0, equal squared amplitudes give the second pair relation, and the same positive determinant gives a contradiction.',
        scope='Coefficient and sign audit used by the necessary D5 valuation model, including empty coincident-both regimes and the previously implicit central deep exclusion. Does not alone establish full LP coverage.')
    path=Path(__file__).resolve().parents[1]/'certificates/d5_source_audit.json';path.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=64+len(records)+len(minors),detail=dict(certificate=path.name,near_coefficient_identities=4,empty_regimes=True,central_deep_exclusion=True,inner_RRR_exclusion=True))],{}
