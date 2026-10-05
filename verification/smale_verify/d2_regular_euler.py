"""Global Euler response exclusion through exact spectral resultants.

The only positive masses where the block or response can degenerate are
permutations of (5,1,1) and (2+2*sqrt(2),1,1). At these masses the required
unit-phase cancellation is impossible by strict amplitude dominance.
"""
import json
from fractions import Fraction as Q
from itertools import permutations
from pathlib import Path
from .symbolic import sympy
from .bernstein import corner_sectors,coefficients,positive_cover
from .d2_mass_rays import path_box
from .interval import R,Z,certify_power_root

def verify(data):
    s=sympy();a,b,x=s.symbols('a b x');records=[];root=Path(__file__).resolve().parents[1]
    def encode(poly,gens):return [dict(coefficient=str(c),powers=list(map(int,mon)))for mon,c in s.Poly(poly,*gens).terms()]
    def shape(j):return sum(s.Rational(t['coefficient'])*a**t['powers'].get('m3',0)*b**t['powers'].get('m4',0)*x**t['powers'].get('x',0)for t in data['resultants'][str(j)]['shape_flat'])
    W={3:-(b*x+1)/a,4:x,5:s.Integer(1)};m={3:a,4:b,5:s.Integer(1)};edges=((3,4),(3,5),(4,5));resultants={}
    for j,signs in ((1,(-1,1,1)),(3,(-1,-1,1)),(4,(1,-1,1))):
        eps=dict(zip(edges,signs));cube=-W[3]/sum(m[k]*eps[(3,k)]/(W[k]-W[3])**2 for k in (4,5))
        eigen=s.cancel(cube*sum((m[i]+m[k])*eps[(i,k)]/(W[k]-W[i])**3 for i,k in edges)-1)
        for target in (s.Integer(0),-s.Rational(1,2)):
            num,den=s.fraction(s.cancel(eigen-target));res=s.factor(s.resultant(shape(j),num,x))
            content,factors=s.factor_list(res);resultants[j,target]=(res,factors)
            rec=dict(sheet=j,eigenvalue=str(target),lambda_numerator=encode(num,(a,b,x)),lambda_denominator=encode(den,(a,b,x)),resultant=encode(res,(a,b)),factor_constant=str(content),factors=[dict(polynomial=encode(p,(a,b)),power=int(k))for p,k in factors])
            if j in(1,4):
                valid=s.Poly(shape(j).subs(b,1),x)
                while valid.eval(-1)==0:valid=valid.exquo(s.Poly(x+1,x))
                line=s.factor(s.resultant(valid.as_expr(),num.subs(b,1),x));rec['saturated_b_equals_1_resultant']=encode(line,(a,))
                allowed=a-5 if target==0 and j==1 else a*a-4*a-4 if target!=0 and j==1 else s.Integer(1)
                for p,k in s.factor_list(line)[1]:
                    if s.Poly(p,a).LC()<0:p=-p
                    if s.rem(allowed,p,a)==0:continue
                    if any(c<0 for c in s.Poly(p,a).all_coeffs()):raise ArithmeticError('Unexpected positive line exception')
            records.append(rec)
            print('D2 Euler',j,'spectral target',target,'exact resultant reconstructed',flush=True)
    R8=next(p for p,k in resultants[1,-s.Rational(1,2)][1]if s.total_degree(p)==8)
    Q1=next(p for p,k in resultants[1,0][1]if s.total_degree(p)==2 and p.has(a))
    if s.expand(Q1-((a-s.Rational(5,2)*(b+1))**2+s.Rational(7,4)*(b-1)**2)):raise ArithmeticError('Response mass SOS identity')
    # All other factors are positive mass units, the saturated b=1 line,
    # or the same exceptional polynomial with permuted homogeneous masses.
    R3=next(p for p,k in resultants[3,-s.Rational(1,2)][1]if s.total_degree(p)==8)
    Q3=next(p for p,k in resultants[3,0][1]if s.total_degree(p)==2 and p.has(a))
    if s.expand(R3-a**8*R8.subs({a:1/a,b:b/a},simultaneous=True)).cancel()!=0:raise ArithmeticError('Spectral mass permutation')
    if s.expand(Q3-a**2*Q1.subs({a:1/a,b:b/a},simultaneous=True)).cancel()!=0:raise ArithmeticError('Response mass permutation')
    exceptions={s.expand(Q1),s.expand(Q3),s.expand(R8),s.expand(R3),b-1}
    for res,factors in resultants.values():
        for p,k in factors:
            if s.expand(p)in exceptions:continue
            if p==b*b-b+1:
                if s.expand(p-((b-s.Rational(1,2))**2+s.Rational(3,4))):raise ArithmeticError('Quadratic positivity')
            elif any(c<=0 for c in s.Poly(p,a,b).coeffs()):raise ArithmeticError('Unexpected nonpositive resultant factor')
    terms=[dict(coefficient=str(c),powers=dict(m3=pa,m4=pb,m5=8-pa-pb))for (pa,pb),c in s.Poly(R8,a,b).terms()]
    covers=[];mass_boxes=[]
    for sector in corner_sectors(terms):
        power=sector.pop('power_terms');cover=positive_cover(coefficients(power),52)
        if cover['sign']!=1:raise ArithmeticError('Spectral cover sign')
        if cover['unresolved'] and sector['largest_mass']!='m3':raise ArithmeticError('Unexpected spectral zero sector')
        for hole in cover['unresolved']:
            (rl,rh),(sl,sh)=path_box(hole['path'])
            if sector['substitution']['m5']=='r':mass_boxes.append([(1/rh,1/rl),(sl,sh)])
            else:mass_boxes.append([(1/(rh*sh),1/(rl*sl)),(1/sh,1/sl)])
        covers.append(dict(sector=sector,cover=cover))
    if not mass_boxes:raise ArithmeticError('Missing critical mass neighborhood')
    hull=[R(min(v[i][0]for v in mass_boxes),max(v[i][1]for v in mass_boxes))for i in range(2)]
    def evaluate(poly):return sum((R(str(c))*hull[0]**pa*hull[1]**pb for (pa,pb),c in s.Poly(poly,a,b).terms()),R(0))
    Haa=evaluate(s.diff(R8,a,2));Hbb=evaluate(s.diff(R8,b,2));Hab=evaluate(s.diff(R8,a,b));Hdet=Haa*Hbb-Hab.square()
    if Haa.lo<=0 or Hdet.lo<=0:raise ArithmeticError('Critical neighborhood is not strictly convex')
    minimal=a*a-4*a-4
    for polynomial in (R8,s.diff(R8,a),s.diff(R8,b)):
        if s.rem(polynomial.subs(b,1),minimal,a):raise ArithmeticError('Critical value or gradient not exactly zero')
    sqrt2,root_record=certify_power_root(Z(2),2,('1.41421356237309504880168872421','0'),'1e-26')
    critical=2+2*sqrt2.re
    if not hull[0].lo<critical.lo<critical.hi<hull[0].hi or not hull[1].lo<=1<=hull[1].hi:raise ArithmeticError('Critical point misses convex hull')
    if critical.lo<=4:raise ArithmeticError('Heavy mass not above four')
    valid=s.Poly(shape(1).subs(b,1),x).exquo(s.Poly(x+1,x))
    coeffs=valid.all_coeffs();A,B,C,D,E=coeffs
    if A!=E or B!=D:raise ArithmeticError('Critical shape not reciprocal')
    q=s.cancel(B/(2*A));num=s.fraction(s.cancel(C/A-2-q*q))[0]
    if s.rem(num,minimal,a):raise ArithmeticError('Critical shape is not a squared reciprocal quadratic')
    def at_critical(poly):return sum((R(str(c))*critical**pa for (pa,),c in s.Poly(poly,a).terms()),R(0))
    qnum,qden=s.fraction(q);qbox=at_critical(qnum)/at_critical(qden)
    if qbox.lower()==0 or qbox.upper()>=2:raise ArithmeticError('Critical shape roots not on unit circle')
    f5=s.Poly(valid.as_expr().subs(a,5),x);factor5=s.factor_list(f5.as_expr())[1]
    for polynomial,power in factor5:
        p=s.Poly(polynomial,x)
        if p.degree()!=2 or p.nth(0)!=p.nth(2) or p.nth(1)**2>=4*p.nth(0)*p.nth(2):raise ArithmeticError('Mass-five shape root not unit modulus')
    # a>4 and |x|=1 imply |W3|<=2/a, while |W4|=|W5|=1.
    # Thus a/sqrt(|W3|)>a>4>2. These strict rational inequalities
    # exclude Phi=0 even with arbitrary unit phases on every body.
    payload=dict(check='D2.6.17.Euler_global_spectral_exclusion',status='PASS',spectral_resultants=records,
        R8_positive_cover=covers,critical_hull=[v.encode()for v in hull],hessian=dict(Haa=Haa.encode(),Hbb=Hbb.encode(),Hab=Hab.encode(),determinant=Hdet.encode()),
        critical_mass_minimal_polynomial=encode(minimal,(a,)),sqrt2_root=root_record,critical_a=critical.encode(),shape_quadratic_q=qbox.encode(),
        mass_five_shape_factors=[dict(polynomial=encode(p,(x,)),power=int(k))for p,k in factor5],
        proof=['Euler L has L*1=0 and L*W=W, and is self-adjoint with respect to the mass matrix.',
            'After distance elimination the homogeneous coordinate equations diagonalize to (I+2L)U=0 and (I-L)V=0. Its third eigenvalue is lambda=trace(L)-1.',
            'If rank(I-L)<=1, mass symmetry and the translation kernel force L=I-P_mean, hence D_edge^3=1/MT; this is a fold Lagrange shape, excluded in the regular Euler chart. Thus rank(I-L)=2.',
            'Unless lambda=-1/2, I+2L is invertible; the only homogeneous coordinate direction is rotation. The omitted FZ5 row follows from the torque identity with coefficient m5*W5 !=0, and W5=1 removes rotation. Hence the selected nine-row block is invertible.',
            'Centring, gauge and constant-translation forcing give Psi1=-m4/2*g*X^2. If g=0, eta=0; the two linear equations imply L*zeta=0 and zeta=r. Since r and 1 are independent, rank(L)<=1 and trace(L)=1, hence lambda=0.',
            'Exact source shape resultants exhaust lambda=0 and lambda=-1/2. The response factor is a sum of squares. The spectral factor is positive outside the displayed Bernstein holes and strictly convex in their convex hull, where its value and gradient vanish at a=2+2*sqrt(2), b=1. Thus this is its unique positive zero.',
            'The extraneous b=1 line is checked after removing every x+1 factor (W3=0 is not a valid chart point). The remaining positive exceptions are only permutations of (5,1,1) and (2+2*sqrt(2),1,1).',
            'At each exception all valid x roots have modulus one. The heavy amplitude a/sqrt(abs(W3)) exceeds 4 while the sum of the other two amplitudes is 2. Phi cannot vanish for any unit phases. Consequently neither Jacobian singularity nor response zero occurs at any compatible positive-mass regular Euler label.'],
        depends_on=['D2.6.17.all_regular_X_unit_circle_exclusions','D2.6.2.chart_and_boundary_fibre'])
    path=root/'certificates/d2_Euler_global.json';path.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=6,detail=dict(certificate=path.name,global_mass_domain=True))],{}
