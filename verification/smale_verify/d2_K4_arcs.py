"""Exact completeness of the four printed positive K4 mass arcs."""
from pathlib import Path
import json
from .symbolic import sympy

def verify(data):
    from ak_resultants import shape_factors
    s=sympy();a,b,c,k=s.symbols('a b c k');mapping={'shape_w3':-a*a,'shape_w4':-b*b,'shape_w5':c*c}
    f=s.expand(sum(s.Rational(q.numerator,q.denominator)*s.prod(mapping[n]**p for n,p in mon)for mon,q in shape_factors()[0][4].terms.items()))
    z=s.symbols('z');cubic=s.Poly(f.subs(c*c,z),z)
    expected=z**3+(a*a+b*b)*z**2+(-a**4-b**4+3*a*a*b*b)*z-(a*a-b*b)**2*(a*a+b*b)
    if s.expand(cubic.as_expr()-expected):raise ArithmeticError('S4 square-coordinate cubic')
    if s.expand(f.subs(c,b)-a*a*(5*b**4-a**4))or s.expand(f.subs(c,a)-b*b*(5*a**4-b**4)):raise ArithmeticError('S4 endpoint signs')
    A=4*k*(1+k**4);B=(1-k*k)*(1+6*k*k+k**4);C=(1+k*k)*(6*k*k-1-k**4)
    if s.expand(f.subs({a:A,b:B,c:C})):raise ArithmeticError('Printed K4 rational shape parametrization')
    positive=k**4+6*k**3+6*k*k+6*k+1
    if s.expand(A-C-(k-1)**2*positive)or s.expand(C-B-2*(5*k**4-1)):raise ArithmeticError('Printed K4 parameter ordering')
    SOS=(k**4-3*k*k+1)**2+15*k**4
    derivative=4*(k*k+1)*SOS/((k*k-1)**2*(k**4+6*k*k+1)**2)
    if s.cancel(s.diff(A/B,k)-derivative):raise ArithmeticError('K4 complete monotone shape-ratio parametrization')
    if s.rem(s.expand(A**4-5*B**4),5*k**4-1,k):raise ArithmeticError('K4 open arc lower endpoint')
    II=[(c**3-b**3)*a,(a**3-c**3)*b,(a**3-b**3)*c]
    III=[(c**3-b**3)*a,(c**3+a**3)*b,(a**3+b**3)*c]
    records=[]
    for name,mass,signs in [('II',II,(1,1,-1)),('III',III,(-1,1,-1))]:
        if s.expand(mass[0]*a*a+mass[1]*b*b-mass[2]*c*c):raise ArithmeticError('K4 positive mass centre')
        if s.cancel(sum(sign*mi/xi for sign,mi,xi in zip(signs,mass,(a,b,c)))):raise ArithmeticError('K4 signed inverse-root cancellation')
        actual=[s.expand(mi.subs({a:A,b:B,c:C}))for mi in mass]
        tangent=[s.expand(s.diff(actual[i],k)*actual[2]-actual[i]*s.diff(actual[2],k))for i in(0,1)]
        gcd=s.gcd(s.Poly(tangent[0],k),s.Poly(tangent[1],k))
        interior=gcd
        while interior.eval(1)==0:interior=interior.exquo(s.Poly(k-1,k))
        roots=interior.count_roots(s.Rational(2,3),1)
        if roots:raise ArithmeticError('K4 arc projective derivative vanishes')
        records.append(dict(name=name,mass=[str(q)for q in mass],printed_parameter_masses=[str(s.factor(q))for q in actual],majority_and_odd_signs=list(signs),derivative_gcd=str(s.factor(gcd.as_expr())),derivative_roots_in_open_domain=0))
    F3,F4,F5=s.symbols('F3 F4 F5');p0=(F3-F4)**2+F5**2+2*F5*(F3+F4)
    if s.expand(p0-((F3+F4+F5)**2-4*F3*F4)):raise ArithmeticError('K4 necessary real signs')
    AA,BB,CC=s.symbols('AA BB CC')
    if s.expand(p0.subs({F3:-AA**2,F4:-BB**2,F5:CC**2})-(CC-AA-BB)*(CC+AA+BB)*(CC-AA+BB)*(CC+AA-BB)):raise ArithmeticError('K4 signed-root factorization')
    # The flat (28) at y=0 has degree four and leading coefficient b^2
    # after normalizing m5=1; a positive K4 zero cannot be a common root
    # at infinity, even when the S4 flat cubic loses its leading term.
    ma,mb,x,y=s.symbols('ma mb x y');terms=data['P_flat']
    p=sum(s.Rational(t['coefficient'])*ma**t['powers'].get('m3',0)*mb**t['powers'].get('m4',0)*x**t['powers'].get('x',0)*y**t['powers'].get('y',0)for t in terms)
    if s.Poly(p.subs(y,0),x).LC()!=mb**2:raise ArithmeticError('K4 no common root at infinity')
    payload=dict(check='D2.6.16.K4_complete_four_arc_classification',status='PASS',shape_sextic=str(f),square_coordinate_cubic=str(expected),parameter_shapes=list(map(str,(A,B,C))),arc_records=records,
        derivative_positive_SOS=str(SOS),domain='5^(-1/4)<k<1; positive homogeneous mass scaling; also exchange m3,m4',
        proof=['The already verified S4 real-root discriminant and nonzero flat P0 leading coefficient make every positive K4 zero a finite real common shape.',
            'P0=(f3+f4+f5)^2-4f3*f4. Since f5>0, an opposite-sign f3,f4 pair makes P0 positive; a positive pair also makes P0 positive. Hence W3,W4<0 and W5>0.',
            'Write W=(-a^2,-b^2,c^2), with a>b after exchanging 3,4. The displayed cubic in c^2 has exactly one positive root by Descartes signs. Evaluating at c=b and c=a shows b<c<a precisely when a/b>5^(1/4).',
            'The factored signed-root equation and positive mass-centre equation have positive solutions only for the two sign choices shown. Their one-dimensional positive kernels are exactly II and III. The remaining sign choices have either all equal signs or one strictly wrong coefficient sign.',
            'The printed (A,B,C) solves S4 identically, obeys B<C<A and has strictly increasing A/B from 5^(1/4) to infinity. By uniqueness of the positive c^2 root it exhausts all positive shapes with b<c<a.',
            'Each positive mass kernel therefore lies on one of the two printed arcs or their body-exchanged images; conversely both displayed mass kernels satisfy the original centre and signed-root equations and hence K4=0. Exact derivative gcds exclude singular projective parametrization points.'],
        depends_on=['D2.6.12.L4_exclusion_and_S4_reality','D2.6.3_6.5.mass_resultants_and_tables'],
        scope='Complete printed positive K4 mass-arc classification. Regularity and qprime!=0 are supplied separately by the stronger global Euler certificate.')
    path=Path(__file__).resolve().parents[1]/'certificates/d2_K4_arcs.json';path.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=4,detail=dict(certificate=path.name,complete_positive_arcs=True))],{}
