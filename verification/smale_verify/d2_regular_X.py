"""Exact reciprocal-polynomial unit-circle exclusion on regular sheets."""
import json
from pathlib import Path
from .symbolic import sympy

def verify():
    from ak_resultants import shape_factors
    s=sympy();z=s.symbols('z');omega=(-1+s.I*s.sqrt(3))/2;records=[]
    for j in (1,3,4,6,8):
        for phase in (1,2):
            mapping={'shape_w3':omega**phase*z,'shape_w4':z,'shape_w5':s.Integer(1)}
            p=sum(s.Rational(c.numerator,c.denominator)*s.prod(mapping[n]**e for n,e in mon)for mon,c in shape_factors()[0][j].terms.items())
            p=s.Poly(p,z,extension=omega)
            reciprocal=s.Poly(sum(s.conjugate(c)*z**k for k,c in enumerate(p.all_coeffs())),z,extension=omega)
            u,v,g=s.gcdex(p,reciprocal)
            if g.degree()!=0 or g.is_zero or u*p+v*reciprocal!=g:raise ArithmeticError('Reciprocal Bezout exclusion failed')
            def encode(poly):
                return [dict(power=int(m[0]),coefficients=[str(q)for q in poly.rep.dom.from_sympy(c).to_list()])for m,c in poly.terms()]
            records.append(dict(sheet=j,ratio_phase=phase,p=encode(p),reciprocal=encode(reciprocal),bezout_u=encode(u),bezout_v=encode(v),gcd=encode(g)))
    payload=dict(check='D2.6.17.all_regular_X_unit_circle_exclusions',status='PASS',cases=records,
        proof='X=0 implies W3^3=W4^3. Distinct W3,W4 give W3/W4=omega or omega^2. Phi=0 and unit-modulus rho then force abs(W3)=abs(W4)=abs(W5)=1. A unit-circle root would also be a root of the conjugate reciprocal polynomial; the exact Bezout identities exclude every such root. S2/S7 follow by exchanging bodies 3 and 4.',
        coefficient_field='Q(omega), ascending ANP lists as returned by SymPy QQ.algebraic_field(omega)')
    path=Path(__file__).resolve().parents[1]/'certificates/d2_regular_X.json';path.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=10,detail=dict(certificate=path.name))],{}
