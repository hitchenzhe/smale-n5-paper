"""Table 5 p5_completion: exact K5 positive-mass ray classification.

All labels are covered by the stronger global regular/fold proof; the
redundant exceptional h^3/h^4 jet orders are not needed by that route.
"""
from fractions import Fraction as F
from smale_verify.symbolic import sympy
from smale_verify.packed_resultant import parameter_resultant
from smale_verify.qpoly import sturm,variations,count_roots,evalp
from smale_verify.d2_positivity import save

def verify_all_K5_labels(d2_mass_rays):
 from smale_verify.d2_K5_labels import verify
 return verify(d2_mass_rays)


def verify_K5_mass_ray(d2_mass_resultants):
 s=sympy();data=d2_mass_resultants;m3,m4,m5,x,y=s.symbols('m3 m4 m5 x y');vs={str(v):v for v in (m3,m4,m5,x,y)}
 def decode(terms):return sum(s.Integer(t['coefficient'])*s.prod(vs[n]**e for n,e in t['powers'].items()) for t in terms)
 shape=decode(data['resultants']['5']['shape_flat']);equation=decode(data['P_flat']).subs(y,0)
 disc=s.discriminant(shape,x)
 if s.Poly(disc+3*m3**2*(m3+m4+m5)**2,m3,m4,m5):raise ArithmeticError('S5 nonreal-root discriminant identity failed')
 remainder=s.rem(equation,shape,x);denominator=(m3**2+m3*m4+m4**2)**3
 A=s.cancel(s.expand(remainder).coeff(x,1)*denominator/(m3*(m3-m4)))
 B=s.cancel(s.expand(remainder).coeff(x,0)*denominator/m3)
 # A,B are polynomials. Exact division certifies the extracted m3-m4 factor.
 A=s.Poly(A,m3,m4,m5).as_expr();B=s.Poly(B,m3,m4,m5).as_expr()
 coefficients,bound=parameter_resultant(s,A.subs(m3,1),B.subs(m3,1),m5,m4)
 R=sum(v*m4**i for i,v in enumerate(coefficients))
 positive=989*m4**8+5138*m4**7+13472*m4**6+22718*m4**5+26822*m4**4+22718*m4**3+13472*m4**2+5138*m4+989
 expected=1024*m4**16*(m4+1)**4*(m4**2-m4+1)**4*(m4**2+m4+1)**21*positive
 if s.Poly(R-expected,m4):raise ArithmeticError('K5 off-diagonal resultant factorization failed')
 # Every factor is positive for m4>0; hence simultaneous remainder zeros
 # require m3=m4. Homogeneity reduces that diagonal to (1,1,rho).
 rho=[F(-12),F(-12),F(-3),F(-6),F(-5),F(-4),F(3),F(2),F(1)]
 polynomial=sum(s.Rational(v.numerator,v.denominator)*m5**i for i,v in enumerate(rho))
 restricted=decode(data['resultants']['5']['K_terms']).subs({m3:1,m4:1})
 if s.Poly(restricted-polynomial**2,m5):raise ArithmeticError('K5 diagonal square identity failed')
 seq=sturm(rho);signs=[1 if row[-1]>0 else -1 for row in seq]
 positive_count=variations(seq,F(0))-sum(a!=b for a,b in zip(signs,signs[1:]))
 lo,hi=F('1.5736409'),F('1.573641')
 if positive_count!=1 or count_roots(seq,lo,hi)!=1:raise ArithmeticError('K5 positive rho root census failed')
 certificate=save('d2_K5_mass_ray.json',dict(S5_discriminant=str(disc),
  off_diagonal_resultant=dict(coefficients=[str(v) for v in coefficients],bound=bound,
      factorization=str(expected),positive_factor_coefficients=[989,5138,13472,22718,26822,22718,13472,5138,989]),
  rho_polynomial=[str(v) for v in rho],positive_roots=positive_count,isolating_interval=[str(lo),str(hi)],
  endpoint_values=[str(evalp(rho,lo)),str(evalp(rho,hi))],ray='positive multiples of (1,1,rho)',
  scope='mass classification only; all-class exceptional jet replay remains separate'))
 return [dict(check='D2.6.16.K5_positive_mass_ray',status='PASS',items=5,
  detail=dict(certificate=certificate,ray='R_+(1,1,rho)',positive_rho_roots=1,
              rho_interval=[str(lo),str(hi)],method='exact off-diagonal resultant, positive factors, diagonal square, rational Sturm isolation'))],{}


if __name__=='__main__':
 from smale_verify.runner import program_cli
 raise SystemExit(program_cli('ak_p5_completion'))
