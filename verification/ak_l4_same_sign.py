"""Table 5 l4_same_sign: D2 S4 reality and same-sign exclusion."""
from smale_verify.d2_positivity import verify_positive_polynomial,save
from smale_verify.symbolic import sympy


def verify_L4_exclusion(d2_mass_resultants):
 data=d2_mass_resultants['resultants']['4'];cover=verify_positive_polynomial(data['L_factors'][-1]['polynomial'])
 s=sympy();m3,m4,m5,x=s.symbols('m3 m4 m5 x');variables={str(v):v for v in (m3,m4,m5,x)}
 polynomial=sum(s.Integer(t['coefficient'])*s.prod(variables[n]**e for n,e in t['powers'].items()) for t in data['shape_flat'])
 discriminant=s.Poly(s.discriminant(polynomial,x),m3,m4,m5)
 if not discriminant or any(value<=0 for monomial,value in discriminant.terms()):raise ArithmeticError('S4 positive discriminant failed')
 encoded=[dict(powers=list(m),coefficient=str(c)) for m,c in discriminant.terms()]
 certificate=save('d2_L4_same_sign.json',dict(L4_cover=cover,S4_discriminant=encoded,
  explanation='All discriminant coefficients positive; the center-of-mass S4 cubic has three distinct real roots. L4 has no positive-mass zero.'))
 return [dict(check='D2.6.12.L4_exclusion_and_S4_reality',status='PASS',items=7,
  detail=dict(certificate=certificate,discriminant_degree=discriminant.total_degree(),
   discriminant_positive_coefficients=len(discriminant.terms()),L4_sectors=len(cover)))],{}


if __name__=='__main__':
 from smale_verify.runner import program_cli
 raise SystemExit(program_cli('ak_l4_same_sign'))
