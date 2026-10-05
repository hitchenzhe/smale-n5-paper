"""Reconstruct D2 mass polynomials from Definitions 6.1 and 6.3.

Uses exact injective integer encoding, with a proved coefficient bound,
rather than relying on sampled modular equality.
"""
from pathlib import Path
from math import gcd
from functools import reduce
import json
from .symbolic import sympy
from .packed_resultant import resultant


def reconstruct():
 from ak_resultants import shape_factors
 s=sympy();m3,m4,m5,x,y=s.symbols('m3 m4 m5 x y');w3,w4,w5=s.symbols('w3 w4 w5')
 mapping={'shape_w3':w3,'shape_w4':w4,'shape_w5':w5}
 shapes={j:sum(s.Rational(c.numerator,c.denominator)*s.prod(mapping[n]**e for n,e in mon)
                 for mon,c in p.terms.items()) for j,p in shape_factors()[0].items()}
 f3,f4,f5=m3**2/w3,m4**2/w4,m5**2/w5
 lam=f3**2+f4**2+f5**2-2*f3*f4-2*f3*f5-2*f4*f5
 theta=3*(f3**2+f4**2+f5**2)+2*f3*f5+2*f4*f5-6*f3*f4
 equation=lam*y**2+2*(f3-f4+f5)*(f3-f4-f5)*y+(theta-lam)/2
 def flat(expr):
  polynomial=s.fraction(s.cancel(expr))[0]
  return s.expand(s.fraction(s.cancel(polynomial.subs({w3:-(m4*x+m5)/m3,w4:x,w5:1})))[0])
 def encode(expr,gens):
  return [dict(coefficient=str(c),powers={str(g):e for g,e in zip(gens,mon) if e})
          for mon,c in s.Poly(expr,*gens).terms()]
 def decode(items):
  variables={str(v):v for v in (m3,m4,m5,x,y)}
  return sum(s.Rational(t['coefficient'])*s.prod(variables[n]**e for n,e in t['powers'].items()) for t in items)
 lf,pf=encode(flat(lam),(x,m3,m4,m5)),encode(flat(equation),(x,y,m3,m4,m5))
 result={}
 for j in range(1,9):
  degree=int(s.total_degree(shapes[j]));sf=encode(flat(shapes[j]),(x,m3,m4,m5))
  lr,lb=resultant(sf,lf,degree,divide_m3=False)
  pr,pb=resultant(sf,pf,degree)
  remaining=s.Poly(decode(lr),m3,m4,m5);factors=[]
  for unit in (m3,m4,m5,m3+m4,m3+m5,m4+m5,m3+m4+m5):
   count=0;divisor=s.Poly(unit,m3,m4,m5)
   while True:
    quotient,remainder=remaining.div(divisor)
    if not remainder.is_zero:break
    remaining=quotient;count+=1
   if count:factors.append(dict(power=count,polynomial=encode(unit,(m3,m4,m5))))
  if j==5:
   L5=m3**2+m4**2+m5**2-m3*m4-m3*m5-m4*m5
   remaining=remaining.exquo(s.Poly(L5,m3,m4,m5))
   factors.append(dict(power=1,polynomial=encode(L5,(m3,m4,m5))))
  content,remaining=remaining.primitive()
  factor=encode(remaining.as_expr(),(m3,m4,m5));factors.append(dict(power=1,polynomial=factor))
  from .exact import P,V,equal
  def sparse(items):
   return P({tuple(sorted(t['powers'].items())):__import__('fractions').Fraction(t['coefficient']) for t in items})
  product=P.coerce(int(content))
  for item in factors:product=product*sparse(item['polynomial'])**item['power']
  equal(product,sparse(lr),'D2 L factorization reconstruction')
  if j==5 and (len(factor)!=91 or remaining.total_degree()!=12 or any(s.Rational(t['coefficient'])<=0 for t in factor)):
   raise ArithmeticError('L5 prime positivity or degree/count mismatch')
  even=[t for t in pr if t['powers'].get('y',0)%2==0]
  odd=[t for t in pr if t['powers'].get('y',0)%2]
  result[str(j)]=dict(shape_degree=degree,shape_flat=sf,L_resultant=lr,L_factors=factors,L_constant=str(content),
      L_bound=lb,P_resultant=pr,P_bound=pb,E_terms=even,O_times_y_terms=odd,
      K_terms=[t for t in even if not t['powers'].get('y',0)])
  print(f'D2 S{j}: exact L/P reconstruction and factorization passed ({len(pr)} P terms)',flush=True)
 for j,k in ((1,2),(6,7),(3,3),(4,4),(5,5),(8,8)):
  def swapped(items):
   return {tuple(sorted((('m4' if n=='m3' else 'm3' if n=='m4' else n),e) for n,e in t['powers'].items())):int(t['coefficient']) for t in items}
  actual=swapped(result[str(j)]['P_resultant'])
  expected={tuple(sorted(t['powers'].items())):int(t['coefficient']) for t in result[str(k)]['P_resultant']}
  # The normalization of S1 and S2 differs by a global minus sign.
  if actual!=expected and actual!={m:-c for m,c in expected.items()}:
   raise ArithmeticError(f'P{j}/P{k} body exchange identity failed')
 return dict(lambda_flat=lf,P_flat=pf,resultants=result)


def verify_mass_resultants():
 payload=reconstruct()
 path=Path(__file__).resolve().parents[1]/'certificates'/'d2_mass_resultants.json'
 path.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
 return [dict(check='D2.6.3_6.5.mass_resultants_and_tables',status='PASS',items=8,
  detail=dict(certificate=path.name,method='exact bounded integer encoding of Sylvester determinants',
   L5_prime=dict(degree=12,positive_coefficients=91),
   decompositions='all eight P=E(y^2)+y O(y^2), K=P(0)',
   exchange='P1/P2, P6/P7 and P3/P4/P5/P8 invariance, up to source sign normalization'))],{'d2_mass_resultants':payload}
