"""Q(i,sqrt(3))(nu)[delta]/(delta**3-1/3), with sparse safe inputs.

The base primitive element is theta=i+sqrt(3), theta**4-4*theta**2+16=0.
Certificates contain rational ANP coefficients and integer powers, not code.
"""
from fractions import Fraction as Rational
from .symbolic import sympy
s=sympy();nu=s.symbols('nu');delta=s.root(s.Rational(1,3),3)
base=s.QQ.algebraic_field(s.I,s.sqrt(3));K=base.frac_field(nu)
old=s.QQ.algebraic_field(s.I,s.sqrt(3),delta);oldK=old.frac_field(nu)
if s.expand(old.ext.as_expr()-(s.I+s.sqrt(3)+delta)):
 raise ArithmeticError('Unexpected primitive field embedding')
def norm(v):
 if not v:return K.zero
 lc=v.denom.LC
 return v.field.raw_new(v.numer.quo_ground(lc),v.denom.quo_ground(lc))
class E:
 def __init__(self,values=()):self.v=tuple(norm(q)for q in list(values)+[K.zero]*(3-len(values)))
 def __bool__(self):return any(self.v)
 def __add__(self,o):
  o=convert(o);return E([a+b for a,b in zip(self.v,o.v)])
 __radd__=__add__
 def __neg__(self):return E([-q for q in self.v])
 def __sub__(self,o):return self+-convert(o)
 def __rsub__(self,o):return convert(o)+-self
 def __mul__(self,o):
  o=convert(o);terms=[K.zero]*3
  for i,a in enumerate(self.v):
   if not a:continue
   for j,b in enumerate(o.v):
    if b:terms[(i+j)%3]=norm(terms[(i+j)%3]+norm(a*b/(3 if i+j>=3 else 1)))
  return E(terms)
 __rmul__=__mul__
 def inverse(self):
  a,b,c=self.v
  determinant=norm(a**3+b**3/3+c**3/9-a*b*c)
  if not determinant:raise ZeroDivisionError('Zero cubic-field divisor')
  return E([norm((a*a-b*c/3)/determinant),norm((c*c/3-a*b)/determinant),norm((b*b-a*c)/determinant)])
 def __truediv__(self,o):return self*convert(o).inverse()
 def __pow__(self,n):
  if n<0:return self.inverse()**(-n)
  out=convert(1)
  for _ in range(n):out=out*self
  return out
 def encode(self):
  def poly(p):return [dict(powers=list(m),coefficients=[str(q)for q in c.to_list()])for m,c in p.items()]
  return [dict(numerator=poly(q.numer),denominator=poly(q.denom))for q in self.v]
def convert(v):
 if isinstance(v,E):return v
 return E([K.convert(v)])
deltaE=E([K.zero,K.one]);theta=E([K.from_sympy(s.I+s.sqrt(3))])+deltaE
def oldcoefficient(c):
 out=convert(0)
 for q in c.to_list():out=out*theta+convert(q)
 return out
def from_old(v):
 def poly(p):
  out=convert(0)
  for (power,),coef in p.items():out+=oldcoefficient(coef)*E([K.from_sympy(nu**power)])
  return out
 return poly(v.numer)/poly(v.denom)
def from_sympy(v):return from_old(oldK.from_sympy(v))
def decode_old(raw):
 def poly(terms):
  out=convert(0)
  for term in terms:
   c=old.new([s.QQ.convert(s.Rational(q))for q in term['coefficients']]);out+=oldcoefficient(c)*E([K.from_sympy(nu**term['powers'][0])])
  return out
 return poly(raw['numerator'])/poly(raw['denominator'])
def decode_new(raw):
 if not isinstance(raw,list)or len(raw)!=3:raise ValueError('Need three cubic-basis coefficients')
 def safeq(v):
  if isinstance(v,(float,bool)):raise ValueError('Rational strings required')
  return s.QQ.convert(s.Rational(Rational(v).numerator,Rational(v).denominator))
 def poly(terms):
  result=K.zero
  used=set()
  for term in terms:
   if set(term)!={'powers','coefficients'}or len(term['powers'])!=1:raise ValueError('Malformed coefficient monomial')
   exponent=term['powers'][0]
   if type(exponent)is not int or exponent<0 or exponent in used:raise ValueError('Invalid or duplicate nu power')
   used.add(exponent)
   if len(term['coefficients'])>4:raise ValueError('Base ANP degree must be less than four')
   c=base.new([safeq(q)for q in term['coefficients']])
   result+=K.convert(c,base)*K.from_sympy(nu**term['powers'][0])
  return norm(result)
 for r in raw:
  if set(r)!={'numerator','denominator'}:raise ValueError('Malformed rational coefficient')
 return E([norm(poly(r['numerator'])/poly(r['denominator']))for r in raw])

_regular_factors=None
_regular_checked=set()

def regular_denominators(value):
 """Strip only factors known nonzero for 0<nu<1."""
 global _regular_factors
 if _regular_factors is None:
  U=nu**4-2*s.I*nu**3-6*nu**2+2*s.I*nu+1
  _regular_factors=[K.from_sympy(expr).numer for expr in (nu,nu-1,nu+1,nu**2+1,U,s.conjugate(U).subs(s.conjugate(nu),nu))]
 for component in value.v:
  p=component.denom
  key=tuple(sorted((power,tuple(c.to_list()))for power,c in p.items()))
  if key in _regular_checked:continue
  for factor in _regular_factors:
   while p.degree()>0:
    common=p.gcd(factor)
    if common.degree()==0:break
    p=p.exquo(common)
  if p.degree()>0:raise ArithmeticError('Unproved positive-domain coefficient denominator')
  _regular_checked.add(key)
 return True
class CF:
 zero=convert(0);one=convert(1)
 convert=staticmethod(convert);from_sympy=staticmethod(from_sympy)

