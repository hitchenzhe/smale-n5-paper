"""Exact Laurent arithmetic in q=M**(1/3), over the split cubic field."""
from .d2_fields import E,CF,decode_new,regular_denominators
class L:
 def __init__(self,terms=None):self.terms={k:v for k,v in (terms or {}).items()if v}
 def __bool__(self):return bool(self.terms)
 def __add__(self,o):
  o=convert(o);out=dict(self.terms)
  for k,v in o.terms.items():out[k]=out.get(k,CF.zero)+v
  return L(out)
 __radd__=__add__
 def __neg__(self):return L({k:-v for k,v in self.terms.items()})
 def __sub__(self,o):return self+-convert(o)
 def __rsub__(self,o):return convert(o)+-self
 def __mul__(self,o):
  o=convert(o);out={}
  for i,a in self.terms.items():
   for j,b in o.terms.items():out[i+j]=out.get(i+j,CF.zero)+a*b
  return L(out)
 __rmul__=__mul__
 def inverse(self):
  if len(self.terms)!=1:raise ArithmeticError('Laurent inverse requires one monomial')
  power,value=next(iter(self.terms.items()));return L({-power:value.inverse()})
 def __truediv__(self,o):return self*convert(o).inverse()
 def __pow__(self,n):
  if n<0:return self.inverse()**(-n)
  out=convert(1)
  for _ in range(n):out=out*self
  return out
 def encode(self):return [dict(q_power=k,coefficient=v.encode())for k,v in sorted(self.terms.items())]
 def size(self):return sum(sum(len(x.numer)+len(x.denom)for x in v.v)for v in self.terms.values())
def convert(v):return v if isinstance(v,L)else L({0:CF.convert(v)})
def from_sympy(v):return L({0:CF.from_sympy(v)})
q=L({1:CF.one})
class LF:
 zero=convert(0);one=convert(1)
 convert=staticmethod(convert);from_sympy=staticmethod(from_sympy)

def decode(raw):
 if not isinstance(raw,list):raise ValueError('Laurent coefficient must be a list')
 terms={}
 for term in raw:
  if set(term)!={'q_power','coefficient'}or type(term['q_power'])is not int or term['q_power']in terms:raise ValueError('Malformed or duplicate Laurent monomial')
  terms[term['q_power']]=decode_new(term['coefficient'])
 return L(terms)

def regular(value):
 for component in value.terms.values():regular_denominators(component)
 return True
