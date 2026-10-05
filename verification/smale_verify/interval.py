"""Outward rounded dyadic rational interval arithmetic (no float operations).

Every operation first uses exact fractions and then enlarges its result to a
2**-120 grid. This controls denominator growth while preserving containment.
"""
from fractions import Fraction as Q
from math import comb,isqrt
GRID=2**120

def rational(v):
    if isinstance(v,(float,bool)):raise ValueError('Use rational strings, never floating inputs')
    return Q(v)

class R:
    def __init__(self,lo,hi=None):
        lo=rational(lo);hi=lo if hi is None else rational(hi)
        if lo>hi:raise ValueError('Reversed interval')
        self.lo=Q((lo*GRID).__floor__(),GRID)
        self.hi=Q((hi*GRID).__ceil__(),GRID)
    def __add__(self,o):
        o=rr(o);return R(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self):return R(-self.hi,-self.lo)
    def __sub__(self,o):return self+-rr(o)
    def __rsub__(self,o):return rr(o)+-self
    def __mul__(self,o):
        o=rr(o);values=[a*b for a in (self.lo,self.hi) for b in (o.lo,o.hi)]
        return R(min(values),max(values))
    __rmul__=__mul__
    def inverse(self):
        if self.lo<=0<=self.hi:raise ArithmeticError('Interval divisor contains zero')
        return R(1/self.hi,1/self.lo)
    def __truediv__(self,o):return self*rr(o).inverse()
    def __rtruediv__(self,o):return rr(o)*self.inverse()
    def square(self):
        return R(0 if self.lo<=0<=self.hi else min(self.lo**2,self.hi**2),max(self.lo**2,self.hi**2))
    def __pow__(self,n):
        if n<0:return self.inverse()**(-n)
        out=R(1)
        for _ in range(n):out=out*self
        return out
    def upper(self):return max(abs(self.lo),abs(self.hi))
    def lower(self):return min(abs(self.lo),abs(self.hi)) if self.lo*self.hi>0 else Q(0)
    def encode(self):return [str(self.lo),str(self.hi)]

def rr(v):return v if isinstance(v,R) else R(v)

class Z:
    def __init__(self,re,im=0):self.re=rr(re);self.im=rr(im)
    def __add__(self,o):
        o=zz(o);return Z(self.re+o.re,self.im+o.im)
    __radd__=__add__
    def __neg__(self):return Z(-self.re,-self.im)
    def __sub__(self,o):return self+-zz(o)
    def __rsub__(self,o):return zz(o)+-self
    def __mul__(self,o):
        o=zz(o);return Z(self.re*o.re-self.im*o.im,self.re*o.im+self.im*o.re)
    __rmul__=__mul__
    def conjugate(self):return Z(self.re,-self.im)
    def inverse(self):
        d=self.re.square()+self.im.square()
        return Z(self.re/d,-self.im/d)
    def __truediv__(self,o):return self*zz(o).inverse()
    def __rtruediv__(self,o):return zz(o)*self.inverse()
    def __pow__(self,n):
        if n<0:return self.inverse()**(-n)
        out=Z(1)
        for _ in range(n):out=out*self
        return out
    def upper(self):return self.re.upper()+self.im.upper()
    def lower(self):return max(self.re.lower(),self.im.lower())
    def modulus_lower(self):
        squared=self.re.lower()**2+self.im.lower()**2
        return Q(isqrt((squared*GRID**2).__floor__()),GRID)
    def midpoint(self):return Z((self.re.lo+self.re.hi)/2,(self.im.lo+self.im.hi)/2)
    def encode(self):return dict(re=self.re.encode(),im=self.im.encode())

def zz(v):return v if isinstance(v,Z) else Z(v)
def decode(raw):return Z(R(*raw['re']),R(*raw['im']))

def certify_power_root(value,n,center,radius):
    """One root of z**n=value in a disk, uniformly for the value box."""
    c=Z(*center);radius=rational(radius)
    if n<2 or radius<=0:raise ValueError('Invalid power-root disk')
    dominant=(n*c**(n-1)).lower()*radius
    error=(c**n-value).upper()+sum((comb(n,k)*(c**(n-k)).upper()*radius**k for k in range(2,n+1)),Q(0))
    if dominant<=error:raise ArithmeticError('Power-root Rouche bound failed')
    # A tiny disk away from zero is separated from its other rotated roots.
    # The n=2/3 applications also supply all disks explicitly when needed.
    return Z(R(c.re.lo-radius,c.re.hi+radius),R(c.im.lo-radius,c.im.hi+radius)),dict(
        degree=n,center=list(map(str,center)),radius=str(radius),dominant_lower=str(dominant),remainder_upper=str(error))

def matrix_product(A,B):
    return [[sum((a*b for a,b in zip(row,col)),Z(0)) for col in zip(*B)] for row in A]

def verified_linear_solve(J,rhs,inverse):
    """Neumann enclosure of J^-1 rhs, with an exact rational preconditioner."""
    n=len(J)
    if len(inverse)!=n or any(len(row)!=n for row in inverse):raise ValueError('Wrong inverse dimension')
    E=matrix_product(inverse,J)
    E=[[Z(int(i==j))-v for j,v in enumerate(row)]for i,row in enumerate(E)]
    contraction=max(sum((v.upper()for v in row),Q(0))for row in E)
    if contraction>=1:raise ArithmeticError('Linear inverse contraction not certified')
    y=[sum((a*b for a,b in zip(row,rhs)),Z(0))for row in inverse]
    correction=contraction*max(v.upper()for v in y)/(1-contraction)
    solution=[Z(v.re+R(-correction,correction),v.im+R(-correction,correction))for v in y]
    return solution,dict(contraction=str(contraction),correction_bound=str(correction))
