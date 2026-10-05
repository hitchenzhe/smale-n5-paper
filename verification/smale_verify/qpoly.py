"""Exact rational univariate polynomials, Sturm chains and rational functions."""
from fractions import Fraction as F
from math import gcd as int_gcd, lcm
def trim(p):
    p=list(p)
    while p and p[-1]==0:p.pop()
    return p

def add(a,b):
    c=[F(0)]*max(len(a),len(b))
    for i,x in enumerate(a):c[i]+=x
    for i,x in enumerate(b):c[i]+=x
    return trim(c)

def neg(a):return [-x for x in a]
def sub(a,b):return add(a,neg(b))

def mul(a,b):
    if not a or not b:return []
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim(c)

def scale(a,c):return trim([F(c)*x for x in a])

def power(a,n):
    out=[F(1)]
    for _ in range(n):out=mul(out,a)
    return out

def evalp(p,x):
    out=F(0)
    for c in reversed(p):out=out*x+c
    return out

def deriv(p):return trim([F(i)*p[i] for i in range(1,len(p))])

def primitive(p):
    p=trim(p)
    if not p:return []
    den=1
    for x in p:den=lcm(den,x.denominator)
    ints=[int(x*den) for x in p]
    content=0
    for x in ints:content=int_gcd(content,abs(x))
    return [F(x//content) for x in ints]

def divmodp(a,b):
    a=trim(a);b=trim(b)
    if not b:raise ZeroDivisionError("zero polynomial")
    q=[F(0)]*max(0,len(a)-len(b)+1)
    while a and len(a)>=len(b):
        k=len(a)-len(b);c=a[-1]/b[-1];q[k]=c
        a=sub(a,[F(0)]*k+scale(b,c))
    return trim(q),a

def gcd(a,b):
    a=primitive(a);b=primitive(b)
    while b:a,b=b,primitive(divmodp(a,b)[1])
    return scale(a,F(1)/a[-1]) if a else []

def squarefree(p):
    if len(p)<2:return trim(p)
    q,r=divmodp(p,gcd(p,deriv(p)))
    if r:raise ArithmeticError("non-exact squarefree quotient")
    return q

def sturm(p):
    p=primitive(squarefree(p))
    if len(p)<2:return [p]
    seq=[p,primitive(deriv(p))]
    while seq[-1]:
        r=divmodp(seq[-2],seq[-1])[1]
        if not r:break
        seq.append(primitive(neg(r)))
    return seq

def variations(seq,x):
    signs=[]
    for p in seq:
        y=evalp(p,x)
        if y:signs.append(1 if y>0 else -1)
    return sum(x!=y for x,y in zip(signs,signs[1:]))

def count_roots(seq,a,b):
    if evalp(seq[0],a)==0 or evalp(seq[0],b)==0:
        raise ValueError("Sturm endpoints must not be roots")
    return variations(seq,a)-variations(seq,b)

def split_nonroot(p,a,b):
    for den in range(2,100):
        for num in range(1,den):
            x=a+(b-a)*F(num,den)
            if evalp(p,x):return x
    raise ArithmeticError("no non-root rational split found")

def isolate(p,a,b):
    p=squarefree(p);seq=sturm(p);boxes=[]
    def rec(lo,hi,n):
        if not n:return
        if n==1:boxes.append([lo,hi]);return
        mid=split_nonroot(p,lo,hi);left=count_roots(seq,lo,mid)
        rec(lo,mid,left);rec(mid,hi,n-left)
    rec(a,b,count_roots(seq,a,b))
    return p,seq,boxes

def refine(p,seq,a,b):
    mid=split_nonroot(p,a,b)
    if count_roots(seq,a,mid):return a,mid
    return mid,b

def iadd(a,b):return (a[0]+b[0],a[1]+b[1])
def imul(a,b):
    v=(a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1])
    return min(v),max(v)

def eval_interval(p,box):
    y=(F(0),F(0))
    for c in reversed(p):y=iadd(imul(y,box),(c,c))
    return y

class Rat:
    def __init__(self,n,d=None):
        self.n=trim(n if isinstance(n,list) else [F(n)])
        self.d=trim(d if d is not None else [F(1)])
        if not self.d:raise ZeroDivisionError("zero denominator")
    @staticmethod
    def coerce(x):return x if isinstance(x,Rat) else Rat(x)
    def __add__(self,o):
        o=Rat.coerce(o);return Rat(add(mul(self.n,o.d),mul(o.n,self.d)),mul(self.d,o.d))
    __radd__=__add__
    def __neg__(self):return Rat(neg(self.n),self.d)
    def __sub__(self,o):return self+-Rat.coerce(o)
    def __rsub__(self,o):return Rat.coerce(o)+-self
    def __mul__(self,o):
        o=Rat.coerce(o);return Rat(mul(self.n,o.n),mul(self.d,o.d))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=Rat.coerce(o)
        if not o.n:raise ZeroDivisionError("rational-function division by zero")
        return Rat(mul(self.n,o.d),mul(self.d,o.n))
    def __pow__(self,n):
        return Rat(power(self.d,-n),power(self.n,-n)) if n<0 else Rat(power(self.n,n),power(self.d,n))

