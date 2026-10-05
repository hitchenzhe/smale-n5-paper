"""Exact nine-dimensional scale algebra, d^3=1/H1 and e^3=1/H2.

All coefficients are rational functions in positive mass ratios A,B,c.
Only monomial inversions are needed by the selected wing Jacobian.
"""
from fractions import Fraction


class ScaleAlgebra:
    def __init__(self, symbolic, lam, kap):
        self.s=symbolic
        self.A,self.B,self.c=symbolic.symbols('A B c')
        p=(symbolic.sqrt(6)+symbolic.I*symbolic.sqrt(2))/2
        self.domain=symbolic.QQ.algebraic_field(symbolic.I,p)
        self.base=self.domain.frac_field(self.A,self.B,self.c)
        self.lam,self.kap=lam,kap
        self.d3=self.norm(self.base.from_sympy(1/(self.c*(2*self.A+lam**3))))
        self.e3=self.norm(self.base.from_sympy(1/(self.c*(2*self.B+kap**3))))
        self.zero=ScaleElement(self,{})
        self.one=self.convert(1)

    def norm(self,v):
        if not v:return self.base.zero
        lc=v.denom.LC
        return v.field.raw_new(v.numer.quo_ground(lc),v.denom.quo_ground(lc))

    def convert(self,v):
        if isinstance(v,ScaleElement):return v
        if isinstance(v,Fraction):v=self.s.Rational(v.numerator,v.denominator)
        return ScaleElement(self,{(0,0):self.norm(self.base.from_sympy(self.s.sympify(v)))})

    def from_sympy(self,expr):
        s=self.s;d,e=s.symbols('d e')
        expr=s.cancel(expr)
        num,den=s.fraction(expr)
        def poly(p):
            out=self.zero
            for (u,v),coef in s.Poly(p,d,e).terms():
                out=out+ScaleElement(self,{(u%3,v%3):self.norm(self.base.from_sympy(coef)*self.d3**(u//3)*self.e3**(v//3))})
            return out
        return poly(num)/poly(den)

    def to_sympy(self,value):
        d,e=self.s.symbols('d e')
        return sum(self.base.to_sympy(v)*d**u*e**w for (u,w),v in value.terms.items())

    def encode(self,value):
        def polynomial(poly):
            return [dict(powers=list(monomial),coefficients=[str(q) for q in coefficient.to_list()])
                    for monomial,coefficient in sorted(poly.items())]
        return [dict(scale_powers=list(scale),numerator=polynomial(coefficient.numer),
                     denominator=polynomial(coefficient.denom))
                for scale,coefficient in sorted(value.terms.items())]

    def decode(self,data):
        def polynomial(items):
            terms={}
            for item in items:
                powers=tuple(item['powers'])
                if len(powers)!=3 or any(type(n) is not int or n<0 for n in powers):raise ValueError('Invalid mass monomial')
                if powers in terms:raise ValueError('Duplicate mass monomial')
                coeff=self.domain.dtype([self.domain.dom.convert(Fraction(q)) for q in item['coefficients']],
                                       self.domain.mod,self.domain.dom)
                terms[powers]=coeff
            return self.base.field.ring.from_dict(terms)
        terms={}
        for item in data:
            scale=tuple(item['scale_powers'])
            if len(scale)!=2 or any(type(n) is not int or not 0<=n<3 for n in scale):raise ValueError('Invalid scale monomial')
            if scale in terms:raise ValueError('Duplicate scale monomial')
            denominator=polynomial(item['denominator'])
            if not denominator:raise ValueError('Zero certificate denominator')
            terms[scale]=self.norm(self.base.field.raw_new(polynomial(item['numerator']),denominator))
        return ScaleElement(self,terms)


class ScaleElement:
    def __init__(self,field,terms):
        self.field=field
        self.terms={k:field.norm(v) for k,v in terms.items() if v}

    def __bool__(self):return bool(self.terms)
    def __neg__(self):return ScaleElement(self.field,{k:-v for k,v in self.terms.items()})
    def __add__(self,other):
        other=self.field.convert(other);result=dict(self.terms)
        for k,v in other.terms.items():result[k]=self.field.norm(result.get(k,self.field.base.zero)+v)
        return ScaleElement(self.field,result)
    __radd__=__add__
    def __sub__(self,other):return self+-self.field.convert(other)
    def __rsub__(self,other):return self.field.convert(other)+-self
    def __mul__(self,other):
        other=self.field.convert(other);f=self.field;result={}
        for (u,v),a in self.terms.items():
            for (w,z),b in other.terms.items():
                q,r=divmod(u+w,3);h,k=divmod(v+z,3)
                value=f.norm(a*b*f.d3**q*f.e3**h)
                result[(r,k)]=f.norm(result.get((r,k),f.base.zero)+value)
        return ScaleElement(f,result)
    __rmul__=__mul__
    def inverse(self):
        if len(self.terms)!=1:raise ArithmeticError('Scale inversion requires a single scale monomial')
        (u,v),a=next(iter(self.terms.items()));f=self.field
        q,r=divmod(-u,3);h,k=divmod(-v,3)
        return ScaleElement(f,{(r,k):f.norm(f.d3**q*f.e3**h/a)})
    def __truediv__(self,other):return self*self.field.convert(other).inverse()
    def __rtruediv__(self,other):return self.field.convert(other)*self.inverse()
    def __pow__(self,n):
        if n<0:return self.inverse()**(-n)
        result=self.field.one;power=self
        while n:
            if n%2:result=result*power
            power=power*power;n//=2
        return result
    def cost(self):return sum(len(v.numer)+len(v.denom) for v in self.terms.values())
