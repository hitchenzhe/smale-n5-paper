
from fractions import Fraction as F
from itertools import combinations, permutations, product
import json, sys, hashlib

class P:
    def __init__(self, terms=None):
        self.terms = {m:F(c) for m,c in (terms or {}).items() if c}
    @staticmethod
    def coerce(x):
        return x if isinstance(x,P) else P({():F(x)})
    def __add__(self, other):
        out = self.terms.copy()
        for m,c in self.coerce(other).terms.items():
            out[m] = out.get(m,F(0)) + c
        return P(out)
    __radd__ = __add__
    def __neg__(self): return P({m:-c for m,c in self.terms.items()})
    def __sub__(self, other): return self + -self.coerce(other)
    def __rsub__(self, other): return self.coerce(other) + -self
    def __mul__(self, other):
        out = {}
        for a,x in self.terms.items():
            for b,y in self.coerce(other).terms.items():
                powers = dict(a)
                for name,exp in b: powers[name] = powers.get(name,0)+exp
                m = tuple(sorted((name,exp) for name,exp in powers.items() if exp))
                out[m] = out.get(m,F(0)) + x*y
        return P(out)
    __rmul__ = __mul__
    def __pow__(self,n):
        if not isinstance(n,int): raise TypeError(n)
        if n < 0:
            if len(self.terms)!=1: raise ValueError('inverse requires a monomial')
            ((m,c),) = self.terms.items()
            return P({tuple((v,e*n) for v,e in m):c**n})
        out=P.coerce(1)
        for _ in range(n): out=out*self
        return out
    def __truediv__(self,x):
        return self * (self.coerce(x)**-1)
    def coefficient(self,var,n):
        out={}
        for m,c in self.terms.items():
            if dict(m).get(var,0)==n:
                r=tuple((v,e) for v,e in m if v!=var)
                out[r]=out.get(r,F(0))+c
        return P(out)
    def exponents(self,var): return sorted({dict(m).get(var,0) for m in self.terms})
    def names(self): return {v for m in self.terms for v,e in m}
    def clear(self,var):
        return self * V(var)**max(0,-min(self.exponents(var),default=0))
    def reduce(self,var,coeffs):
        # coeffs are low-to-high coefficients of a monic polynomial.
        out=self
        degree=len(coeffs)-1
        if coeffs[-1]!=1: raise ValueError('monic modulus required')
        if min(out.exponents(var),default=0)<0: raise ValueError('clear denominators first')
        while max(out.exponents(var),default=0)>=degree:
            n=max(out.exponents(var))
            lead=out.coefficient(var,n)
            relation=sum((V(var)**k)*c for k,c in enumerate(coeffs))
            out=out-lead*V(var)**(n-degree)*relation
        return out
def V(name): return P({((name,1),):F(1)})
def zero(p,label):
    p=P.coerce(p)
    if p.terms: raise AssertionError((label,len(p.terms),list(p.terms.items())[:3]))
def equal(a,b,label): zero(a-b,label)
def edge(i,j): return tuple(sorted((i,j)))
def mon_order(p,order):
    powers=dict(p)
    return (sum(powers.values()),tuple(powers.get(v,0) for v in order))
def lead_term(poly,order):
    m=max(poly.terms,key=lambda p:mon_order(p,order))
    return m,poly.terms[m]
def mon_divides(a,b):
    da,db=dict(a),dict(b)
    return all(db.get(v,0)>=e for v,e in da.items())
def mon_quotient(a,b):
    da,db=dict(a),dict(b)
    return tuple(sorted((v,db.get(v,0)-e) for v,e in da.items() if db.get(v,0)!=e)+[(v,e) for v,e in db.items() if v not in da and e])
def mon_lcm(a,b):
    da,db=dict(a),dict(b)
    return tuple(sorted((v,max(da.get(v,0),db.get(v,0))) for v in da.keys()|db.keys() if max(da.get(v,0),db.get(v,0))))
def divide_poly(poly,basis,order):
    f=P(poly.terms); rem=P()
    while f.terms:
        mf,cf=lead_term(f,order); reduced=False
        for g in basis:
            mg,cg=lead_term(g,order)
            if mon_divides(mg,mf):
                qmon=mon_quotient(mg,mf)
                f=f-P({qmon:cf/cg})*g
                reduced=True;break
        if not reduced:
            rem=rem+P({mf:cf});f=f-P({mf:cf})
    return rem
def diff_poly(poly,var):
    out=P()
    for mon,c in poly.terms.items():
        powers=dict(mon);n=powers.get(var,0)
        if n:
            powers[var]-=1
            m=tuple(sorted((v,e) for v,e in powers.items() if e))
            out.terms[m]=out.terms.get(m,F(0))+n*c
    return P(out.terms)
def groebner(polys,order):
    G=[]
    for p in polys:
        r=divide_poly(p,G,order)
        if r.terms:
            _,c=lead_term(r,order);G.append(r*(1/c))
    pairs=[(i,j) for i in range(len(G)) for j in range(i)]
    while pairs:
        i,j=pairs.pop(0)
        mi,ci=lead_term(G[i],order);mj,cj=lead_term(G[j],order)
        l=mon_lcm(mi,mj)
        qi,qj=mon_quotient(mi,l),mon_quotient(mj,l)
        # Buchberger product criterion.
        if all(not (dict(mi).get(v,0) and dict(mj).get(v,0)) for v in set(dict(mi))|set(dict(mj))):
            continue
        s=P({qi:F(1,ci)})*G[i]-P({qj:F(1,cj)})*G[j]
        r=divide_poly(s,G,order)
        if r.terms:
            _,c=lead_term(r,order);r=r*(1/c)
            k=len(G);G.append(r);pairs.extend((k,h) for h in range(k))
            if len(G)>80:raise RuntimeError('Groebner basis exceeded safety limit')
    return G

def equations(m,z,w,D):
    def rows(x):
        return {j:x[j]-sum(m[k]*D[edge(j,k)]**3*(x[j]-x[k]) for k in m if k!=j) for j in m}
    ties={e:D[e]**2*(z[e[1]]-z[e[0]])*(w[e[1]]-w[e[0]])-1 for e in D}
    return rows(z),rows(w),ties

def substitute(poly,mapping):
    out=P()
    for mon,coef in poly.terms.items():
        term=P.coerce(coef)
        for name,power in mon:
            term=term*mapping.get(name,V(name))**power
        out=out+term
    return out

def det(matrix):
    from functools import lru_cache
    n=len(matrix)
    if any(len(row)!=n for row in matrix):raise ValueError('determinant needs a square matrix')
    entries=[[P.coerce(x) for x in row] for row in matrix]
    @lru_cache(None)
    def minor(rows,columns):
        if not rows:return P.coerce(1)
        row=min(rows,key=lambda r:sum(bool(entries[r][c].terms) for c in columns))
        position=rows.index(row)
        remaining=tuple(r for r in rows if r!=row)
        result=P()
        for j,col in enumerate(columns):
            if entries[row][col].terms:
                result=result+(-1)**(position+j)*entries[row][col]*minor(remaining,tuple(c for c in columns if c!=col))
        return result
    return minor(tuple(range(n)),tuple(range(n)))

def solve_rational(matrix,rhs):
    n=len(rhs)
    if len(matrix)!=n or any(len(row)!=n for row in matrix):raise ValueError('square rational system required')
    a=[list(map(F,row))+[F(value)] for row,value in zip(matrix,rhs)]
    for j in range(n):
        pivot=next((i for i in range(j,n) if a[i][j]),None)
        if pivot is None:raise ValueError('singular rational system')
        a[j],a[pivot]=a[pivot],a[j]
        lead=a[j][j];a[j]=[x/lead for x in a[j]]
        for i in range(n):
            if i!=j:
                q=a[i][j];a[i]=[x-q*y for x,y in zip(a[i],a[j])]
    return [row[-1] for row in a]

def normalized_system():
    m={j:V('m'+str(j)) for j in range(1,6)}
    z={j:V('z'+str(j)) for j in m}
    w={j:V('w'+str(j)) for j in m}
    D={e:V('D%d%d'%e) for e in combinations(m,2)}
    Ez,Ew,ties=equations(m,z,w,D)
    return dict(m=m,z=z,w=w,D=D,Ez=Ez,Ew=Ew,ties=ties,pair=(1,2),tri=(3,4,5))

def quotient_norm(poly,variable='k',modulus=(4,0,-2,0,1)):
    degree=len(modulus)-1
    k=V(variable)
    columns=[(poly*k**j).reduce(variable,modulus) for j in range(degree)]
    return det([[columns[j].coefficient(variable,i) for j in range(degree)] for i in range(degree)])
