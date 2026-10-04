
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
m={j:V('m'+str(j)) for j in range(1,6)}
z={j:V('z'+str(j)) for j in m}
w={j:V('w'+str(j)) for j in m}
D={e:V('D%d%d'%e) for e in combinations(m,2)}
def equations(m,z,w,D):
    def rows(x):
        return {j:x[j]-sum(m[k]*D[edge(j,k)]**3*(x[j]-x[k]) for k in m if k!=j) for j in m}
    ties={e:D[e]**2*(z[e[1]]-z[e[0]])*(w[e[1]]-w[e[0]])-1 for e in D}
    return rows(z),rows(w),ties

Ez,Ew,ties=equations(m,z,w,D)
checks=[]
def record(name,count,detail=None):
    checks.append({'check':name,'status':'PASS','items':count,'detail':detail})
equal(sum(m[j]*Ez[j] for j in m),sum(m[j]*z[j] for j in m),'mass center z')
equal(sum(m[j]*Ew[j] for j in m),sum(m[j]*w[j] for j in m),'mass center w')
zero(sum(m[j]*(w[j]*Ez[j]-z[j]*Ew[j]) for j in m),'total torque')
I=sum(m[j]*z[j]*w[j] for j in m)
U=sum(m[i]*m[j]*D[(i,j)] for i,j in D)
vir=sum(m[j]*w[j]*Ez[j] for j in m)+sum(m[i]*m[j]*D[(i,j)]*ties[(i,j)] for i,j in D)
equal(I-U,vir,'I-U ideal certificate')
record('common.mass_center_torque_virial',4)

pair=(1,2); tri=(3,4,5)
equal(sum(m[i]*Ew[i] for i in pair),
      sum(m[i]*w[i] for i in pair)-sum(m[i]*m[k]*D[edge(i,k)]**3*(w[i]-w[k]) for i in pair for k in tri),'D1.5.2.w')
equal(sum(m[k]*Ez[k] for k in tri),
      sum(m[k]*z[k] for k in tri)-sum(m[i]*m[k]*D[edge(i,k)]**3*(z[k]-z[i]) for i in pair for k in tri),'D1.5.2.z')
TP=sum(m[i]*m[k]*D[edge(i,k)]**3*(z[i]*w[k]-w[i]*z[k]) for i in pair for k in tri)
equal(TP,-sum(m[i]*(w[i]*Ez[i]-z[i]*Ew[i]) for i in pair),'D1.5.2.torque')
record('D1.5.2.polynomial_syzygies',3)

tz={k:V('tz'+str(k)) for k in (3,4,5)}
tw={k:V('tw'+str(k)) for k in (3,4,5)}
tD={edge(i,j):V('tD%d%d'%(i,j)) for i,j in combinations((3,4,5),2)}
tphi={k:sum(m[l]*tD[edge(k,l)]**3*(tz[k]-tz[l]) for l in (3,4,5) if k!=l) for k in (3,4,5)}
tpsi={k:sum(m[l]*tD[edge(k,l)]**3*(tw[k]-tw[l]) for l in (3,4,5) if k!=l) for k in (3,4,5)}
tchi={e:tD[e]**2*(tz[e[1]]-tz[e[0]])*(tw[e[1]]-tw[e[0]])-1 for e in tD}
tU=sum(m[i]*m[j]*tD[(i,j)] for i,j in tD)
equal(sum(m[k]*tz[k]*tpsi[k] for k in (3,4,5)),
      tU+sum(m[i]*m[j]*tD[(i,j)]*tchi[(i,j)] for i,j in tD),'D1.5.11 virial identity')
for name in sorted({'tz3','tz4','tz5','tw3','tw4','tw5','tD34','tD35','tD45'}):
    lhs=2*diff_poly(tU,name)+sum(m[k]*(tpsi[k]*diff_poly(tz[k],name)+tphi[k]*diff_poly(tw[k],name)) for k in (3,4,5))
    rhs=sum(m[i]*m[j]*(tD[(i,j)]*diff_poly(tchi[(i,j)],name)-2*tchi[(i,j)]*diff_poly(tD[(i,j)],name)) for i,j in tD)
    equal(lhs,rhs,'D1.5.11 differential coefficient '+name)
record('D1.5.11.triangle_virial_and_differential',10)

# D1 exponents are stored as integers equal to twice the true exponents.
tv=('L','M','R','S')
def T(v): return P({tuple((x,n) for x,n in zip(tv,v) if n):F(1)})
MP=m[1]+m[2]; q=MP*V('MT_inverse')
A={1:m[2]*V('pi'),2:-m[1]*V('pi')}
B={1:m[2]*V('varpi'),2:-m[1]*V('varpi')}
zh,wh=V('Zhat'),V('What')
zz={i:T((0,2,2,0))*zh+T((-2,0,0,0))*A[i] for i in pair}
ww={i:T((0,0,0,2))*wh+T((2,0,0,0))*B[i] for i in pair}
zet={k:V('zeta'+str(k)) for k in tri}; eta={k:V('eta'+str(k)) for k in tri}
for k in tri:
    zz[k]=-q*T((0,2,2,0))*zh+T((0,2,0,0))*zet[k]
    ww[k]=-q*T((0,0,0,2))*wh+eta[k]
DD={e:D[e]*(T((0,0,0,0)) if e==(1,2) else T((1,0,0,0)) if e[0] in pair else T((0,-1,0,0))) for e in D}
ze,we,te=equations(m,zz,ww,DD)
weighted={}
for i in pair:
    weighted['z'+str(i)]=T((2,0,0,0))*ze[i]
    weighted['w'+str(i)]=T((-2,0,0,0))*we[i]
for k in tri:
    weighted['z'+str(k)]=T((0,1,0,0))*ze[k]
    weighted['w'+str(k)]=T((0,3,0,0))*we[k]
for e in D: weighted['tie%d%d'%e]=te[e]
def groups(p):
    out={}
    for mon,c in p.terms.items():
        powers=dict(mon); v=tuple(powers.get(x,0) for x in tv)
        rest=tuple((x,n) for x,n in mon if x not in tv)
        out[v]=out.get(v,P())+P({rest:c})
    return {v:c for v,c in out.items() if c.terms}
g={name:groups(p) for name,p in weighted.items()}
for name,gg in g.items():
    if (0,0,0,0) not in gg: raise AssertionError(('missing c0',name))
    if {'Zhat','What'} & gg[(0,0,0,0)].names(): raise AssertionError(('center in c0',name))
for i in pair:
    equal(g['z'+str(i)][(0,0,0,0)],A[i]*(1-MP*D[(1,2)]**3),'D1 boundary z')
    equal(g['w'+str(i)][(0,0,0,0)],B[i]*(1-MP*D[(1,2)]**3),'D1 boundary w')
for k in tri:
    phi=sum(m[l]*D[edge(k,l)]**3*(zet[k]-zet[l]) for l in tri if l!=k)
    psi=sum(m[l]*D[edge(k,l)]**3*(eta[k]-eta[l]) for l in tri if l!=k)
    equal(g['z'+str(k)][(0,0,0,0)],-phi,'D1 boundary triangle z')
    equal(g['w'+str(k)][(0,0,0,0)],-psi,'D1 boundary triangle w')
for e in D:
    i,j=e
    expect=MP**2*D[e]**2*V('pi')*V('varpi')-1 if e==(1,2) else -D[e]**2*A[i]*eta[j]-1 if i in pair else D[e]**2*(zet[j]-zet[i])*(eta[j]-eta[i])-1
    equal(g['tie%d%d'%e][(0,0,0,0)],expect,'D1 boundary tie')

cone=[(F(1),0,0,0),(0,F(1),0,0),(-1,0,0,F(1)),(0,F(3,2),F(1),0),(F(1),F(1),F(1),0)]
def solve_square(a,b):
    n=len(b); a=[list(map(F,row))+[F(rhs)] for row,rhs in zip(a,b)]
    for j in range(n):
        pivot=next((k for k in range(j,n) if a[k][j]),None)
        if pivot is None: return None
        a[j],a[pivot]=a[pivot],a[j]
        p=a[j][j]; a[j]=[v/p for v in a[j]]
        for k in range(n):
            if k!=j:
                h=a[k][j]; a[k]=[x-h*y for x,y in zip(a[k],a[j])]
    return [row[-1] for row in a]
certs={}
for e in sorted(set(v for gg in g.values() for v in gg if any(v))):
    target=tuple(F(v,2) for v in e); witness=None
    for n in range(1,5):
        for cols in combinations(range(5),n):
            for rr in combinations(range(4),n):
                c=solve_square([[cone[k][r] for k in cols] for r in rr],[target[r] for r in rr])
                if c is not None and all(v>=0 for v in c) and all(sum(c[h]*cone[k][r] for h,k in enumerate(cols))==target[r] for r in range(4)):
                    witness=[F(0)]*5
                    for k,v in zip(cols,c): witness[k]=v
                    break
            if witness is not None: break
        if witness is not None: break
    if witness is None: raise AssertionError(('no cone certificate',target))
    if not any(witness): raise AssertionError('zero certificate')
    # Separate certificate verification.
    if any(v<0 for v in witness): raise AssertionError('negative cone multiplier')
    if tuple(sum(witness[k]*cone[k][r] for k in range(5)) for r in range(4))!=target: raise AssertionError('invalid cone identity')
    certs[','.join(map(str,target))]=list(map(str,witness))
record('D1.5.1.chart_and_cone',20,{'distinct_nonzero_exponents':len(certs),'certificates':certs})

# Exact check of the explicit twisted Diagram-3 family in Lemma 7.12.
tau,pp,qq,yy=map(V,('tau','P3','Q3','Y3'))
alpha,beta,eps,sgn=map(V,('alpha','beta','eps','eta'))
d3_relations=[(1-pp)*(1+qq)-alpha*tau**2,
              (1+pp)*(1-qq)-beta*tau**2,
              4*pp*qq*yy**2+1,alpha**3-1,beta**3-1,eps**2-1,sgn**2-1]
d3_order=('alpha','beta','eps','eta','tau','P3','Q3','Y3')
d3_basis=groebner(d3_relations,d3_order)
d3_z={1:P.coerce(1),2:pp,3:P.coerce(-1),4:-pp,5:P()}
d3_w={1:P.coerce(1),2:-qq,3:P.coerce(-1),4:qq,5:P()}
d3_D={(1,2):alpha,(3,4):alpha,(1,4):beta,(2,3):beta,
       (1,3):eps*tau/2,(1,5):-eps*tau,(3,5):-eps*tau,
       (2,4):tau*yy,(2,5):-2*tau*yy,(4,5):-2*tau*yy}
d3_m={1:F(1,2),2:F(1,2),3:F(1,2),4:F(1,2),5:F(1,8)}
d3_rows=[]
for j in range(1,6):
    d3_rows.append(d3_z[j]-sum(d3_m[l]*d3_D[edge(j,l)]**3*(d3_z[j]-d3_z[l]) for l in range(1,6) if l!=j))
    d3_rows.append(d3_w[j]-sum(d3_m[l]*d3_D[edge(j,l)]**3*(d3_w[j]-d3_w[l]) for l in range(1,6) if l!=j))
for e in combinations(range(1,6),2):
    i,j=e
    d3_rows.append(d3_D[e]**2*(d3_z[j]-d3_z[i])*(d3_w[j]-d3_w[i])-tau**2)
zero(d3_z[3]-d3_z[1]-d3_w[3]+d3_w[1],'D3 twisted gauge')
for idx,row in enumerate(d3_rows):
    rem=divide_poly(row,d3_basis,d3_order)
    zero(rem,'D3.7.12 row '+str(idx+1))
record('D3.7.12.twisted_family',len(d3_rows),{'relations':len(d3_relations),'groebner_basis_size':len(d3_basis),'method':'exact ideal reduction over Q'})

# D6 chart, constructed directly from the original equations.
t=V('t'); Z={k:V('Z'+str(k)) for k in m}; W={k:V('W'+str(k)) for k in m}
az={1:-1,2:-1,3:1,4:1,5:1}; aw={1:1,2:1,3:-1,4:-1,5:1}
cross={(1,3),(1,4),(2,3),(2,4)}
zd={k:t**az[k]*Z[k] for k in m}; wd={k:t**aw[k]*W[k] for k in m}
dd={e:(t*D[e] if e in cross else D[e]) for e in D}
rz,rw,rt=equations(m,zd,wd,dd)
wz={1:1,2:1,3:-1,4:-1,5:1}; wy={1:-1,2:-1,3:1,4:1,5:1}
rz={k:t**wz[k]*rz[k] for k in m}; rw={k:t**wy[k]*rw[k] for k in m}
rows=list(rz.values())+list(rw.values())+list(rt.values())+[Z[2]-Z[1]-1,W[4]-W[3]-1]
degrees=[]
for row in rows:
    if min(row.exponents('t'),default=0)<0 or max(row.exponents('t'),default=0)>5: raise AssertionError('D6 polynomial degree')
    degrees.append(max(row.exponents('t'),default=0))
equal(m[1]*rz[1]+m[2]*rz[2]+m[5]*rz[5]+t**2*(m[3]*rz[3]+m[4]*rz[4]),
      m[1]*Z[1]+m[2]*Z[2]+t**2*(m[3]*Z[3]+m[4]*Z[4]+m[5]*Z[5]),'D6 z weighted center')
equal(m[3]*rw[3]+m[4]*rw[4]+m[5]*rw[5]+t**2*(m[1]*rw[1]+m[2]*rw[2]),
      m[3]*W[3]+m[4]*W[4]+t**2*(m[1]*W[1]+m[2]*W[2]+m[5]*W[5]),'D6 w weighted center')
for i,j in cross:
    equal(rt[(i,j)].coefficient('t',0),-D[(i,j)]**2*Z[i]*W[j]-1,'D6 cross boundary')
wingz={1,2,5}; wingw={3,4,5}
for e in set(D)-cross:
    i,j=e
    target=D[e]**2*((Z[j]-Z[i]) if e==(1,2) else (Z[j]-Z[i]) if e not in {(1,5),(2,5)} else -Z[i])*((W[j]-W[i]) if e!=(3,5) and e!=(4,5) else -W[i])-1
    equal(rt[e].coefficient('t',0),target,'D6 wing boundary')
zwing_names={'Z1','Z2','W1','W2','W5','D12','D15','D25'}|{str('m'+str(i)) for i in m}
wwing_names={'Z3','Z4','Z5','W3','W4','D34','D35','D45'}|{str('m'+str(i)) for i in m}
for row in [rz[1],rz[2],rz[5],rw[1],rw[2],rt[(1,2)],rt[(1,5)],rt[(2,5)],rows[-2]]:
    if not row.coefficient('t',0).names()<=zwing_names: raise AssertionError('z wing separation')
for row in [rw[3],rw[4],rw[5],rz[3],rz[4],rt[(3,4)],rt[(3,5)],rt[(4,5)],rows[-1]]:
    if not row.coefficient('t',0).names()<=wwing_names: raise AssertionError('w wing separation')
record('D6.9.2.chart_polynomial_boundary_center_split',22,{'t_degrees':degrees,'note':'rotation substitution item ii is not included yet'})

Sigma=sum(m[i]*m[j]*D[(i,j)]**3*Z[i]*W[j] for i,j in cross)
lhs=sum(m[i]*(W[i]*rz[i]-Z[i]*rw[i]) for i in (1,2))+m[5]*W[5]*rz[5]
Q=m[5]*(W[5]*Z[5]*(1-sum(m[k]*D[edge(k,5)]**3 for k in (1,2,3,4)))+
        sum(m[k]*D[edge(k,5)]**3*W[k]*Z[5] for k in (1,2))+
        sum(m[k]*D[edge(k,5)]**3*W[5]*Z[k] for k in (3,4)))+t**3*sum(m[i]*m[j]*D[(i,j)]**3*W[i]*Z[j] for i,j in cross)
equal(lhs,-t*Sigma+t**2*Q,'D6 9.4 explicit Q')
record('D6.9.4.torque_explicit_Q',1,{'Q_terms':len(Q.terms),'Q_t_degrees':Q.exponents('t')})

# Norms in Q[A,B,k]/(k^4-2k^2+4), using a 4-by-4 multiplication matrix.
k=V('k'); AA=V('A'); BB=V('B'); modulus=[4,0,-2,0,1]
def det(matrix):
    n=len(matrix); out=P()
    for perm in permutations(range(n)):
        inv=sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
        term=P.coerce((-1)**inv)
        for i,j in enumerate(perm): term=term*matrix[i][j]
        out=out+term
    return out
def norm(poly):
    columns=[(poly*k**j).reduce('k',modulus) for j in range(4)]
    return det([[columns[j].coefficient('k',i) for j in range(4)] for i in range(4)])
pol=AA*k**3-4*AA*k+BB*k**3-BB*k**2-4*BB*k+BB
equal(norm(pol),9*(8*AA**2+16*AA*BB+9*BB**2)**2,'D6 primitive mixed norm')
equal(norm(1-2*k**3),1089,'D6 primitive same norm')
equal(norm(2*BB+k**3),16*(BB**2+2)**2,'D6 primitive opposite denominator norm')
equal(norm(2*BB),16*BB**4,'D6 primitive opposite numerator norm')
record('D6.9.9.table_norms_only',4,{'note':'does not verify jet derivation or Jacobian invertibility'})

# Verify the z-wing parametrization in D6 Lemma 9.6.  Eliminate omega using
# omega=lambda^2/2-1; omega^3=1 then gives the monic lambda polynomial below.
lam=V('lam'); ds=V('scale'); cc=V('c'); bb=V('b')
om=lam**2/2-1
aa=(ds**-3-cc*lam**3)/2
wing_sub={'m1':aa,'m2':aa,'m3':bb,'m4':bb,'m5':cc,
          'Z1':F(-1,2),'Z2':F(1,2),
          'W1':-ds*(aa+2*cc*lam),'W2':ds*(aa+2*cc*om*lam),
          'W5':aa*ds*(4*lam**-2-1),
          'D12':ds,'D15':lam*ds,'D25':om*lam*ds}
lambda_mod=[-16,0,12,0,-6,0,1]
def substitute(poly,mapping):
    out=P()
    for mon,coef in poly.terms.items():
        term=P.coerce(coef)
        for name,power in mon:
            base=mapping.get(name,V(name))
            term=term*base**power
        out=out+term
    return out
z_wing_rows=[rz[1].coefficient('t',0),rz[2].coefficient('t',0),rz[5].coefficient('t',0),
             rw[1].coefficient('t',0),rw[2].coefficient('t',0),
             rt[(1,2)].coefficient('t',0),rt[(1,5)].coefficient('t',0),
             rt[(2,5)].coefficient('t',0),Z[2]-Z[1]-1]
for idx,row in enumerate(z_wing_rows):
    reduced=substitute(row,wing_sub).clear('lam').reduce('lam',lambda_mod)
    zero(reduced,'D6 atlas z-wing row '+str(idx+1))
record('D6.9.6.atlas_z_wing_exact_substitution',len(z_wing_rows),
       {'parameters':'rational functions in scale,c,b modulo lam^6-6lam^4+12lam^2-16',
        'w_wing':'follows by z<->w and exchanging equal-mass pairs'})

print(json.dumps({'python':sys.version.split()[0],'arithmetic':'exact rational sparse Laurent polynomials; no third-party libraries','checks':checks,'scope':'partial algebra verification only; no complete theorem verification'},indent=2))
