"""Reconstruction of ak_diagram1_separated_center_cone from Table 5, paper page 56.

Only the checks reported by this module are implemented; see program_index.json.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from smale_verify.exact import (P, V, zero, equal, edge, equations, substitute,
    diff_poly, divide_poly, groebner, det, quotient_norm)

def verify_polynomial_syzygies(D,Ew,Ez,m,w,z):
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    pair=(1,2); tri=(3,4,5)
    equal(sum(m[i]*Ew[i] for i in pair),
          sum(m[i]*w[i] for i in pair)-sum(m[i]*m[k]*D[edge(i,k)]**3*(w[i]-w[k]) for i in pair for k in tri),'D1.5.2.w')
    equal(sum(m[k]*Ez[k] for k in tri),
          sum(m[k]*z[k] for k in tri)-sum(m[i]*m[k]*D[edge(i,k)]**3*(z[k]-z[i]) for i in pair for k in tri),'D1.5.2.z')
    TP=sum(m[i]*m[k]*D[edge(i,k)]**3*(z[i]*w[k]-w[i]*z[k]) for i in pair for k in tri)
    equal(TP,-sum(m[i]*(w[i]*Ez[i]-z[i]*Ew[i]) for i in pair),'D1.5.2.torque')
    record('D1.5.2.polynomial_syzygies',3)
    return checks, {'pair':pair, 'tri':tri}

def verify_chart_and_cone(D,m,pair,tri):
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
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
    return checks, {'A':A, 'B':B, 'q':q, 'solve_square':solve_square, 'wh':wh, 'zet':zet, 'zh':zh}

if __name__ == "__main__":
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram1_separated_center_cone'))
