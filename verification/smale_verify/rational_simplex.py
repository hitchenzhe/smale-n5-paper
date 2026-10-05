"""Rational two-phase simplex producing independently checked witnesses.

Free variables are split into two nonnegative variables. Strict rows use
a common free slack s<=1, maximized exactly. No tolerance or float decision.
"""
from fractions import Fraction as Q
from verify_lp_certificates import check_case

def tableau(A,b,c,pivot_limit=10000):
    m=len(b);n=len(c);B=[n+i for i in range(m)];N=list(range(n))+[-1]
    D=[list(row)+[Q(-1),rhs]for row,rhs in zip(A,b)]
    D+=[[-q for q in c]+[Q(0),Q(0)],[Q(0)]*n+[Q(1),Q(0)]];pivots=0
    def pivot(r,j):
        nonlocal pivots
        pivots+=1
        if pivots>pivot_limit:raise RuntimeError('Rational simplex pivot budget exhausted')
        p=D[r][j]
        for i in range(m+2):
            if i!=r:
                for k in range(n+2):
                    if k!=j:D[i][k]-=D[r][k]*D[i][j]/p
        for k in range(n+2):
            if k!=j:D[r][k]/=p
        for i in range(m+2):
            if i!=r:D[i][j]/=-p
        D[r][j]=1/p;B[r],N[j]=N[j],B[r]
    def run(phase):
        objective=m+1 if phase==1 else m
        while True:
            eligible=[j for j in range(n+1)if (phase==1 or N[j]!=-1)and D[objective][j]<0]
            if not eligible:return True
            j=min(eligible,key=lambda j:N[j])
            rows=[i for i in range(m)if D[i][j]>0]
            if not rows:return False
            r=min(rows,key=lambda i:(D[i][-1]/D[i][j],B[i]));pivot(r,j)
    if m and min(b)<0:
        r=min(range(m),key=lambda i:D[i][-1]);pivot(r,n)
        if not run(1)or D[m+1][-1]<0:
            y=[Q(0)]*m
            for j,label in enumerate(N):
                if label>=n:y[label-n]=D[m+1][j]
            return None,None,y,pivots
        if D[m+1][-1]!=0:raise ArithmeticError('Nonzero artificial optimum')
        if -1 in B:
            r=B.index(-1);eligible=[j for j in range(n+1)if D[r][j]!=0]
            if eligible:pivot(r,min(eligible,key=lambda j:N[j]))
    if not run(2):raise ArithmeticError('Unexpected unbounded common slack')
    x=[Q(0)]*n;y=[Q(0)]*m
    for i,label in enumerate(B):
        if 0<=label<n:x[label]=D[i][-1]
    for j,label in enumerate(N):
        if label>=n:y[label-n]=D[m][j]
    return D[m][-1],x,y,pivots

def solve(A,b,strict=(),**unused):
    n=len(A[0]);AA=[list(map(Q,row))for row in A];bb=list(map(Q,b));strict=set(strict)
    case=dict(variables=n,A=[[str(q)for q in row]for row in AA],b=list(map(str,bb)),strict=sorted(strict))
    rows=[(row[:],rhs,i in strict,{i:Q(1)})for i,(row,rhs)in enumerate(zip(AA,bb))];history=[]
    def certificate(weights):
        yy=[Q(0)]*len(AA)
        for i,q in weights.items():yy[i]=q
        result=dict(case,kind='infeasible',y=list(map(str,yy)));check_case(result);return result
    def clean(rows):
        best={}
        for aa,rhs,flag,weights in rows:
            nonzero=next((abs(q)for q in aa if q),None)
            if nonzero is None:
                if rhs<0 or(rhs==0 and flag):return None,certificate(weights)
                continue
            aa=[q/nonzero for q in aa];rhs/=nonzero;weights={i:q/nonzero for i,q in weights.items()}
            key=tuple(aa);old=best.get(key)
            if old is None or rhs<old[1]or(rhs==old[1]and flag and not old[2]):best[key]=(aa,rhs,flag,weights)
        return list(best.values()),None
    while True:
        rows,bad=clean(rows)
        if bad:return bad
        lookup={tuple(aa):(rhs,flag,weights)for aa,rhs,flag,weights in rows};chosen=None
        for aa,rhs,flag,weights in rows:
            other=lookup.get(tuple(-q for q in aa))
            if not flag and other is not None and not other[1]and other[0]==-rhs:
                chosen=(aa,rhs,weights,other[2]);break
        if chosen is None:break
        ea,er,ep,en=chosen
        j=max((j for j,q in enumerate(ea)if q),key=lambda j:sum(bool(r[0][j])for r in rows))
        history.append((j,ea[:],er));updated=[]
        for aa,rhs,flag,weights in rows:
            if aa[j]:
                multiplier=aa[j]/ea[j]
                use=en if multiplier>0 else ep
                updated_weights=weights.copy()
                for i,q in use.items():updated_weights[i]=updated_weights.get(i,Q(0))+abs(multiplier)*q
                updated.append(([q-multiplier*p for q,p in zip(aa,ea)],rhs-multiplier*er,flag,updated_weights))
            else:updated.append((aa,rhs,flag,weights))
        rows=updated
    active=[j for j in range(n)if any(aa[j]for aa,rhs,flag,w in rows)];nr=len(active)
    extended=[[aa[j]for j in active]+[Q(flag)]for aa,rhs,flag,w in rows]+[[Q(0)]*nr+[Q(1)]]
    split=[row+[-q for q in row]for row in extended];objective=[Q(0)]*nr+[Q(1)]+[Q(0)]*nr+[Q(-1)]
    optimum,x,y,pivots=tableau(split,[rhs for aa,rhs,flag,w in rows]+[Q(1)],objective)
    if optimum is not None and optimum>0:
        witness=[Q(0)]*n
        for k,j in enumerate(active):witness[j]=x[k]-x[nr+1+k]
        for j,aa,rhs in reversed(history):witness[j]=(rhs-sum(q*witness[k]for k,q in enumerate(aa)if k!=j))/aa[j]
        case.update(kind='feasible',x=list(map(str,witness)))
    else:
        weights={}
        for yy,(_,_,_,w)in zip(y[:-1],rows):
            for i,q in w.items():weights[i]=weights.get(i,Q(0))+yy*q
        return certificate(weights)
    check_case(case) # Always independently check against the ORIGINAL rows.
    return case
