"""Exact strict rational feasibility with checkable Motzkin multipliers.

Fourier--Motzkin elimination chooses the least expensive next variable.
Each derived row retains its nonnegative combination of original rows.
No floating LP decision is used. Row-budget exhaustion raises an exception.
"""
from fractions import Fraction as Q
from math import gcd,lcm


def solve(A,b,strict=(),row_limit=50000):
    from verify_lp_certificates import check_case
    n=len(A[0]) if A else 1;strict=set(strict)
    rows=[(tuple(map(Q,a)),Q(rhs),i in strict,{i:Q(1)}) for i,(a,rhs) in enumerate(zip(A,b))]
    original=dict(variables=n,A=[[str(v) for v in a] for a in A],b=[str(v) for v in b],strict=sorted(strict))
    history=[];active=list(range(n))
    def contradiction(row):
        a,rhs,is_strict,y=row
        if any(a) or (rhs>=0 and (rhs!=0 or not is_strict)):return None
        witness=[Q(0)]*len(A)
        for i,value in y.items():witness[i]=value
        case=dict(original,kind='infeasible',y=[str(v) for v in witness])
        check_case(case)
        return case
    def simplify(candidates):
        best={}
        for row in candidates:
            bad=contradiction(row)
            if bad:return None,bad
            a,rhs,flag,y=row
            if not any(a):continue
            denominator=1
            for value in a:denominator=lcm(denominator,value.denominator)
            integers=[int(value*denominator) for value in a];content=0
            for value in integers:content=gcd(content,abs(value))
            scale=Q(denominator,content)
            normalized=tuple(value*scale for value in a)
            row=(normalized,rhs*scale,flag,{i:v*scale for i,v in y.items()})
            old=best.get(normalized)
            if old is None or row[1]<old[1] or (row[1]==old[1] and flag and not old[2]):best[normalized]=row
        return list(best.values()),None
    rows,bad=simplify(rows)
    if bad:return bad
    while active:
        variable=min(active,key=lambda j:sum(a[j]>0 for a,r,s,y in rows)*sum(a[j]<0 for a,r,s,y in rows))
        positive=[r for r in rows if r[0][variable]>0];negative=[r for r in rows if r[0][variable]<0]
        history.append((variable,positive,negative))
        count=len(positive)*len(negative)
        if count>row_limit:raise RuntimeError('Exact LP elimination row budget exceeded; no result claimed')
        candidates=[r for r in rows if not r[0][variable]]
        for pa,pb,ps,py in positive:
            for na,nb,ns,ny in negative:
                pscale=-na[variable];nscale=pa[variable]
                aa=tuple(pscale*a+nscale*c for a,c in zip(pa,na));bb=pscale*pb+nscale*nb
                yy={i:pscale*v for i,v in py.items()}
                for i,v in ny.items():yy[i]=yy.get(i,Q(0))+nscale*v
                candidates.append((aa,bb,ps or ns,yy))
        rows,bad=simplify(candidates)
        if bad:return bad
        active.remove(variable)
    solution=[Q(0)]*n
    for variable,positive,negative in reversed(history):
        def bound(row):
            a,b,flag,y=row
            return (b-sum(v*solution[j] for j,v in enumerate(a) if j!=variable))/a[variable],flag
        upper=[bound(row) for row in positive];lower=[bound(row) for row in negative]
        hi=min((v for v,s in upper),default=None);lo=max((v for v,s in lower),default=None)
        if hi is None and lo is None:value=Q(0)
        elif hi is None:value=lo+1
        elif lo is None:value=hi-1
        elif lo<hi:value=(lo+hi)/2
        elif lo==hi:value=lo
        else:raise ArithmeticError('Exact LP back-substitution found contradictory bounds')
        solution[variable]=value
    case=dict(original,kind='feasible',x=[str(v) for v in solution]);check_case(case)
    return case
