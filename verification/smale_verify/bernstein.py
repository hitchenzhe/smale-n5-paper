"""Exact tensor Bernstein coefficients on the rational unit square."""
from fractions import Fraction as F
from math import comb


def coefficients(power_terms):
    n=max((u for u,v in power_terms),default=0)
    m=max((v for u,v in power_terms),default=0)
    first={}
    for (a,b),value in power_terms.items():
        for i in range(a,n+1):
            first[(i,b)]=first.get((i,b),F(0))+F(comb(i,a),comb(n,a))*value
    result=[[F(0) for _ in range(m+1)] for _ in range(n+1)]
    for (i,b),value in first.items():
        for j in range(b,m+1):result[i][j]+=F(comb(j,b),comb(m,b))*value
    return result


def corner_sectors(terms):
    """Three maximal-mass charts, each with two sectors for its zero corner."""
    for largest in ('m3','m4','m5'):
        others=[m for m in ('m3','m4','m5') if m!=largest]
        for order in (others,others[::-1]):
            power={}
            for term in terms:
                masses=term['powers'];a=masses.get(order[0],0);b=masses.get(order[1],0)
                mon=(a+b,b)
                power[mon]=power.get(mon,F(0))+F(term['coefficient'])
            lowest=min(u for u,v in power)
            power={(u-lowest,v):value for (u,v),value in power.items() if value}
            yield dict(largest_mass=largest,substitution={largest:'1',order[0]:'r',order[1]:'r*s'},
                       removed_r_power=lowest,power_terms=power)


def split_axis(matrix,axis):
    def split_vector(values):
        left=[values[0]];right=[values[-1]]
        while len(values)>1:
            values=[(a+b)/2 for a,b in zip(values,values[1:])]
            left.append(values[0]);right.append(values[-1])
        return left,list(reversed(right))
    if axis==1:
        pairs=[split_vector(row) for row in matrix]
        return [[a for a,b in pairs][i] for i in range(len(matrix))],[[b for a,b in pairs][i] for i in range(len(matrix))]
    transpose=[list(column) for column in zip(*matrix)]
    first,second=split_axis(transpose,1)
    return [list(column) for column in zip(*first)],[list(column) for column in zip(*second)]


def positive_cover(matrix,max_depth=20):
    sign=1 if matrix[-1][-1]>0 else -1
    matrix=[[sign*v for v in row] for row in matrix]
    pending=[('',matrix)];leaves=[];unresolved=[]
    while pending:
        path,b=pending.pop()
        flat=[value for row in b for value in row]
        if min(flat)>=0 and b[-1][-1]>0:
            leaves.append(dict(path=path,minimum=str(min(flat)),upper_corner=str(b[-1][-1])))
        elif len(path)>=max_depth:
            unresolved.append(dict(path=path,minimum=str(min(flat))))
        else:
            axis=len(path)%2;a,c=split_axis(b,axis)
            pending.extend([(path+'0',a),(path+'1',c)])
    return dict(sign=sign,leaves=leaves,unresolved=unresolved)


def blowup_critical_corner(power_terms,mode):
    # r=1-h, s=h*k (mode 0) or r=1-h*k, s=h (mode 1).
    result={}
    for (i,j),value in power_terms.items():
        for a in range(i+1):
            monomial=(a+j,j if mode==0 else a)
            result[monomial]=result.get(monomial,F(0))+value*comb(i,a)*(-1)**a
    result={m:v for m,v in result.items() if v}
    lowest=min(a for a,b in result)
    return {(a-lowest,b):v for (a,b),v in result.items()},lowest
