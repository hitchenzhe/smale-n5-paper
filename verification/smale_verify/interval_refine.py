"""Refine an already certified real polynomial root by interval Newton images."""
from .interval import Z,R,rational,matrix_product

def refine_real_root(case,image,iterations=4):
    n=case['variables'];box=[Z(R(*r['re']))for r in image]
    x0=[Z(*r)for r in case['center']]
    pre=[[Z(*v)for v in row]for row in case['preconditioner']]
    def eval_poly(terms,values,column=None):
        out=Z(0)
        for term in terms:
            powers=list(term['powers']);coefficient=Z(*term['coefficient'])
            if column is not None:
                coefficient=coefficient*powers[column];powers[column]-=1
                if not coefficient.upper():continue
            value=coefficient
            for z,p in zip(values,powers):value=value*z**p
            out=out+value
        return out
    f0=[eval_poly(p,x0)for p in case['polynomials']]
    center=[x0[i]-sum((a*b for a,b in zip(row,f0)),Z(0))for i,row in enumerate(pre)]
    for _ in range(iterations):
        J=[[eval_poly(poly,box,j)for j in range(n)]for poly in case['polynomials']]
        E=matrix_product(pre,J);E=[[Z(int(i==j))-v for j,v in enumerate(row)]for i,row in enumerate(E)]
        delta=[v-c for v,c in zip(box,x0)]
        new=[center[i]+sum((a*b for a,b in zip(row,delta)),Z(0))for i,row in enumerate(E)]
        box=[Z(R(max(old.re.lo,v.re.lo),min(old.re.hi,v.re.hi)))for old,v in zip(box,new)]
    return box
