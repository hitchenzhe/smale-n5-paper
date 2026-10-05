#!/usr/bin/env python3
"""Exact rational interval verifier for complex polynomial-box certificates.

No floating-point operations are used.  A successful Krawczyk case checks a
polynomial root enclosure using rational rectangular complex intervals,
verifies the supplied rational preconditioner inverse, and verifies a strict
contraction bound.  A cover case verifies every split and each leaf.
"""
from fractions import Fraction
import argparse
import json
from pathlib import Path
import sys
from math import comb


def q(value, where="rational"):
    if isinstance(value, (float, bool)):
        raise ValueError(f"{where}: decimal floats are forbidden; use an integer or rational string")
    try:
        return Fraction(value)
    except (ValueError, TypeError, ZeroDivisionError) as exc:
        raise ValueError(f"{where}: invalid rational {value!r}") from exc


class I:
    def __init__(self, lo, hi=None):
        self.lo=q(lo)
        self.hi=q(lo if hi is None else hi)
        if self.lo>self.hi: raise ValueError("interval endpoints reversed")
    def __add__(self,o):
        o=ii(o); return I(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self): return I(-self.hi,-self.lo)
    def __sub__(self,o): return self+-ii(o)
    def __rsub__(self,o): return ii(o)+-self
    def __mul__(self,o):
        o=ii(o); p=(self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi)
        return I(min(p),max(p))
    __rmul__=__mul__
    def __pow__(self,n):
        if n<0: raise ValueError("negative powers are not supported in polynomial certificates")
        out=I(1)
        for _ in range(n):out=out*self
        return out
    def contains_strictly(self,o): return self.lo<o.lo and o.hi<self.hi
    def abs_upper(self): return max(abs(self.lo),abs(self.hi))
    def __repr__(self): return f"[{self.lo},{self.hi}]"


def ii(x): return x if isinstance(x,I) else I(x)


class C:
    def __init__(self,re,im=0): self.re=ii(re); self.im=ii(im)
    def __add__(self,o):
        o=cc(o); return C(self.re+o.re,self.im+o.im)
    __radd__=__add__
    def __neg__(self): return C(-self.re,-self.im)
    def __sub__(self,o): return self+-cc(o)
    def __rsub__(self,o): return cc(o)+-self
    def __mul__(self,o):
        o=cc(o); return C(self.re*o.re-self.im*o.im,self.re*o.im+self.im*o.re)
    __rmul__=__mul__
    def __pow__(self,n):
        if n<0: raise ValueError("negative powers are not supported in polynomial certificates")
        out=C(1)
        for _ in range(n):out=out*self
        return out
    def abs_upper(self): return self.re.abs_upper()+self.im.abs_upper()
    def strictly_inside(self,box):
        return box.re.contains_strictly(self.re) and box.im.contains_strictly(self.im)
    def __repr__(self): return f"({self.re})+i({self.im})"


def cc(x):
    if isinstance(x,C):return x
    return C(x)


def point_complex(raw, where):
    if not isinstance(raw,list) or len(raw)!=2: raise ValueError(f"{where}: expected [real,imag]")
    return C(q(raw[0],where),q(raw[1],where))


def box_complex(raw, where):
    return C(I(*raw["re"]),I(*raw["im"]))


def parse_polys(raw,n):
    polys=[]
    for pi,terms in enumerate(raw):
        poly=[]
        for ti,term in enumerate(terms):
            if set(term)!={"coefficient","powers"}: raise ValueError(f"polynomial {pi} term {ti}: malformed")
            powers=term["powers"]
            if len(powers)!=n or any(not isinstance(k,int) or k<0 for k in powers):
                raise ValueError(f"polynomial {pi} term {ti}: invalid powers")
            coef=point_complex(term["coefficient"],f"polynomial {pi} coefficient")
            poly.append((coef,powers))
        polys.append(poly)
    if len(polys)!=n: raise ValueError("square system required: one polynomial per variable")
    return polys


def eval_poly(poly,x):
    out=C(0)
    for coef,powers in poly:
        term=coef
        for v,p in zip(x,powers):term=term*(v**p)
        out=out+term
    return out


def eval_derivative(poly,x,j):
    out=C(0)
    for coef,powers in poly:
        p=powers[j]
        if not p:continue
        term=coef*p
        for k,(v,e) in enumerate(zip(x,powers)):
            term=term*(v**(e-1 if k==j else e))
        out=out+term
    return out


def matrix_mul(A,B):
    if not A or not B or len(A[0])!=len(B):raise ValueError("matrix dimension mismatch")
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))),C(0))
             for j in range(len(B[0]))] for i in range(len(A))]


def matrix_vec(A,x):
    return [sum((a*v for a,v in zip(row,x)),C(0)) for row in A]


def check_inverse(A,B):
    n=len(A)
    if len(B)!=n or any(len(row)!=n for row in A+B):raise ValueError("preconditioner and inverse must be square")
    ab=matrix_mul(A,B)
    ba=matrix_mul(B,A)
    for i in range(n):
        for j in range(n):
            expected=C(int(i==j))
            if ab[i][j].re.lo!=expected.re.lo or ab[i][j].im.lo!=0:
                raise ValueError("preconditioner times supplied inverse is not identity")
            if ba[i][j].re.lo!=expected.re.lo or ba[i][j].im.lo!=0:
                raise ValueError("supplied inverse times preconditioner is not identity")


def verify_krawczyk(case,box_override=None):
    n=case["variables"]
    if not isinstance(n,int) or n<1:raise ValueError("variables must be positive")
    polys=parse_polys(case["polynomials"],n)
    box=[box_complex(b,f"box[{i}]") for i,b in enumerate(box_override or case["box"])]
    x0=[point_complex(v,f"center[{i}]") for i,v in enumerate(case["center"])]
    R=[[point_complex(v,f"preconditioner[{i}]") for v in row] for i,row in enumerate(case["preconditioner"])]
    Rinv=[[point_complex(v,f"inverse[{i}]") for v in row] for i,row in enumerate(case["preconditioner_inverse"])]
    if len(box)!=n or len(x0)!=n or len(R)!=n or any(len(r)!=n for r in R):raise ValueError("dimension mismatch")
    for i,(x,b) in enumerate(zip(x0,box)):
        if not b.re.lo<=x.re.lo<=b.re.hi or not b.im.lo<=x.im.lo<=b.im.hi:raise ValueError(f"center[{i}] outside box")
    check_inverse(R,Rinv)

    f0=[eval_poly(p,x0) for p in polys]
    J=[[eval_derivative(polys[i],[box[k] for k in range(n)],j) for j in range(n)] for i in range(n)]
    RJ=matrix_mul(R,J)
    B=[[C(int(i==j))-RJ[i][j] for j in range(n)] for i in range(n)]
    delta=[box[i]-x0[i] for i in range(n)]
    center=[x0[i]-v for i,v in enumerate(matrix_vec(R,f0))]
    correction=matrix_vec(B,delta)
    image=[center[i]+correction[i] for i in range(n)]
    if any(not image[i].strictly_inside(box[i]) for i in range(n)):
        raise ValueError("Krawczyk image is not strictly inside the supplied box")
    contraction=max((sum((entry.abs_upper() for entry in row),Fraction(0)) for row in B),default=Fraction(0))
    if contraction>=1:raise ValueError(f"contraction bound is not below one: {contraction}")
    return {"status":"PASS","kind":"krawczyk","variables":n,
            "krawczyk_image":[{"re":[str(z.re.lo),str(z.re.hi)],"im":[str(z.im.lo),str(z.im.hi)]} for z in image],
            "contraction_bound":str(contraction),"preconditioner_inverse":"exactly verified"}


def verify_exclusion(case,box_override=None):
    n=case["variables"]
    polys=parse_polys(case["polynomials"],n)
    idx=case["polynomial"]
    box=[box_complex(b,f"box[{i}]") for i,b in enumerate(box_override or case["box"])]
    val=eval_poly(polys[idx],box)
    if val.re.lo<=0<=val.re.hi and val.im.lo<=0<=val.im.hi:
        raise ValueError("polynomial enclosure contains zero")
    return {"status":"PASS","kind":"exclude_zero","polynomial":idx,
            "range":{"re":[str(val.re.lo),str(val.re.hi)],"im":[str(val.im.lo),str(val.im.hi)]}}


def verify_uniform_rouche(case):
    """Every parameter-box polynomial root lies in one of n disjoint disks."""
    count=case['parameters']
    params=[box_complex(raw,'parameter box') for raw in case['parameter_box']]
    if type(count) is not int or count<1 or len(params)!=count:raise ValueError('Invalid parameter box dimension')
    terms=[]
    for term in case['polynomial']:
        powers=term['powers']
        if len(powers)!=count+1 or any(type(n) is not int or n<0 for n in powers):raise ValueError('Invalid Rouche monomial')
        terms.append((point_complex(term['coefficient'],'Rouche coefficient'),powers))
    degree=max((powers[-1] for coefficient,powers in terms),default=0)
    disks=case['disks']
    if degree<1 or len(disks)!=degree:raise ValueError('Need exactly degree-many disks')
    def distance_lower(interval):
        return min(abs(interval.lo),abs(interval.hi)) if interval.lo*interval.hi>0 else Fraction(0)
    leading=C(0)
    for coefficient,powers in terms:
        if powers[-1]!=degree:continue
        term=coefficient
        for parameter,power in zip(params,powers):term=term*parameter**power
        leading=leading+term
    if not max(distance_lower(leading.re),distance_lower(leading.im)):
        raise ValueError('Shape degree is not certified constant throughout the mass box')
    results=[];centers=[];radii=[]
    for disk in disks:
        center=point_complex(disk['center'],'disk center');radius=q(disk['radius'])
        if radius<=0:raise ValueError('Rouche disk radius must be positive')
        for old,r in zip(centers,radii):
            if max(abs(center.re.lo-old.re.lo),abs(center.im.lo-old.im.lo))<=radius+r:
                raise ValueError('Rouche disks are not certified disjoint')
        centers.append(center);radii.append(radius)
        coefficients=[C(0) for _ in range(degree+1)]
        for coefficient,powers in terms:
            term=coefficient
            for parameter,power in zip(params,powers):term=term*parameter**power
            for k in range(powers[-1]+1):
                coefficients[k]=coefficients[k]+term*comb(powers[-1],k)*center**(powers[-1]-k)
        dominant=max(distance_lower(coefficients[1].re),distance_lower(coefficients[1].im))*radius
        remainder=coefficients[0].abs_upper()+sum((c.abs_upper()*radius**k for k,c in enumerate(coefficients) if k>=2),Fraction(0))
        if dominant<=remainder:raise ValueError('Rouche strict linear-dominance inequality failed')
        results.append(dict(center=disk['center'],radius=str(radius),dominant_lower=str(dominant),remainder_upper=str(remainder)))
    return dict(status='PASS',kind='uniform_rouche',degree=degree,disks=results,
                scope='one root in each disk for every parameter in the box; disks exhaust all roots')


def split_box(parent,coordinate,cut):
    flat=[]
    for b in parent:flat.extend([[q(v) for v in b["re"]],[q(v) for v in b["im"]]])
    if not isinstance(coordinate,int) or coordinate<0 or coordinate>=2*len(parent):raise ValueError("bad split coordinate")
    cut=q(cut,"split")
    if not flat[coordinate][0]<cut<flat[coordinate][1]:raise ValueError("split must lie strictly inside parent")
    left=[list(v) for v in flat];right=[list(v) for v in flat]
    left[coordinate][1]=cut;right[coordinate][0]=cut
    def regroup(values):
        return [{"re":values[2*i],"im":values[2*i+1]} for i in range(len(parent))]
    return regroup(left),regroup(right)


def same_box(a,b):
    for i,(x,y) in enumerate(zip(a,b)):
        for part in ("re","im"):
            if [q(t) for t in x[part]]!=[q(t) for t in y[part]]:return False
    return len(a)==len(b)


def verify_cover_node(node,parent,base,leaves):
    box=node.get("box")
    if parent is not None and not same_box(box,parent):raise ValueError("child box does not match its parent subdivision")
    if "split" in node:
        if "left" not in node or "right" not in node:raise ValueError("split node needs left and right children")
        left,right=split_box(box,node["split"]["coordinate"],node["split"]["at"])
        verify_cover_node(node["left"],left,base,leaves)
        verify_cover_node(node["right"],right,base,leaves)
        return
    cert=node.get("certificate",{})
    local=dict(base);local["box"]=box
    local.update(cert)
    if cert.get("kind")=="krawczyk": result=verify_krawczyk(local)
    elif cert.get("kind")=="exclude_zero":
        result=verify_exclusion(local)
    else:raise ValueError("leaf certificate kind must be krawczyk or exclude_zero")
    leaves.append(result)


def check(case):
    kind=case.get("kind")
    if kind=="krawczyk":return verify_krawczyk(case)
    if kind=="exclude_zero":return verify_exclusion(case)
    if kind=='uniform_rouche':return verify_uniform_rouche(case)
    if kind=="cover":
        verify_cover_node(case["tree"],case["domain"],case,leaves:=[])
        return {"status":"PASS","kind":"cover","leaf_count":len(leaves),"leaves":leaves}
    raise ValueError("kind must be krawczyk, exclude_zero, or cover")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate",nargs="?",type=Path)
    parser.add_argument("--stdin",action="store_true")
    args=parser.parse_args()
    if bool(args.certificate)==args.stdin:parser.error("provide exactly one path or --stdin")
    obj=json.load(sys.stdin) if args.stdin else json.loads(args.certificate.read_text(encoding="utf-8"))
    cases=obj if isinstance(obj,list) else obj.get("cases",[obj])
    results=[dict({"case":c.get("id","unnamed")},**check(c)) for c in cases]
    print(json.dumps({"checker":"exact rational complex box certificates","file":"<stdin>" if args.stdin else args.certificate.name,
                      "case_count":len(results),"results":results,
                      "warning":"Only supplied box and polynomial claims are certified; exhaustive coverage requires a verified complete split tree and complete manuscript model."},indent=2))

if __name__=="__main__":
    try:main()
    except (KeyError,TypeError,ValueError,json.JSONDecodeError) as exc:
        print(f"FAIL: {exc}",file=sys.stderr);raise SystemExit(1)
