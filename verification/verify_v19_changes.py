#!/usr/bin/env python3
"""Independent reconstruction of v19 algebra, plus an honest source/package audit.

This is a supplementary check, not a proof of the analytic coverage statements.
It does not import any of the ak_* formula builders or their certificates.
"""
import argparse
from collections import Counter
import hashlib
from itertools import combinations
import json
from pathlib import Path
import re
from smale_verify.symbolic import sympy

s = sympy()
checks = []


def equal(actual, expected, name):
    if s.cancel(s.expand(actual - expected)) != 0:
        raise ArithmeticError(name + ": identity failed")


def record(name, **detail):
    checks.append(dict(check=name, status="PASS", detail=detail))


def global_products():
    pairs = list(combinations(range(1, 6), 2))
    E = s.Matrix([[int(a == b or set(a).isdisjoint(b)) for b in pairs] for a in pairs])
    equal(E.det(), 128, "Petersen determinant")
    if any(not x.is_Integer for x in 128 * E.inv()):
        raise ArithmeticError("Petersen inverse has unexpected denominators")
    if any(sum(E[:, j]) != 4 for j in range(10)):
        raise ArithmeticError("Total-product multiplicities")
    orders = {p: (2 if set(p) <= {1, 2, 3} else -1 if p in
                       ((1, 5), (2, 4), (3, 4), (3, 5)) else 0) for p in pairs}
    products = {p: sum(E[i, j] * orders[q] for j, q in enumerate(pairs))
                for i, p in enumerate(pairs)}
    equal(sum(orders.values()), 2, "Diagram 9 total-product order")
    equal(products[(1, 5)], -1, "Diagram 9 p15 pole")
    equal(products[(2, 4)], -1, "Diagram 9 p24 pole")
    record("V19.global_products", determinant=128, distance_power=128,
           D9_squared_distance_orders={"".join(map(str, p)): n for p, n in orders.items()},
           D9_Delta_order=2, scope="Monomial and valuation arithmetic; AK distance estimates remain analytic inputs")


def normalized_and_inertia():
    m = dict(zip(range(1, 6), s.symbols("m1:6")))
    z = dict(zip(range(1, 6), s.symbols("z1:6")))
    w = dict(zip(range(1, 6), s.symbols("w1:6")))
    D = {p: s.Symbol("d" + "".join(map(str, p))) for p in combinations(range(1, 6), 2)}
    def dij(i, j):
        return D[tuple(sorted((i, j)))]
    Fz = {i: z[i] - sum(m[j]*dij(i,j)**3*(z[i]-z[j]) for j in m if j != i) for i in m}
    Fw = {i: w[i] - sum(m[j]*dij(i,j)**3*(w[i]-w[j]) for j in m if j != i) for i in m}
    I = sum(m[i]*z[i]*w[i] for i in m)
    U = sum(m[i]*m[j]*dij(i,j) for i,j in D)
    T = {p: D[p]**2*(z[p[1]]-z[p[0]])*(w[p[1]]-w[p[0]])-1 for p in D}
    equal(I-U, sum(m[i]*w[i]*Fz[i] for i in m)+sum(m[i]*m[j]*D[i,j]*T[i,j] for i,j in D), "Virial polynomial identity")
    S5 = sum(m[j]*m[l]*dij(j,l)**3*((z[j]-z[5])*(w[l]-w[5])-(w[j]-w[5])*(z[l]-z[5])) for j in (3,4) for l in (1,2))
    rhs = sum(m[k]*z[k]*(w[k]-w[5]) for k in (1,2))
    rhs += sum(m[j]*w[j]*(z[j]-z[5]) for j in (3,4))+S5
    rhs -= sum(m[j]*((z[j]-z[5])*Fw[j]-(w[j]-w[5])*Fz[j]) for j in (3,4))
    rhs += w[5]*sum(m[k]*Fz[k] for k in m)
    equal(I, rhs, "Diagram 7 inertia identity")
    a,b,t,X1,X2,Y1,Y2 = s.symbols("a b t X1 X2 Y1 Y2")
    M=a+b+t
    zc=-(a*X1+b*X2)/M
    wc=-(a*Y1+b*Y2)/M
    outer=a*(X1+zc)*(Y1+wc)+b*(X2+zc)*(Y2+wc)+t*zc*wc
    expected=(a*b*(X1-X2)*(Y1-Y2)+a*t*X1*Y1+b*t*X2*Y2)/M
    equal(outer,expected,"Diagram 5 inertia decomposition")
    C=s.Matrix([[1+a/t,b/t],[a/t,1+b/t]])
    Cinv=s.Matrix([[t+b,-b],[-a,t+a]])/M
    if s.simplify(C*Cinv-s.eye(2))!=s.zeros(2):
        raise ArithmeticError("Diagram 5 force matrix inverse")
    equal(C.det(), M/t, "Diagram 5 force matrix determinant")
    record("V19.polynomial_inertia_and_outer_equations", identities=5,
           scope="Original coordinate equations reconstructed from scratch, without a gauge or certificate assumptions")


def triangle_hessian():
    A,B,rho,ell=s.symbols("A B rho ell", nonzero=True)
    K=ell*s.Matrix([[3/A**2,1/(A*B)],[1/(A*B),3/B**2]])/(4*rho)
    equal(K.det(),ell**2/(2*rho**2*A**2*B**2),"D5 transverse Hessian determinant")
    equal(K.inv()[0,0],3*A**2*rho/(2*ell),"D5 inverse Hessian 11")
    equal(K.inv()[1,1],3*B**2*rho/(2*ell),"D5 inverse Hessian 22")
    # Use rho^2=A B only after computing the inverse/determinant.
    equal((K.det()-ell**2/(2*rho**6)).subs(B,rho**2/A),0,"D5 rho determinant")
    u,v,h1,h2,b1,b2,eta=s.symbols("u v h1 h2 b1 b2 eta")
    equal(9+3*u+6*v-3*(u+(1+3*v)/2),(15+3*v)/2,"D5 mixed lower bound")
    equal((3+b2-b1)/3+b1,(3+b2+2*b1)/3,"D5 single contact balance")
    substitutions={h1:h2-eta,b1:-h2+3*eta,b2:-h2}
    equal((-(h1+3*b1)/2).subs(substitutions),h2-4*eta,"D5 double contact F1")
    equal((-(h2+3*b2)/2).subs(substitutions),h2,"D5 double contact F2")
    equal((h1-s.Rational(3,2)*(h1+b2)).subs(substitutions),h2+eta/2,"D5 double contact pair force")
    record("V19.D5_hessian_and_contact_arithmetic", identities=9,
           scope="Exact transverse Hessian and order identities; finite far/near case coverage reviewed in the audit report")


def d8_upper():
    X,q,r,Y=s.symbols("X q r Y", nonzero=True)
    h,mu,e=s.symbols("h mu e")
    a,b,c,d,f,g=s.symbols("a b c d f g")
    x=a*h+b*mu+c*h*h
    v=d*h+f*mu+g*h*h
    closure= q**-2-(Y*r*r)**-1+X**-1
    E2=-q**3-Y*r**3/4+mu*Y
    variables={X:-s.Rational(1,2)+x,Y:s.Rational(1,2)+x,q:-s.Rational(1,2)+v,r:1+h}
    def jet(expr,substitution):
        value=expr.subs(substitution, simultaneous=True).subs({h:e*h,mu:e*e*mu}, simultaneous=True)
        return s.expand(s.series(value,e,0,3).removeO().subs(e,1))
    equations=[]
    for expr in (closure,E2):
        J=jet(expr,variables)
        equations.extend([s.expand(J).coeff(h,1).coeff(mu,0),s.expand(J).coeff(mu,1),s.expand(J).coeff(h,2)])
    solutions=s.solve(equations,(a,b,c,d,f,g),dict=True)
    expected={a:-s.Rational(3,4),b:2,c:0,d:-s.Rational(1,4),f:0,g:s.Rational(3,8)}
    if solutions!=[expected]:
        raise ArithmeticError("D8 independent upper implicit expansion: "+str(solutions))
    substitutions={key:s.expand(val.subs(expected)) for key,val in variables.items()}
    E0=q**3-X/4+mu*X
    beta=2*q-(X**-1-r/Y)/4
    equal(jet(E0,substitutions),-mu+3*h*h/16,"D8 first balance")
    equal(jet(beta,substitutions),15*h*h/4,"D8 beta second-order coefficient")
    equal(s.Matrix([[-1,s.Rational(3,16)],[-4,s.Rational(15,4)]]).det(),-3,"D8 final pairing determinant")
    alpha=s.symbols("alpha",positive=True)
    N=alpha*alpha+1
    m5=alpha**3/N**2
    boundary={X:-1/N,Y:alpha**2/N,q:-alpha/N,r:1}
    dq=alpha/3
    dr=-N*N/(3*alpha*alpha)
    tangent={X:1,Y:1,q:dq,r:dr}
    for expr in (q**3-m5*X,-alpha**2*q**3-m5*Y*r**3):
        equal(sum(s.diff(expr,var).subs(boundary)*component for var,component in tangent.items()),0,"D8 upper force tangent")
    equal(sum(s.diff(closure,var).subs(boundary)*component for var,component in tangent.items()),N*N*(1-alpha**4)/(3*alpha**4),"D8 reduced closure derivative")
    record("V19.D8_direct_upper_expansion", implicit_coefficients={str(k):str(v) for k,v in expected.items()},
           first_balance="-mu+3h^2/16", beta="15h^2/4", final_determinant=-3,
           scope="Independent scalar elimination from the printed exact force and distance equations")


def printed_resultants():
    k,A,B=s.symbols("k A B")
    poly=k**4-2*k*k+4
    mixed=A*k**3-4*A*k+B*k**3-B*k*k-4*B*k+B
    equal(s.resultant(poly,mixed,k),9*(8*A*A+16*A*B+9*B*B)**2,"D6 mixed norm")
    equal(s.resultant(poly,1-2*k**3,k),1089,"D6 same root norm")
    equal((2*B)**4/s.resultant(poly,2*B+k**3,k),B**4/(B*B+2)**2,"D6 opposite root norm")
    equal(B**4-8*B*B+64,(B*B-4)**2+48,"D6 determinant positivity")
    n=s.symbols("n",real=True)
    equal(n**4-2*s.I*n**3-6*n*n+2*s.I*n+1,(3*(n-s.I)**4+(n+s.I)**4)/4,"D2 Upsilon")
    t=s.symbols("t",positive=True)
    equal(1+t**15-t**4-t**13,(1-t**4)*(1-t**13)+t**15-t**17,"D3 minus inequality identity")
    equal((1-t**3+t**6-t**9+t**12)*(1+t**9),(1+t**15)*(1-t**3+t**6),"D3 minus polynomial ratio")
    record("V19.printed_norms_and_nonvanishing", identities=7,
           scope="Printed D6 norm factors and elementary D2/D3 nonvanishing identities")


def source_audit(root):
    path=root/"Smale_5-CT-v19.tex"
    text=path.read_text(encoding="utf-8-sig")
    clean=re.sub(r"(?<!\\)%[^\n]*","",text)
    labels=re.findall(r"\\label\{([^}]+)\}",clean)
    references=re.findall(r"\\(?:eqref|ref|pageref|autoref)\{([^}]+)\}",clean)
    missing=sorted(set(references)-set(labels))
    duplicates={k:v for k,v in Counter(labels).items() if v>1}
    if missing or duplicates:
        raise ArithmeticError("Source labels: "+str((missing,duplicates)))
    here=Path(__file__).resolve().parent
    index=json.loads((here/"program_index.json").read_text(encoding="utf-8"))
    missing_program_labels={p["program"]:sorted(set(p.get("source_labels",[]))-set(labels)) for p in index["programs"]}
    missing_program_labels={k:v for k,v in missing_program_labels.items() if v}
    appendix=clean.split(r"\section{Computational statements used in the proof}",1)[1].split(r"\section*{Declaration of AI use}",1)[0]
    identifiers=[x.replace(r"\_","_") for x in re.findall(r"\\code\{([^}]+)\}",appendix)]
    functions=sorted({x for x in identifiers if x.startswith("verify_")})
    function_modules={name:sorted({p["program"] for p in index["programs"]
        if p.get("module") and (here/p["module"]).is_file() and name in p.get("public_functions",[])})
        for name in functions}
    unknown_functions=[name for name,modules in function_modules.items() if not modules]
    if unknown_functions:
        raise ArithmeticError("Appendix function not implemented: "+str(unknown_functions))
    modules=[x if x.startswith("ak_") else "ak_"+x for x in identifiers
             if x!="ak_" and not x.startswith("verify_")]
    unavailable=[x for x in modules if not(here/(x+".py")).exists()]
    replacements={p["program"]:p.get("scope","") for p in index["programs"] if p.get("status")=="COVERED_BY_REPLACEMENT"}
    return dict(manuscript=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest().upper(),
                labels=len(labels),missing_references=missing,duplicate_labels=duplicates,
                outdated_active_code_labels=missing_program_labels,
                appendix_verification_functions=function_modules,
                appendix_programs_without_same_named_file=unavailable,
                documented_replacements={x:replacements.get(x,"UNMAPPED") for x in unavailable},
                legacy_version_mentions=[dict(line=i+1,text=line) for i,line in enumerate(text.splitlines()) if "v10" in line or "7389c875" in line],
                scope="Program aliases and documentary claims are reported, not silently accepted as exact formula replays")


def verify_reconstructed_identities():
    checks.clear()
    for fn in (global_products,normalized_and_inertia,triangle_hessian,d8_upper,printed_resultants):
        fn()
    return list(checks), {}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root",type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    verify_reconstructed_identities()
    for check in checks:
        print(check["check"]+": PASS",flush=True)
    report=dict(status="PASS",checks=checks,source_audit=source_audit(args.source_root),
                scope="Supplementary exact reconstructions; not a formal proof of the manuscript or its analytic coverage")
    destination=args.output or args.source_root/"verification"/"v19-records"/"v19_changes.json"
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print("Record:",destination)


if __name__=="__main__":
    main()
