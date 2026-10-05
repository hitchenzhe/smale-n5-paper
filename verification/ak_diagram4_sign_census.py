"""Table 5: ak_diagram4_sign_census, all 64 signs by exact Sturm chains."""
import argparse, json, sys
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from smale_verify.qpoly import *
from ak_diagram4_first_order_formula import vectors
from ak_diagram4_zero_uniqueness import positive_denominator_numerator, sign_at_root
from ak_diagram4_third_order_formula import second_obstruction
from ak_diagram4_sign_orbits import positive_variants
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary',action='store_true',help='print only the verified counts')
    parser.add_argument('--output',type=Path,help='write the full rational root-count record to this JSON file')
    args=parser.parse_args()
    n=[F(0),F(1)];n2=power(n,2);n4=power(n,4)
    b=scale(add(n4,[F(1)]),2)
    c=add(add(neg(n4),scale(n2,6)),[-F(1)])
    bc=add(b,c);n2p1=add(n2,[F(1)])
    m=[Rat(power(mul(b,c),2)),Rat(power(mul(b,bc),2)),Rat(power(mul(c,bc),2))]
    g4=Rat(mul(bc,sub(n2,[F(1)])),mul(c,n2p1))
    g5=Rat(mul(bc,sub(n2,[F(1)])),scale(mul(n,b),2))
    left=F(1);lo=F(2414213,1000000);hi=F(2414214,1000000)
    # Verify 1+sqrt(2) is inside the bracket by squaring positive endpoints.
    if not ((lo-1)**2<2 and (hi-1)**2>2):raise AssertionError("bad sqrt(2) bracket")
    theorem_set={s for s in product((1,-1),repeat=6) if s[0]*s[4]==1 and s[1]*s[3]==-1}
    positive_counts={s:0 for s in product((1,-1),repeat=6)}
    records=[];second_inputs=[]
    # Normalize the two triples by their first signs.  Flipping either entire
    # triple preserves the root set of H; it only multiplies H and R by the
    # product of the two flip signs.  This reduces 64 Sturm runs to 16 while
    # reconstructing all 64 sign cases exactly.
    for s14,s15,s24,s25 in product((1,-1),repeat=4):
        signs=(1,s14,s15,1,s24,s25)
        X1,Y1=vectors(signs,m,g4,g5,1)
        H=X1[0]*Y1[0]+X1[1]*Y1[1]
        R=X1[0]*Y1[1]-X1[1]*Y1[0]
        h=positive_denominator_numerator(H);r=positive_denominator_numerator(R)
        if not h:raise AssertionError((signs,"identically-zero eliminant"))
        if evalp(h,lo)==0 or evalp(h,hi)==0:raise AssertionError((signs,"rational bracket endpoint root"))
        hs,seq,boxes=isolate(h,left,lo)
        if count_roots(sturm(h),lo,hi):raise AssertionError((signs,"root in endpoint bracket"))
        rs=[sign_at_root(r,hs,seq,*box) for box in boxes]
        root_records=[]
        for box,sg in zip(boxes,rs):
            variants=positive_variants(signs,sg)
            for full in variants:positive_counts[full]+=1
            root_records.append({"interval":[str(a) for a in box],"R_sign":sg,
                                 "positive_sign_vectors":[list(v) for v in variants]})
            if variants:second_inputs.append((signs,box))
        records.append({"normalized_signs":list(signs),"H_numerator_coefficients_ascending":[str(x) for x in h],
                        "R_numerator_coefficients_ascending":[str(x) for x in r],
                        "eliminant_root_count":len(boxes),
                        "roots":root_records})

    expected={s for s,nroots in positive_counts.items() if nroots}
    if any(nroots>1 for nroots in positive_counts.values()) or expected!=theorem_set:
        raise AssertionError({"found":sorted(expected),"expected":sorted(theorem_set),
                              "root_multiplicities":{str(k):v for k,v in positive_counts.items() if v>1}})

    second=[]
    for signs,root_box in second_inputs:
        X1,Y1=vectors(signs,m,g4,g5,1);X3,Y3=vectors(signs,m,g4,g5,3)
        H=X1[0]*Y1[0]+X1[1]*Y1[1];R=X1[0]*Y1[1]-X1[1]*Y1[0]
        h=positive_denominator_numerator(H);hs,seq,boxes=isolate(h,left,lo)
        obstruction=second_obstruction(signs,m,g4,g5)
        real=gcd(hs,positive_denominator_numerator(obstruction[0]))
        imag=gcd(hs,positive_denominator_numerator(obstruction[1]))
        box=root_box
        rz=len(real)>1 and count_roots(sturm(real),*box)>0
        iz=len(imag)>1 and count_roots(sturm(imag),*box)>0
        if rz and iz:raise AssertionError((signs,"second obstruction may vanish at root"))
        second.append({"signs":list(signs),"real_component_may_vanish":rz,
                       "imag_component_may_vanish":iz,"root_interval":[str(x) for x in box],
                       "real_gcd_coefficients_ascending":[str(x) for x in real],
                       "imag_gcd_coefficients_ascending":[str(x) for x in imag]})

    result={"check":"D4 Lemma 8.10 exact sign census","arithmetic":"rational polynomials and Sturm sequences; no floats",
            "domain":"1<n<1+sqrt(2); rational endpoint bracket checked",
            "normalized_sturm_cases":len(records),"vectors":64,"with_unique_positive_pair":len(expected),
            "without_positive_pair":64-len(expected),"second_obstruction_checks":len(second),
            "records":records,"second_obstruction":second}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+"\n",encoding='utf-8')
    if args.summary:
        result={k:v for k,v in result.items() if k not in {'records','second_obstruction'}}
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    try:main()
    except (ArithmeticError,AssertionError,ValueError,ZeroDivisionError) as exc:
        print(f"FAIL: {exc}",file=sys.stderr);raise SystemExit(1)
