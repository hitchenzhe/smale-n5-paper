"""Reconstruction of ak_diagram1_edge_triangle from Table 5, paper page 56.

Only the checks reported by this module are implemented; see program_index.json.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from smale_verify.exact import (P, V, zero, equal, edge, equations, substitute,
    diff_poly, divide_poly, groebner, det, quotient_norm)

def verify_morse_bott_triangle():
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    # D1 Lemma 5.9: specialize to the triangle-mass parametrization b=x,c=1,
    # kappa=1.  Translation is left free, so the 9-variable Jacobian has two
    # translation directions and two scaling directions in its kernel.
    xpar,rpar=V('xparam'),V('rparam')
    mm3=xpar**2;mm4=xpar**2*rpar**-2;mm5=rpar**-2
    xx={3:V('x3'),4:V('x4'),5:V('x5')}; yyv={3:V('y3'),4:V('y4'),5:V('y5')}
    ddv={e:V('d%d%d'%e) for e in combinations((3,4,5),2)}
    mmtri={3:mm3,4:mm4,5:mm5}
    ph={j:sum(mmtri[l]*ddv[edge(j,l)]**3*(xx[j]-xx[l]) for l in (3,4,5) if j!=l) for j in (3,4,5)}
    ps={j:sum(mmtri[l]*ddv[edge(j,l)]**3*(yyv[j]-yyv[l]) for l in (3,4,5) if j!=l) for j in (3,4,5)}
    ch={e:ddv[e]**2*(xx[e[1]]-xx[e[0]])*(yyv[e[1]]-yyv[e[0]])-1 for e in ddv}
    mb_rows=[ph[3],ph[4],ps[3],ps[4],ch[(3,4)],ch[(3,5)],ch[(4,5)]]
    mb_vars=['x3','x4','x5','y3','y4','y5','d34','d35','d45']
    mb_base={'x3':P(),'x4':rpar,'x5':-rpar/xpar,'y3':P(),'y4':rpar,'y5':-rpar/xpar,
              'd34':rpar**-1,'d35':xpar/rpar,'d45':-xpar}
    rp_relation=xpar*rpar+rpar-1
    rp_order=('rparam','xparam')
    for i,row in enumerate(mb_rows):
        zero(divide_poly(substitute(row,mb_base).clear('xparam').clear('rparam'),[rp_relation],rp_order),
             'D1.5.9 base point row '+str(i))
    mb_jac=[[substitute(diff_poly(row,v),mb_base) for v in mb_vars] for row in mb_rows]
    for rownum,row in enumerate(mb_jac):
        for which,vector in enumerate([
            [1,1,1,0,0,0,0,0,0],[0,0,0,1,1,1,0,0,0],
            [mb_base['x3'],mb_base['x4'],mb_base['x5'],0,0,0,-mb_base['d34']/2,-mb_base['d35']/2,-mb_base['d45']/2],
            [0,0,0,mb_base['y3'],mb_base['y4'],mb_base['y5'],-mb_base['d34']/2,-mb_base['d35']/2,-mb_base['d45']/2]]):
            residual=sum((entry*vec for entry,vec in zip(row,vector)),P()).clear('xparam').clear('rparam')
            zero(divide_poly(residual,[rp_relation],rp_order),'D1.5.9 kernel row '+str(rownum)+'.'+str(which))
    positive_poly=None;chosen_minor=None
    for row_ids in combinations(range(7),5):
        for col_ids in combinations(range(9),5):
            candidate=det([[mb_jac[i][j] for j in col_ids] for i in row_ids]).clear('xparam').clear('rparam')
            deg_r=max(candidate.exponents('rparam'),default=0);converted=P()
            for mon,coef in candidate.terms.items():
                ex=dict(mon);kpow=ex.get('rparam',0)
                rest=tuple((v,n) for v,n in mon if v!='rparam')
                converted=converted+P({rest:coef})*(V('xparam')+1)**(deg_r-kpow)
            coeff={dict(mon).get('xparam',0):c for mon,c in converted.terms.items()}
            if coeff and (all(c>0 for c in coeff.values()) or all(c<0 for c in coeff.values())):
                positive_poly=converted;chosen_minor=(row_ids,col_ids);break
        if positive_poly is not None:break
    if positive_poly is None:raise AssertionError('D1 Morse-Bott: no sign-definite rank-five minor found')
    zero(divide_poly(rpar**deg_r*positive_poly-candidate,[rp_relation],rp_order),
         'D1.5.9 cleared minor equals substituted polynomial modulo r(x+1)=1')
    pair_moment=mm3*mm4*mb_base['x4']**2+mm3*mm5*mb_base['x5']**2+mm4*mm5*(mb_base['x4']-mb_base['x5'])**2
    closed_moment=mm3*mm4/mm5+mm3*mm5/mm4+mm4*mm5/mm3
    zero(divide_poly((pair_moment-closed_moment).clear('xparam').clear('rparam'),[rp_relation],rp_order),'D1.5.9 inertia formula')
    coeff_list=[[dict(mon).get('xparam',0),str(c)] for mon,c in sorted(positive_poly.terms.items(),key=lambda it:dict(it[0]).get('xparam',0))]
    record('D1.5.9.morse_bott_rank_kernel_inertia',37,{'minor_rows':[i+1 for i in chosen_minor[0]],'minor_columns':[mb_vars[j] for j in chosen_minor[1]],'sign':'constant for x>0','minor_coefficients_x':coeff_list,'inertia_formula':'verified by exact pairwise-distance identity'})
    return checks, {}

if __name__ == "__main__":
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram1_edge_triangle'))
