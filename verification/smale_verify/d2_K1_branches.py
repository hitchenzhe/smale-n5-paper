"""Complete K1 common-root and positive real-branch certificate.

Joint Bernstein exclusions prove the first linear subresultant never
vanishes identically on the positive K1 locus. The two branches are given
by an exact algebraic parametrization with a unique positive shape root;
the proof does not depend on a floating root tracker.
"""
from pathlib import Path
import json
from fractions import Fraction as Q
from .symbolic import sympy
from .bernstein import corner_sectors,coefficients,split_axis,blowup_critical_corner
from .d2_mass_rays import path_box

def joint_cover(matrices,face=False,max_depth=32):
    pending=[('',matrices)];leaves=[]
    while pending:
        path,ms=pending.pop();chosen=None;sign=None
        k_lo=path_box(path)[1][0]if face else Q(1)
        for i,matrix in enumerate(ms):
            ss=1 if matrix[-1][-1]>0 else -1
            minimum=min(ss*q for row in matrix for q in row)
            if minimum>=0 and ss*matrix[-1][-1]>0 and(k_lo>0 or ss*matrix[-1][0]>0):
                chosen=i;sign=ss;break
        if chosen is not None:
            mat=ms[chosen];leaves.append(dict(path=path,nonzero_polynomial=chosen,sign=sign,minimum=str(min(sign*q for row in mat for q in row)),upper_corner=str(mat[-1][-1]),positive_mass_face=str(mat[-1][0])if k_lo==0 else None))
        elif len(path)>=max_depth:raise ArithmeticError('K1 joint exclusion has an unresolved box')
        else:
            splits=[split_axis(matrix,len(path)%2)for matrix in ms]
            pending.extend([(path+'0',[q[0]for q in splits]),(path+'1',[q[1]for q in splits])])
    return leaves

def verify(data):
    from ak_resultants import shape_factors
    s=sympy();a,b,x=s.symbols('a b x');root=Path(__file__).resolve().parents[1]
    def decode(ts):return sum(s.Rational(t['coefficient'])*a**t['powers'].get('m3',0)*b**t['powers'].get('m4',0)*x**t['powers'].get('x',0)for t in ts)
    F=decode(data['resultants']['1']['shape_flat']);G=decode([t for t in data['P_flat']if not t['powers'].get('y',0)]);K1=decode(data['resultants']['1']['K_terms'])
    subs=s.subresultants(F,G,x)
    if s.degree(subs[-2],x)!=1 or s.expand(subs[-1]-a**10*K1):raise ArithmeticError('K1 source subresultant/resultant reconstruction')
    linear=s.Poly(subs[-2],x);C1=s.cancel(linear.nth(1)/a**4);C0=s.cancel(linear.nth(0)/a**4)
    if s.Poly(G,x).LC()!=b*b:raise ArithmeticError('K1 common root at infinity is not excluded')
    polys=[C1,C0,K1];sectors=[];encoded=[]
    for p in polys:
        degree=int(s.total_degree(p));terms=[dict(coefficient=str(c),powers=dict(m3=int(pa),m4=int(pb),m5=degree-int(pa)-int(pb)))for (pa,pb),c in s.Poly(p,a,b).terms()]
        encoded.append(terms);sectors.append(list(corner_sectors(terms)))
    covers=[]
    for sector in range(6):
        metadata={name:value for name,value in sectors[0][sector].items()if name!='power_terms'}
        if sector in(0,2):
            matrices=[coefficients(sectors[i][sector]['power_terms'])for i in range(3)]
            leaves=joint_cover(matrices);covers.append(dict(sector=metadata,mode='direct',leaves=leaves))
        else:
            for mode in(0,1):
                transformed=[blowup_critical_corner(sectors[i][sector]['power_terms'],mode)for i in range(3)]
                matrices=[coefficients(p)for p,power in transformed];leaves=joint_cover(matrices,face=mode==1)
                covers.append(dict(sector=metadata,mode=mode,removed_h_powers=[int(power)for p,power in transformed],leaves=leaves))
        print('D2 K1: joint subresultant cover sector',sector+1,'/6 complete',flush=True)
    # Real common roots must have W3,W4<0 by the real P0 sign identity.
    r,c,t,z=s.symbols('r c t z');mapping={'shape_w3':-r*r,'shape_w4':-s.Integer(1),'shape_w5':c*c}
    shape=s.expand(sum(s.Rational(q.numerator,q.denominator)*s.prod(mapping[n]**power for n,power in mon)for mon,q in shape_factors()[0][1].terms.items()))
    normal=s.Poly(shape.subs(c*c,z).subs(r*r,t),z)
    # No shape with |W3|>=|W4|: the translated polynomial has
    # strictly negative coefficients, including its pure z coefficient.
    u=s.symbols('u')
    ordered=s.Poly(s.expand(normal.as_expr().subs(t,1+u)),u,z)
    if any(q>=0 for q in ordered.coeffs()) or ordered.coeff_monomial(z)>=0:
        raise ArithmeticError('K1 strict shape ordering')
    expected_coefficients={0:(1-t)**3*(t*t+t+1),3:-t*t-4*t-1,4:-2*(t+1),5:-1}
    if any(s.expand(normal.nth(j)-q) for j,q in expected_coefficients.items()):
        raise ArithmeticError('K1 Descartes endpoint/high coefficient identities')
    Csecond=normal.nth(1);Cthird=normal.nth(2)
    positive=t*(6-3*t+6*t*t-2*t**3)
    if s.expand(Csecond-2*Cthird-positive):raise ArithmeticError('K1 Descartes coefficient implication')
    if min(row[0]for row in coefficients({(int(power[0]),0):Q(c)for power,c in s.Poly(6-3*t+6*t*t-2*t**3,t).terms()}))<=0:raise ArithmeticError('K1 coefficient implication positivity')
    critical=1-2*t*t-8*t**3+t**4-8*t**5
    if s.expand(shape.subs(c,r)-critical.subs(t,r*r)):raise ArithmeticError('K1 parameter threshold polynomial')
    if s.expand(s.diff(critical,t)+4*t*(1-t*t+6*t+10*t**3)):raise ArithmeticError('K1 threshold strict monotonicity')
    if s.Poly(critical,t).count_roots(0,1)!=1 or s.Poly(critical,t).count_roots(s.Rational(2,5),s.Rational(9,20))!=1:raise ArithmeticError('K1 threshold root isolation')
    if s.expand(shape.subs(c,1)+r*r*((r**4-1)**2+16)):raise ArithmeticError('K1 upper shape bound c<1')
    mass_branches=[[(1-c**3)*r,c**3-r**3,(1-r**3)*c],[(1+c**3)*r,c**3-r**3,(1+r**3)*c]]
    branch_records=[];immersion_records=[]
    for mass,signs in zip(mass_branches,((1,1,-1),(1,-1,-1))):
        if s.expand(mass[0]*r*r+mass[1]-mass[2]*c*c)or s.cancel(signs[0]*mass[0]/r+signs[1]*mass[1]+signs[2]*mass[2]/c):raise ArithmeticError('K1 positive mass branch equations')
        determinant=r*r*signs[1]-signs[0]/r
        branch_records.append(dict(mass=list(map(str,mass)),signs=list(signs),mass_solve_determinant=str(s.factor(determinant))))
    # Prove immersion directly in the algebraic shape coordinates. This
    # does not assume that the Euler normalization admits a finite scale.
    # Along shape=0, D=shape_c*d/dr-shape_r*d/dc is a nonzero tangent.
    # A stationary projective mass map would annul both ratio numerators.
    for branch,(mass,sign) in enumerate(zip(mass_branches,(-1,1))):
        derivatives=[];resultants=[]
        for m in mass[:2]:
            h=s.expand(s.diff(shape,c)*(s.diff(m,r)*mass[2]-m*s.diff(mass[2],r))
                -s.diff(shape,r)*(s.diff(m,c)*mass[2]-m*s.diff(mass[2],c)))
            derivatives.append(h);resultants.append(s.resultant(shape,h,c))
        common=s.Poly(s.gcd(*resultants),r).monic().as_expr()
        left,right=(30,28) if sign==-1 else (28,30)
        quadratic=r*r-sign*r+1
        expected=r**10*(r-1)**left*(r+1)**right*quadratic**2
        if s.expand(common-expected):
            raise ArithmeticError('K1 projective mass immersion resultant gcd')
        if s.expand(quadratic-((r-s.Rational(sign,2))**2+s.Rational(3,4))):
            raise ArithmeticError('K1 immersion quadratic positivity')
        immersion_records.append(dict(branch=branch,
            tangent_numerators=[str(q)for q in derivatives],
            resultant_coefficients=[[str(q)for q in s.Poly(v,r).all_coeffs()]for v in resultants],
            monic_resultant_gcd=str(s.factor(common)),
            nonzero_domain='0<r<1; the quadratic is (r-sign/2)^2+3/4'))
    payload=dict(check='D2.6.16.K1_unique_real_common_root_and_two_branches',status='PASS',subresultant_polynomials=encoded,
        joint_positive_mass_covers=covers,shape_equation=str(normal.as_expr()),threshold_polynomial=str(critical),threshold_interval=['2/5','9/20'],branches=branch_records,
        projective_immersion=immersion_records,
        parametrization='0<r<sqrt(t0); c is the positive square root of the unique positive root z of the displayed degree-five polynomial. Then r<c<1. Positive mass scaling is free.',
        proof=['The nonzero flat P0 leading coefficient excludes common roots at infinity. At K1=0 a finite common root exists by the exact mass-resultant definition.',
            'The six ordered mass sectors are covered directly or by two complementary blow-up charts. At every positive mass point at least one of C1,C0,K1 is certified nonzero. The k=0 genuine positive-mass faces of the second blow-up chart are checked explicitly. Thus on K1=0 the first linear subresultant is not identically zero.',
            'A finite common root must annul that linear polynomial, so C1 cannot be zero there and the common root is unique. Real mass coefficients force this unique common root to be real.',
            'The real P0 identity (f3+f4+f5)^2-4f3*f4 forces W3,W4<0. Put W=(-a^2,-b^2,c^2). The S1 equation is negative when a>=b, so normalize b=1 and set 0<r=a/b<1.',
            'The square-coordinate shape polynomial has exactly one positive root by Descartes signs: its constant is positive, its three highest coefficients are negative, and Csecond-2*Cthird is strictly positive on 0<t<1. Hence Cthird>=0 implies Csecond>0. This also proves the positive shape root is simple.',
            'Evaluation at c=1 proves c<1. Positive solutions of the signed cancellation and mass-centre equations require c>r. This occurs exactly when the strictly decreasing threshold polynomial is positive, namely r^2<t0. The two displayed positive mass kernels exhaust the four relative sign choices.',
            'Both mass-solve determinants are nonzero on 0<r<1. The shape curve is smooth because its positive z root is simple. The two projective mass-ratio tangent numerators cannot vanish together on the shape curve: their exact resultants have gcd r^10(r-1)^30(r+1)^28(r^2+r+1)^2 on branch I and r^10(r-1)^28(r+1)^30(r^2-r+1)^2 on branch II. Every factor is nonzero on 0<r<1. This directly proves immersion without assuming a finite normalized Euler scale. Unique common-root elimination, local mass/shape reconstruction, and the single vanishing signed P0 factor show the positive K1 branches are smooth.',
            'The K2 case is the verified body exchange. Regularity and qprime!=0 on both entire branches come from the stronger global Euler response certificate, so the original one-parameter response interval cover is unnecessary.'],
        depends_on=['D2.6.17.Euler_global_spectral_exclusion','D2.6.17.all_regular_X_unit_circle_exclusions','D2.6.3_6.5.mass_resultants_and_tables'],
        scope='Complete unique real common-root and two-branch classification with an exact algebraic parameter, plus a separate stronger global response proof. This replaces the unspecified original rational parameter and branch interval cover.')
    output=root/'certificates/d2_K1_branches.json';output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=2,detail=dict(certificate=output.name,unique_real_common_root=True,complete_positive_branches=True))],{}
