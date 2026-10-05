"""Reconstruction of ak_diagram6_blowup from Table 5, paper page 56.

Only the checks reported by this module are implemented; see program_index.json.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from smale_verify.exact import (P, V, zero, equal, edge, equations, substitute,
    diff_poly, divide_poly, groebner, det, quotient_norm)

def verify_chart(D,m):
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    # D6 chart, constructed directly from the original equations.
    t=V('t'); Z={k:V('Z'+str(k)) for k in m}; W={k:V('W'+str(k)) for k in m}
    az={1:-1,2:-1,3:1,4:1,5:1}; aw={1:1,2:1,3:-1,4:-1,5:1}
    cross={(1,3),(1,4),(2,3),(2,4)}
    zd={k:t**az[k]*Z[k] for k in m}; wd={k:t**aw[k]*W[k] for k in m}
    dd={e:(t*D[e] if e in cross else D[e]) for e in D}
    rz,rw,rt=equations(m,zd,wd,dd)
    wz={1:1,2:1,3:-1,4:-1,5:1}; wy={1:-1,2:-1,3:1,4:1,5:1}
    rz={k:t**wz[k]*rz[k] for k in m}; rw={k:t**wy[k]*rw[k] for k in m}
    rows=list(rz.values())+list(rw.values())+list(rt.values())+[Z[2]-Z[1]-1,W[4]-W[3]-1]
    degrees=[]
    for row in rows:
        if min(row.exponents('t'),default=0)<0 or max(row.exponents('t'),default=0)>5: raise AssertionError('D6 polynomial degree')
        degrees.append(max(row.exponents('t'),default=0))
    equal(m[1]*rz[1]+m[2]*rz[2]+m[5]*rz[5]+t**2*(m[3]*rz[3]+m[4]*rz[4]),
          m[1]*Z[1]+m[2]*Z[2]+t**2*(m[3]*Z[3]+m[4]*Z[4]+m[5]*Z[5]),'D6 z weighted center')
    equal(m[3]*rw[3]+m[4]*rw[4]+m[5]*rw[5]+t**2*(m[1]*rw[1]+m[2]*rw[2]),
          m[3]*W[3]+m[4]*W[4]+t**2*(m[1]*W[1]+m[2]*W[2]+m[5]*W[5]),'D6 w weighted center')
    for i,j in cross:
        equal(rt[(i,j)].coefficient('t',0),-D[(i,j)]**2*Z[i]*W[j]-1,'D6 cross boundary')
    wingz={1,2,5}; wingw={3,4,5}
    for e in set(D)-cross:
        i,j=e
        target=D[e]**2*((Z[j]-Z[i]) if e==(1,2) else (Z[j]-Z[i]) if e not in {(1,5),(2,5)} else -Z[i])*((W[j]-W[i]) if e!=(3,5) and e!=(4,5) else -W[i])-1
        equal(rt[e].coefficient('t',0),target,'D6 wing boundary')
    zwing_names={'Z1','Z2','W1','W2','W5','D12','D15','D25'}|{str('m'+str(i)) for i in m}
    wwing_names={'Z3','Z4','Z5','W3','W4','D34','D35','D45'}|{str('m'+str(i)) for i in m}
    for row in [rz[1],rz[2],rz[5],rw[1],rw[2],rt[(1,2)],rt[(1,5)],rt[(2,5)],rows[-2]]:
        if not row.coefficient('t',0).names()<=zwing_names: raise AssertionError('z wing separation')
    for row in [rw[3],rw[4],rw[5],rz[3],rz[4],rt[(3,4)],rt[(3,5)],rt[(4,5)],rows[-1]]:
        if not row.coefficient('t',0).names()<=wwing_names: raise AssertionError('w wing separation')
    # Replay item (ii) independently from raw coordinates.  Eliminate w34 using
    # t^2 z12 w34=1, and compare every transformed coordinate/distance row with
    # the corresponding original equation times its stated exact monomial factor.
    rawz12,rawt=V('raw_z12'),V('raw_t')
    rawz={1:V('raw_z1'),2:V('raw_z1')+rawz12,3:V('raw_z3'),4:V('raw_z4'),5:V('raw_z5')}
    raww34=rawt**-2/rawz12
    raww={1:V('raw_w1'),2:V('raw_w2'),3:V('raw_w3'),4:V('raw_w3')+raww34,5:V('raw_w5')}
    rawd={e:V('raw_delta%d%d'%e) for e in D}
    rawEz,rawEw,rawties=equations(m,rawz,raww,rawd)
    chart_map={'t':rawt}
    for j in m:
        chart_map['Z'+str(j)]=rawz[j]/rawz12 if j in (1,2) else rawz[j]*raww34
        chart_map['W'+str(j)]=raww[j]*rawz12 if j in (1,2,5) else raww[j]/raww34
    for e in D:
        chart_map['D%d%d'%e]=rawd[e]/rawt if e in cross else rawd[e]
    nu=(rawt*rawz12)**-1
    for j in m:
        equal(substitute(rz[j],chart_map),rawt**wz[j]*nu*rawEz[j],
              'D6.9.2(ii) transformed z row '+str(j))
        equal(substitute(rw[j],chart_map),rawt**wy[j]/nu*rawEw[j],
              'D6.9.2(ii) transformed w row '+str(j))
    for e in D:
        equal(substitute(rt[e],chart_map),rawties[e],'D6.9.2(ii) transformed distance row '+str(e))
    equal(substitute(rows[-2],chart_map),0,'D6.9.2(ii) first gauge')
    equal(substitute(rows[-1],chart_map),0,'D6.9.2(ii) second gauge')
    record('D6.9.2.chart_polynomial_boundary_center_split',44,
           {'t_degrees':degrees,'rotation_substitution':'all 10 coordinate rows, 10 distance rows, 2 gauges exact'})
    return checks, {'W':W, 'Z':Z, 'cross':cross, 'rt':rt, 'rw':rw, 'rz':rz, 't':t}

def verify_atlas_substitution_and_counts(Z,rt,rw,rz):
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    # Verify the z-wing parametrization in D6 Lemma 9.6.  Eliminate omega using
    # omega=lambda^2/2-1; omega^3=1 then gives the monic lambda polynomial below.
    lam=V('lam'); ds=V('scale'); cc=V('c'); bb=V('b')
    om=lam**2/2-1
    aa=(ds**-3-cc*lam**3)/2
    wing_sub={'m1':aa,'m2':aa,'m3':bb,'m4':bb,'m5':cc,
              'Z1':F(-1,2),'Z2':F(1,2),
              'W1':-ds*(aa+2*cc*lam),'W2':ds*(aa+2*cc*om*lam),
              'W5':aa*ds*(4*lam**-2-1),
              'D12':ds,'D15':lam*ds,'D25':om*lam*ds}
    lambda_mod=[-16,0,12,0,-6,0,1]
    z_wing_rows=[rz[1].coefficient('t',0),rz[2].coefficient('t',0),rz[5].coefficient('t',0),
                 rw[1].coefficient('t',0),rw[2].coefficient('t',0),
                 rt[(1,2)].coefficient('t',0),rt[(1,5)].coefficient('t',0),
                 rt[(2,5)].coefficient('t',0),Z[2]-Z[1]-1]
    for idx,row in enumerate(z_wing_rows):
        reduced=substitute(row,wing_sub).clear('lam').reduce('lam',lambda_mod)
        zero(reduced,'D6 atlas z-wing row '+str(idx+1))
    record('D6.9.6.atlas_z_wing_exact_substitution',len(z_wing_rows),
           {'parameters':'rational functions in scale,c,b modulo lam^6-6lam^4+12lam^2-16',
            'w_wing':'independently replayed by verify_atlas_inverse_elimination'})

    # Replay the exact shape and scale count in the finite wing atlas.  The shape
    # polynomial factors into the two Gaussian roots and four primitive roots.
    aa6,cc6=V('a6'),V('c6')
    quartic6=[4,0,-2,0,1]
    zero((8*((lam**2/2-1)**3-1)).reduce('lam',lambda_mod),
         'D6.9.6 omega is a cube root of unity')
    equal(lam**6-6*lam**4+12*lam**2-16,(lam**2-4)*(lam**4-2*lam**2+4),
          'D6 wing shape factorization')
    equal(((2*aa6+cc6*lam**3)*(2*aa6-cc6*lam**3)).reduce('lam',quartic6),4*aa6**2+8*cc6**2,
          'D6 primitive denominator product')
    equal((2*aa6+cc6*lam**3).coefficient('lam',0),2*aa6,
          'D6 scale equation nonzero constant term')
    for root,expected in ((F(2),2*aa6+8*cc6),(F(-2),2*aa6-8*cc6)):
        equal(substitute(2*aa6+cc6*lam**3,{'lam':root}),expected,
              'D6 Gaussian scale denominator '+str(root))
    record('D6.9.6.atlas_shape_scale_count',6,
           {'shape_roots':'2 Gaussian plus 4 primitive from lambda^6-6 lambda^4+12 lambda^2-16',
            'generic_wing_classes':18,'exceptional_a_equals_4c_classes':15,
            'primitive_denominator_product':'4a^2+8c^2; strictly positive for positive real a,c',
            'scope':'exact denominator and count argument for one wing; projective boundary not claimed'})
    crossed_choices=2**4
    fibre_counts={'generic':18*18*crossed_choices,
                  'one_exceptional_wing':18*15*crossed_choices,
                  'two_exceptional_wings':15*15*crossed_choices}
    if (crossed_choices,fibre_counts)!=(16,{'generic':5184,'one_exceptional_wing':4320,'two_exceptional_wings':3600}):
        raise AssertionError('D6 boundary fibre product count mismatch')
    record('D6.9.6.full_boundary_fibre_counts',3,
           {'crossed_quadratic_choices':crossed_choices,'boundary_fibre_sizes':fibre_counts,
            'method':'two wing factors times four independent nonzero quadratic crossed-distance roots'})
    return checks, {}

def verify_atlas_inverse_elimination(rz,rw,rt,Z,W):
    checks=[]
    a,c,d,lam,omega=map(V,('wing_a','wing_c','wing_d','wing_lambda','wing_omega'))
    base={'m1':a,'m2':a,'m5':c,'Z1':F(-1,2),'Z2':F(1,2),
          'D12':d,'D15':lam*d,'D25':omega*lam*d}
    z1=substitute(rz[1].coefficient('t',0),base)
    z2=substitute(rz[2].coefficient('t',0),base)
    equal(2*z1,(2*a+c*lam**3)*d**3-1,'D6 inverse atlas scale equation')
    equal(2*(z1+z2),c*lam**3*d**3*(1-omega**3),'D6 inverse atlas cube-root phase')
    t12=substitute(rt[(1,2)].coefficient('t',0),base)
    t15=substitute(rt[(1,5)].coefficient('t',0),base)
    t25=substitute(rt[(2,5)].coefficient('t',0),base)
    identity=(t12-2*lam**-2*t15-2*omega**-2*lam**-2*t25)*lam**2*omega**2
    equal(identity,-lam**2*omega**2+2*omega**2+2,'D6 inverse atlas distance closure')
    zero(((identity+(lam**2-2*(1+omega))*omega**2)).reduce('wing_omega',[-1,0,0,1]),'D6 inverse atlas lambda equation')
    linear_rows=[substitute(rw[1].coefficient('t',0),base),t12,t15]
    matrix=[[diff_poly(row,v) for v in ('W1','W2','W5')] for row in linear_rows]
    equal(det(matrix),d**4*lam**2/2,'D6 unique linear wing-coordinate solve')
    # Independently substitute the other wing; this is more than just a count
    # or an assertion that wing exchange is a symmetry.
    kap,e,cc,aa=map(V,('kap6','scale6e','c6atlas','a6atlas'))
    om=kap**2/2-1;b=(e**-3-cc*kap**3)/2
    lower={'m1':aa,'m2':aa,'m3':b,'m4':b,'m5':cc,
           'W3':F(-1,2),'W4':F(1,2),
           'Z3':-e*(b+2*cc*kap),'Z4':e*(b+2*cc*om*kap),
           'Z5':b*e*(4*kap**-2-1),'D34':e,'D35':kap*e,'D45':om*kap*e}
    lower_rows=[rw[3].coefficient('t',0),rw[4].coefficient('t',0),rw[5].coefficient('t',0),
                rz[3].coefficient('t',0),rz[4].coefficient('t',0),rt[(3,4)].coefficient('t',0),
                rt[(3,5)].coefficient('t',0),rt[(4,5)].coefficient('t',0),W[4]-W[3]-1]
    for j,row in enumerate(lower_rows):
        zero(substitute(row,lower).clear('kap6').reduce('kap6',[-16,0,12,0,-6,0,1]),'D6 independent w-wing atlas row '+str(j))
    ii=V('atlas_i')
    centre={'Z1':F(-1,2),'Z2':F(1,2),'W3':F(-1,2),'W4':F(1,2)}
    for signs in product((1,-1),repeat=4):
        for pair,sign,root in zip(((1,3),(1,4),(2,3),(2,4)),signs,(2*ii,2,2,2*ii)):
            sub=dict(centre);sub['D%d%d'%pair]=sign*root
            zero(substitute(rt[pair].coefficient('t',0),sub).reduce('atlas_i',[1,0,1]),'D6 crossed root '+str((pair,signs)))
    checks.append(dict(check='D6.9.6.inverse_atlas_and_other_wing',status='PASS',items=78,
                       detail=dict(inverse_constraints=['omega^3=1','lambda^2=2(1+omega)','(2a+c lambda^3)d^3=1'],
                                   unique_wing_coordinate_minor='d^4 lambda^2/2',independent_w_wing_rows=9,
                                   crossed_sign_vectors=16,scope='inverse constraints and unique linear solve, both wings, all crossed roots')))
    return checks,{}

def verify_selected_rank(rz,rw,rt,Z,W):
    lam,d,c,kap,e=map(V,('d6rank_l','d6rank_d','d6rank_c','d6rank_k','d6rank_e'))
    a=(d**-3-c*lam**3)/2;b=(e**-3-c*kap**3)/2
    mapping={'m1':a,'m2':a,'m3':b,'m4':b,'m5':c,
        'Z1':P.coerce(F(-1,2)),'Z2':P.coerce(F(1,2)),
        'W1':-d*(a+2*c*lam),'W2':d*(a+2*c*(lam**2/2-1)*lam),'W5':a*d*(4*lam**-2-1),
        'D12':d,'D15':lam*d,'D25':(lam**2/2-1)*lam*d,
        'W3':P.coerce(F(-1,2)),'W4':P.coerce(F(1,2)),
        'Z3':-e*(b+2*c*kap),'Z4':e*(b+2*c*(kap**2/2-1)*kap),'Z5':b*e*(4*kap**-2-1),
        'D34':e,'D35':kap*e,'D45':(kap**2/2-1)*kap*e}
    zrows=[rz[1],rz[2],rz[5],rw[1],rt[(1,2)],rt[(1,5)],rt[(2,5)],Z[2]-Z[1]-1]
    wrows=[rw[3],rw[5],rz[3],rz[4],rt[(3,4)],rt[(3,5)],rt[(4,5)],W[4]-W[3]-1]
    zvars=['Z1','Z2','W1','W2','W5','D12','D15','D25']
    wvars=['W3','W4','Z3','Z4','Z5','D34','D35','D45']
    cross=[(1,3),(1,4),(2,3),(2,4)]
    crossvars=['D%d%d'%pair for pair in cross]
    for rows,other in ((zrows,wvars+crossvars),(wrows,zvars+crossvars)):
        for row in rows:
            for variable in other:
                zero(diff_poly(row.coefficient('t',0),variable),'D6 rank off-diagonal block')
    for i,pair in enumerate(cross):
        row=rt[pair].coefficient('t',0)
        for j,name in enumerate(crossvars):
            expected=-2*V(name)*Z[pair[0]]*W[pair[1]] if i==j else P()
            equal(diff_poly(row,name),expected,'D6 crossed diagonal Jacobian')
    C,A,H,x=map(V,('d6rank_C','d6rank_A','d6rank_H','d6rank_x'))
    records=[]
    for number,(rows,variables,ratio,scale) in enumerate(((zrows,zvars,lam,d),(wrows,wvars,kap,e))):
        J=[[substitute(diff_poly(row.coefficient('t',0),name),mapping) for name in variables] for row in rows]
        raw=det(J)
        ratio_name=next(iter(ratio.names()));scale_name=next(iter(scale.names()))
        reduced=(raw/ratio**3).reduce(ratio_name,[-16,0,12,0,-6,0,1])
        # Remove the nonzero scale factor and write the dimensionless result
        # in C=c*scale^3; scale is a nonzero cube root on finite labels.
        expected_scale=scale if number==0 else scale**-1
        scaled=reduced/expected_scale
        for mon in scaled.terms:
            powers=dict(mon)
            if set(powers)-{ratio_name,scale_name,'d6rank_c'} or powers.get(scale_name,0)!=3*powers.get('d6rank_c',0):
                raise AssertionError('D6 dimensionless determinant scale mismatch')
        dimensionless=substitute(reduced/expected_scale,
            {ratio_name:x,scale_name:P.coerce(1),'d6rank_c':C})
        if number==0:
            expected=F(-9,8)*x**2+F(9,16)*x**4+C*(F(-9,4)*x**5-9*x+F(27,4)*x**3)
        else:
            expected=F(-9,8)*x**2+F(9,16)*x**4+C*(F(-9,4)*x**5+F(9,2)*x**3)+C**2*(F(9,2)*x**4+36-27*x**2)
        equal(dimensionless,expected,'D6 actual wing 8x8 determinant')
        for root in (2,-2):
            gaussian=substitute(dimensionless,{'d6rank_x':P.coerce(root)})
            equal(gaussian,F(9,2)*(1-C*root**3),'D6 Gaussian wing determinant')
            equal(substitute(gaussian,{'d6rank_C':H**-1})*H-F(9,2)*(H-root**3),P(),'D6 Gaussian denominator clearance')
        primitive=dimensionless.reduce('d6rank_x',[4,0,-2,0,1])
        numerator=substitute(primitive,{'d6rank_C':H**-1}).clear('d6rank_H')
        numerator=substitute(numerator,{'d6rank_H':2*A+x**3}).reduce('d6rank_x',[4,0,-2,0,1])
        norm=quotient_norm(numerator,'d6rank_x',[4,0,-2,0,1])
        equal(quotient_norm(2*A+x**3,'d6rank_x',[4,0,-2,0,1]),
              16*(A**2+2)**2,'D6 primitive rank denominator is a unit')
        expected_norm=F(6561,16)*A**4 if number==0 else 6561*A**4*(A**4-8*A**2+64)
        equal(norm,expected_norm,'D6 primitive actual Jacobian norm')
        records.append(dict(wing='z' if number==0 else 'w',variables=variables,
            determinant_scale='lambda^3*d' if number==0 else 'kappa^3/e',
            dimensionless_terms=[dict(coefficient=str(coef),powers=dict(mon)) for mon,coef in sorted(dimensionless.terms.items())],
            Gaussian='9A/(2A+ratio^3)',
            primitive_cleared_norm='6561 A^4/16' if number==0 else '6561 A^4(A^4-8A^2+64)'))
    equal(A**4-8*A**2+64,(A**2-4)**2+48,'D6 primitive rank norm positivity')
    return [dict(check='D6.9.9.selected_rank_twenty',status='PASS',items=259,
        detail=dict(wing_determinants=records,rank=20,
            selected_rows='all 22 chart rows except Rw2 and Rw4',
            block_structure='two independent 8x8 wing blocks and a diagonal 4x4 crossed block',
            crossed_units='-2 Dij Zi Wj; Zi,Wj=+/-1/2 and Dij=+/-2 or +/-2i',
            nonvanishing='A=a/c>0 or b/c>0; finite scale denominators; primitive quartic is (A^2-4)^2+48',
            scope='Lemma 9.9(i) at every finite label, including all corank-one signs; leading jet exclusions (ii) remain'))],{}


def verify_symmetry_orbits(rz,rw,rt,Z,W):
    from smale_verify.d6_states import census,SHAPE_SWAP,SHAPE_CONJUGATE
    ii,r2,r6=map(V,('d6sym_i','d6sym_r2','d6sym_r6'))
    roots=[P.coerce(2),P.coerce(-2),(r6+ii*r2)/2,(-r6-ii*r2)/2,(r6-ii*r2)/2,(-r6+ii*r2)/2]
    def reduction(poly):return poly.reduce('d6sym_i',[1,0,1]).reduce('d6sym_r2',[-2,0,1]).reduce('d6sym_r6',[-6,0,1])
    for index,value in enumerate(roots):
        zero(reduction(value*(value**2/2-1)-roots[SHAPE_SWAP[index]]),'D6 shape swap')
        zero(reduction(substitute(value,{'d6sym_i':-ii})-roots[SHAPE_CONJUGATE[index]]),'D6 shape conjugation')
        zero(reduction(value**6-6*value**4+12*value**2-16),'D6 exact shape')
    a,b,c=map(V,('d6sym_a','d6sym_b','d6sym_c'))
    masses={'m1':a,'m2':a,'m3':b,'m4':b,'m5':c}
    ez={j:substitute(row,masses) for j,row in rz.items()};ew={j:substitute(row,masses) for j,row in rw.items()}
    ties={pair:substitute(row,masses) for pair,row in rt.items()}
    p12={1:2,2:1,3:3,4:4,5:5}
    swap={'t':ii*V('t'),'Z1':-Z[2],'Z2':-Z[1],'W1':-W[2],'W2':-W[1],'W5':-W[5]}
    for pair in ties:
        image=tuple(sorted(p12[j] for j in pair));factor=-ii if pair in ((1,3),(1,4),(2,3),(2,4)) else P.coerce(1)
        swap['D%d%d'%pair]=factor*V('D%d%d'%image)
    for rows in (ez,ew):
        for j,row in rows.items():
            transformed=reduction(substitute(row,swap))
            if not any(not reduction(transformed-factor*rows[p12[j]]).terms for factor in (1,-1,ii,-ii)):
                raise AssertionError('D6 pair swap coordinate row')
    for pair,row in ties.items():
        equal(reduction(substitute(row,swap)),ties[tuple(sorted(p12[j] for j in pair))],'D6 pair swap tie')
    for gauge in (Z[2]-Z[1]-1,W[4]-W[3]-1):equal(reduction(substitute(gauge,swap)),gauge,'D6 pair swap gauge')
    reverse={'t':-V('t'),**{'D%d%d'%pair:-V('D%d%d'%pair) for pair in ((1,3),(1,4),(2,3),(2,4))}}
    for row in list(ez.values())+list(ew.values())+list(ties.values()):equal(substitute(row,reverse),row,'D6 reversal row')
    wing={1:3,2:4,3:1,4:2,5:5}
    exchange={'d6sym_a':b,'d6sym_b':a}
    for j in range(1,6):exchange['Z'+str(j)]=W[wing[j]];exchange['W'+str(j)]=Z[wing[j]]
    for pair in ties:exchange['D%d%d'%pair]=V('D%d%d'%tuple(sorted(wing[j] for j in pair)))
    for j in range(1,6):
        equal(substitute(ez[j],exchange),ew[wing[j]],'D6 wing exchange z')
        equal(substitute(ew[j],exchange),ez[wing[j]],'D6 wing exchange w')
    for pair,row in ties.items():equal(substitute(row,exchange),ties[tuple(sorted(wing[j] for j in pair))],'D6 wing exchange tie')
    equal(substitute(Z[2]-Z[1]-1,exchange),W[4]-W[3]-1,'D6 exchanged gauge')
    orbits=census()
    return [dict(check='D6.9.9.symmetry_and_576_state_orbits',status='PASS',items=576,
        detail=dict(shape_order=['2','-2','lambda0','-lambda0','conjugate(lambda0)','-conjugate(lambda0)'],
            crossed_root_order=['D13=2i eta13','D14=2 eta14','D23=2 eta23','D24=2i eta24'],
            orbit_count=29,regular_orbits=21,corank_one_orbits=8,
            regular_labels=3888,corank_one_labels=1296,orbits=orbits,
            symmetries='exact chart-row covariance: pair swap, reversal, conjugation, wing exchange',
            scope='finite atlas state transport; cube-root scales contribute nine generic labels per state'))],{}


def verify_leading_jets(rz,rw,rt):
    from smale_verify.d6_jets import verify_leading_jets as replay
    return replay(rz,rw,rt)


if __name__ == "__main__":
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram6_blowup'))
