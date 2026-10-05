"""Reconstruction of ak_diagram3_blowup from Table 5, paper page 56.

Only the checks reported by this module are implemented; see program_index.json.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from functools import lru_cache
from smale_verify.exact import (P, V, zero, equal, edge, equations, substitute,
    diff_poly, divide_poly, groebner, det, quotient_norm)

def verify_necessary_condition_identities():
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    # D3 first-condition syzygies and the exact mass eliminant from the order-3
    # obstruction.  These are polynomial identities before any specialization.
    ma,mb,md,me,ms=map(V,('d3_m1','d3_m2','d3_m4','d3_m5','d3_m3'))
    AA,BB,CC,DD=map(V,('d3_A','d3_B','d3_C','d3_D'))
    d14,d23,d12,d34,d13=map(V,('d3_D14','d3_D23','d3_D12','d3_D34','d3_D13'))
    fz1=AA-md*d14**3*(AA-BB); fz4=BB-ma*d14**3*(BB-AA)
    fz2=AA-ms*d23**3*(AA-BB); fz3=BB-mb*d23**3*(BB-AA)
    equal((ma+md)*(mb*fz2+ms*fz3)-(mb+ms)*(ma*fz1+md*fz4),
          (BB-AA)*(ma*ms-mb*md),'D3 first-condition z syzygy')
    fw1=CC-mb*d12**3*(CC-DD); fw2=DD-ma*d12**3*(DD-CC)
    fw3=DD-md*d34**3*(DD-CC); fw4=CC-ms*d34**3*(CC-DD)
    equal((ma+mb)*(ms*fw3+md*fw4)-(md+ms)*(ma*fw1+mb*fw2),
          (DD-CC)*(ma*ms-mb*md),'D3 first-condition w syzygy')
    g3=(BB-AA)-(DD-CC); T13=d13**2*(BB-AA)*(DD-CC)-1
    equal(d13**2*(BB-AA)**2-1,T13+d13**2*(BB-AA)*g3,'D3 gauge-distance identity')
    m3rel=mb*md/ma; M3=ma+mb+m3rel+md; Q3=me**2*M3
    Nre=((ma*mb*md)**3-Q3*(ma**6+mb**3*md**3))**2-4*Q3**2*ma**6*mb**3*md**3
    Nim=((ma*mb*md)**3-Q3*ma**3*(mb**3+md**3))**2-4*Q3**2*ma**6*mb**3*md**3
    S3=(ma**3-mb**3)*(ma**3-md**3)*(2*(ma*mb*md)**3-Q3*(ma**3+mb**3)*(ma**3+md**3))
    zero((Nre-Nim+Q3*S3).clear('d3_m1'),'D3 order-3 eliminant relation')
    equal((ma**3*(m3rel+md)-md**3*(ma+mb)).clear('d3_m1'),
          (ma*(ma+md)*(ma-md)*(m3rel+md)).clear('d3_m1'),'D3 order-4 cube identity 1')
    equal((ma**3*(mb+m3rel)-mb**3*(ma+md)).clear('d3_m1'),
          (ma*(ma+mb)*(ma-mb)*(mb+m3rel)).clear('d3_m1'),'D3 order-4 cube identity 2')
    record('D3.7.3_7.6.syzygies_and_mass_eliminant',6,
           {'items':['z and w first-condition syzygies','gauge-distance identity',
                     'order-3 mass eliminant relation','two order-4 cube identities'],
            'scope':'necessary-condition identities; companion stages verify the label Jacobian, jets and family census'})
    return checks, {}

def boundary_data():
    # m3=m2*m4/m1 is the first necessary mass condition.  Keep reciprocal
    # sums as monomial variables until denominators have been cleared.
    ma,mb,md,me,s12,s14=map(V,('m31','m32','m34','m35','s312','s314'))
    m={1:ma,2:mb,3:mb*md/ma,4:md,5:me}
    t=V('t3chart')
    names=('A','B','x','y','zhat5','C','D','u','v','what5')
    c={n:V('d3c_'+n) for n in names}
    A,B,x,y,zh,C,DD,u,v,wh=(c[n] for n in names)
    deltas={pair:V('d3c_D%d%d'%pair) for pair in combinations(range(1,6),2)}
    sides={(1,2),(1,4),(2,3),(3,4)}
    z={1:A/t,2:A/t+t*x,3:B/t,4:B/t+t*y,5:t**2*zh}
    w={1:C/t,2:DD/t,3:DD/t+t*v,4:C/t+t*u,5:t**2*wh}
    rawd={pair:value if pair in sides else t*value for pair,value in deltas.items()}
    rz,rw,ties=equations(m,z,w,rawd)
    rz={j:(t if j!=5 else t**-2)*row for j,row in rz.items()}
    rw={j:(t if j!=5 else t**-2)*row for j,row in rw.items()}
    gauge=t*(z[3]-z[1]-w[3]+w[1])
    point={'d3c_A':-md/s14,'d3c_B':ma/s14,'d3c_C':-mb/s12,'d3c_D':ma/s12,
           'd3c_x':s12*deltas[(1,2)],'d3c_y':-md*s12/ma*deltas[(3,4)],
           'd3c_u':s14*deltas[(1,4)],'d3c_v':mb*s14/ma*deltas[(2,3)]}
    point['d3c_zhat5']=-(ma*deltas[(1,5)]**3+mb*deltas[(2,5)]**3)*point['d3c_A']-(m[3]*deltas[(3,5)]**3+md*deltas[(4,5)]**3)*point['d3c_B']
    point['d3c_what5']=-(ma*deltas[(1,5)]**3+md*deltas[(4,5)]**3)*point['d3c_C']-(mb*deltas[(2,5)]**3+m[3]*deltas[(3,5)]**3)*point['d3c_D']
    root_powers={(1,2):(3,s12**-1),(1,4):(3,s14**-1),
                 (2,3):(3,ma/(mb*s14)),(3,4):(3,ma/(md*s12)),
                 (1,3):(2,P.coerce(1)),(2,4):(2,P.coerce(-1)),
                 (1,5):(2,s12*s14/(mb*md)),(2,5):(2,-s12*s14/(ma*md)),
                 (3,5):(2,s12*s14/ma**2),(4,5):(2,-s12*s14/(ma*mb))}
    def normal_form(poly):
        roots={'d3c_D%d%d'%pair:(power,value) for pair,(power,value) in root_powers.items()}
        result=P()
        for mon,coefficient in poly.terms.items():
            term=P.coerce(coefficient)
            for name,exponent in mon:
                if name in roots:
                    power,value=roots[name]
                    quotient,remainder=divmod(exponent,power)
                    term=term*V(name)**remainder*value**quotient
                else:term=term*V(name)**exponent
            result=result+term
        return result
    def reduce_label(poly):
        poly=normal_form(poly)
        poly=poly.clear('s312').clear('s314')
        return substitute(poly,{'s312':ma+mb,'s314':ma+md})
    def normalize(poly):
        poly=normal_form(poly)
        powers={name:max(0,-min(poly.exponents(name),default=0)) for name in ('s312','s314')}
        for name,n in powers.items():poly=poly*V(name)**n
        poly=substitute(poly,{'s312':ma+mb,'s314':ma+md})
        for name,n in powers.items():poly=poly*V(name)**-n
        return poly
    return dict(m=m,c=c,deltas=deltas,rz=rz,rw=rw,ties=ties,gauge=gauge,
                point=point,normal_form=normal_form,normalize=normalize,reduce_label=reduce_label,s12=s12,s14=s14,ma=ma,mb=mb,md=md,me=me)

def verify_chart_and_labels():
    data=boundary_data();m,c,D=data['m'],data['c'],data['deltas']
    rz,rw,ties=data['rz'],data['rw'],data['ties']
    A,B,CC,DD=c['A'],c['B'],c['C'],c['D']
    expected_z={1:A-m[4]*D[(1,4)]**3*(A-B),4:B-m[1]*D[(1,4)]**3*(B-A),
                2:A-m[3]*D[(2,3)]**3*(A-B),3:B-m[2]*D[(2,3)]**3*(B-A),
                5:c['zhat5']+(m[1]*D[(1,5)]**3+m[2]*D[(2,5)]**3)*A+(m[3]*D[(3,5)]**3+m[4]*D[(4,5)]**3)*B}
    expected_w={1:CC-m[2]*D[(1,2)]**3*(CC-DD),2:DD-m[1]*D[(1,2)]**3*(DD-CC),
                3:DD-m[4]*D[(3,4)]**3*(DD-CC),4:CC-m[3]*D[(3,4)]**3*(CC-DD),
                5:c['what5']+(m[1]*D[(1,5)]**3+m[4]*D[(4,5)]**3)*CC+(m[2]*D[(2,5)]**3+m[3]*D[(3,5)]**3)*DD}
    expected_ties={(1,2):D[(1,2)]**2*c['x']*(DD-CC)-1,(3,4):D[(3,4)]**2*c['y']*(CC-DD)-1,
                   (1,4):D[(1,4)]**2*c['u']*(B-A)-1,(2,3):D[(2,3)]**2*c['v']*(B-A)-1,
                   (1,3):D[(1,3)]**2*(B-A)*(DD-CC)-1,(2,4):D[(2,4)]**2*(B-A)*(CC-DD)-1,
                   (1,5):D[(1,5)]**2*A*CC-1,(2,5):D[(2,5)]**2*A*DD-1,
                   (3,5):D[(3,5)]**2*B*DD-1,(4,5):D[(4,5)]**2*B*CC-1}
    for j in m:
        equal(rz[j].coefficient('t3chart',0),expected_z[j],'D3 chart z boundary '+str(j))
        equal(rw[j].coefficient('t3chart',0),expected_w[j],'D3 chart w boundary '+str(j))
    for pair,row in ties.items():equal(row.coefficient('t3chart',0),expected_ties[pair],'D3 chart tie boundary '+str(pair))
    equal(data['gauge'].coefficient('t3chart',0),(B-A)-(DD-CC),'D3 boundary gauge')
    rows=list(rz.values())+list(rw.values())+list(ties.values())+[data['gauge']]
    for n,row in enumerate(rows):
        degrees=row.exponents('t3chart')
        if min(degrees,default=0)<0 or 1 in degrees or max(degrees,default=0)>6:
            raise AssertionError(('D3 polynomial chart degree',n,degrees))
        zero(data['reduce_label'](substitute(row.coefficient('t3chart',0),data['point'])),'D3 label row '+str(n))
    return [dict(check='D3.7.3.chart_and_label_substitution',status='PASS',items=42,
                 detail=dict(rows=21,degree_bound=6,no_linear_t=True,
                             root_classes=3**4*2**6,scope='explicit label substitution; inverse completeness follows by the linear boundary equations'))],{}

def selected_system(data):
    rz=data['rz'];rw=data['rw']
    m=data['m']
    rows=[m[1]*rz[1]+m[4]*rz[4],rz[4]-rz[1],rz[3]-rz[2],
          m[1]*rw[1]+m[2]*rw[2],rw[2]-rw[1],rw[3]-rw[4],rz[5],rw[5]]
    rows += [data['ties'][pair] for pair in sorted(data['ties'])]
    rows += [data['gauge']]
    return [substitute(row,{'d3c_B':data['c']['A']+1}) for row in rows]

def selected_boundary_system(data):
    return [row.coefficient('t3chart',0) for row in selected_system(data)]

def verify_label_jacobian():
    data=boundary_data();rows=selected_boundary_system(data)
    variables=['d3c_'+name for name in ('A','x','y','zhat5','C','D','u','v','what5')]+['d3c_D%d%d'%pair for pair in sorted(data['deltas'])]
    J=[[substitute(diff_poly(row,v),data['point']) for v in variables] for row in rows]
    actual=det(J)
    m,D=data['m'],data['deltas'];ma,mb,md=data['ma'],data['mb'],data['md']
    expected=-5184*ma**2*mb**3*md**3*(D[(1,2)]*D[(1,4)]*D[(2,3)]*D[(3,4)])**4
    expected=expected*D[(1,3)]*D[(2,4)]*D[(1,5)]*D[(2,5)]*D[(3,5)]*D[(4,5)]/(data['s12']*data['s14'])
    orientation=None
    for sign in (1,-1):
        if not data['reduce_label'](actual-sign*expected).terms:orientation=sign;break
    if orientation is None:raise AssertionError('D3 19x19 label Jacobian determinant mismatch')
    return [dict(check='D3.7.3.label_jacobian',status='PASS',items=1,
                 detail=dict(rank=19,variables=variables,
                             rows=['m1Fz1+m4Fz4','Fz4-Fz1','Fz3-Fz2','m1Fw1+m2Fw2','Fw2-Fw1','Fw3-Fw4','Fz5','Fw5']+['T%d%d'%pair for pair in sorted(D)]+['gauge'],
                             printed_formula_orientation=orientation,
                             nonvanishing='positive masses and all cubic/quadratic label roots nonzero'))],{}

def verify_order_three():
    data=boundary_data();m,D=data['m'],data['deltas'];rz,rw=data['rz'],data['rw']
    ellz=rz[1]-rz[2]-m[3]/m[2]*rz[3]+m[4]/m[1]*rz[4]
    ellw=rw[1]+m[2]/m[1]*rw[2]-m[2]/m[1]*rw[3]-rw[4]
    zero(ellz.coefficient('t3chart',0),'D3 left kernel z identity')
    zero(ellw.coefficient('t3chart',0),'D3 left kernel w identity')
    Omega=m[1]*m[2]*m[4]*(D[(1,3)]+D[(2,4)])+m[5]*(m[1]**2*D[(1,5)]+m[2]*m[4]*D[(3,5)]+m[1]*m[2]*D[(2,5)]+m[1]*m[4]*D[(4,5)])
    for name,row,factor in (('z',ellz,data['s12']/(m[1]**2*m[2])),('w',ellw,data['s14']/(m[1]**2*m[4]))):
        zero(data['reduce_label'](substitute(row.coefficient('t3chart',2),data['point'])),'D3 G2 '+name)
        zero(data['reduce_label'](substitute(row.coefficient('t3chart',3),data['point'])-factor*Omega),'D3 G3 '+name)
    return [dict(check='D3.7.5.order_three_obstruction',status='PASS',items=6,
                 detail=dict(G2='0',G3='((m1+m2)/(m1^2m2),(m1+m4)/(m1^2m4))*Omega',
                             scope='all cube/square-root sign classes by independent monic reductions'))],{}

@lru_cache(maxsize=1)
def compute_second_jet():
    data=boundary_data();m,D=data['m'],data['deltas'];rows=selected_system(data)
    variables=['d3c_'+n for n in ('A','x','y','zhat5','C','D','u','v','what5')]+['d3c_D%d%d'%pair for pair in sorted(D)]
    normalize=data['normalize']
    J=[[normalize(substitute(diff_poly(row.coefficient('t3chart',0),v),data['point'])) for v in variables] for row in rows]
    denominator=-5184*m[1]**2*m[2]**3*m[4]**3*(D[(1,2)]*D[(1,4)]*D[(2,3)]*D[(3,4)])**4
    denominator=denominator*D[(1,3)]*D[(2,4)]*D[(1,5)]*D[(2,5)]*D[(3,5)]*D[(4,5)]/(data['s12']*data['s14'])
    rhs=[normalize(-substitute(row.coefficient('t3chart',2),data['point'])) for row in rows]
    X2=[]
    for col in range(19):
        replaced=[[rhs[i] if j==col else J[i][j] for j in range(19)] for i in range(19)]
        X2.append(normalize(det(replaced)/denominator))
    for i in range(19):
        zero(data['reduce_label'](sum((J[i][j]*X2[j] for j in range(19)),P())-rhs[i]),'D3 second jet linear row '+str(i))
    return data,variables,X2

def verify_order_four():
    data,variables,X2=compute_second_jet();m,D=data['m'],data['deltas'];rz,rw=data['rz'],data['rw']
    ellz=rz[1]-rz[2]-m[3]/m[2]*rz[3]+m[4]/m[1]*rz[4]
    ellw=rw[1]+m[2]/m[1]*rw[2]-m[2]/m[1]*rw[3]-rw[4]
    Pi=(m[1]*D[(1,2)]-m[4]*D[(3,4)])*(m[1]*D[(1,4)]-m[2]*D[(2,3)])
    for name,ell,factor in (('z',ellz,m[4]*data['s12']),('w',ellw,m[2]*data['s14'])):
        ell=substitute(ell,{'d3c_B':data['c']['A']+1})
        zero(ell.coefficient('t3chart',0),'D3 vanishing left Hessian identity '+name)
        actual=substitute(ell.coefficient('t3chart',4),data['point'])
        actual=actual+sum((substitute(diff_poly(ell.coefficient('t3chart',2),v),data['point'])*q for v,q in zip(variables,X2)),P())
        zero(data['reduce_label'](actual+Pi/m[1]**2*factor),'D3 G4 '+name)
    witnesses={name:[dict(coefficient=str(coef),powers=dict(mon)) for mon,coef in sorted(value.terms.items())] for name,value in zip(variables,X2)}
    return [dict(check='D3.7.6.second_jet_and_order_four',status='PASS',items=23,
                 detail=dict(method='Cramer determinants with the certified monomial nonzero denominator; all 19 solved rows replayed',
                             G4='-Pi/m1^2 * (m4(m1+m2),m2(m1+m4))',
                             second_jet=witnesses,relations='same independent root powers and sum reciprocals as boundary_data'))],{}

def verify_mass_families_and_census():
    a,b,me,c,h=map(V,('d3f_a','d3f_b','d3f_m5','d3f_c','d3f_h'))
    results=[]
    for kind in ('E12','E14'):
        ma=a**2;mb=a**2 if kind=='E12' else b**2
        md=b**2 if kind=='E12' else a**2
        M=ma+mb+mb*md/ma+md;Q=me**2*M
        nr=((ma*mb*md)**3-Q*(ma**6+mb**3*md**3))**2-4*Q**2*ma**6*mb**3*md**3
        ni=((ma*mb*md)**3-Q*ma**3*(mb**3+md**3))**2-4*Q**2*ma**6*mb**3*md**3
        factor=a**12*(a**6*b**6-me**2*2*(a**2+b**2)*(a**3+b**3)**2)*(a**6*b**6-me**2*2*(a**2+b**2)*(a**3-b**3)**2)
        equal(nr,factor,'D3 mass-family real factorization '+kind)
        equal(ni,factor,'D3 mass-family imaginary factorization '+kind)
        for branch in (1,-1):
            substituted=substitute(nr,{'d3f_m5':a**3*b**3/(c*h)})
            substituted=substituted.clear('d3f_c').clear('d3f_h')
            substituted=substituted.reduce('d3f_c',[-2*(a**2+b**2),0,1])
            substituted=substitute(substituted,{'d3f_h':a**3+branch*b**3})
            zero(substituted,'D3 positive mass family converse '+kind+str(branch))
            # On each real/imaginary triple, the same two-term expression
            # applies, after dividing the imaginary component by i.
            for nu in ((1,) if branch==1 else (1,-1)):
                signs=[]
                for sigma,s1,s2 in product((1,-1),repeat=3):
                    ca=sigma*nu+s1;cb=branch*sigma*nu+s2
                    positive_me=a**3*b**3/(c*nu*h)
                    omega_real=ma*mb*md*sigma+positive_me*(ma**2*s1*c/b+mb*md*s2*c/a)
                    d25,d45=(s1*c/b,s2*c/a) if kind=='E12' else (s2*c/a,s1*c/b)
                    omega_imag=ma*mb*md*sigma+positive_me*(ma*mb*d25+ma*md*d45)
                    expected=a**4*b**2/(nu*h)*(sigma*nu*h+s1*a**3+s2*b**3)
                    equal(omega_real,expected,'D3 real sign form '+kind)
                    equal(omega_imag,expected,'D3 imaginary sign form '+kind)
                    signs.append(dict(signs=[sigma,s1,s2],coefficients=[ca,cb],vanishes_identically=ca==cb==0))
                    if ca and cb and ca*cb<0 and branch==1:
                        raise AssertionError('unexpected positive-family cancellation')
                    # In the minus family opposite nonzero coefficients
                    # have equal magnitudes, so a positive zero forces a=b.
                    if ca*cb<0 and abs(ca)!=abs(cb):
                        raise AssertionError('unexpected unequal cancellation')
                if sum(row['vanishes_identically'] for row in signs)!=2:
                    raise AssertionError('family sign census')
                results.append(dict(family=kind,branch=branch,nu=nu,sign_records=signs,labels=2*2*27))
    unequal_cubes=[q for q in product(range(3),repeat=4) if q[1]==q[2]]
    equal_cubes=[q for q in product(range(3),repeat=4) if q[0]==q[3] or q[1]==q[2]]
    both_equal=[q for q in equal_cubes if q[0]==q[3] and q[1]==q[2]]
    if (len(unequal_cubes),len(equal_cubes),len(both_equal))!=(27,45,9):
        raise AssertionError('D3 cube-label census')
    return [dict(check='D3.7.8.mass_family_factorizations',status='PASS',items=8,
        detail=dict(patterns=['E12+','E12-','E14+','E14-'],
            factorization='a^12 product_+,- (a^6b^6-m5^2 c^2(a^3 +/- b^3)^2)',
            c_square='2(a^2+b^2)',positive_m5='a^3b^3/(c |a^3 +/- b^3|)',
            scope='factorization gives both positive possibilities; minus branch is impossible at a=b')),
        dict(check='D3.7.9.all_family_label_counts',status='PASS',items=5184,
        detail=dict(families=results,unequal_labels=108,equal_labels=180,
            unequal_cube_indices=unequal_cubes,equal_cube_indices=equal_cubes,
            equal_sign_relations=['D15=D35=-2D13','D25=D45=-2D24'],
            justification='two vanishing sign choices in each triple; unequal-mass alternative cube equality is excluded by the checked cube-difference identities'))],{}


if __name__ == "__main__":
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram3_blowup'))
