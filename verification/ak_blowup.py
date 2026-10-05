"""Reconstruction of ak_blowup from Table 5, paper page 56.

Only the checks reported by this module are implemented; see program_index.json.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from smale_verify.exact import (P, V, zero, equal, edge, equations, substitute,
    diff_poly, divide_poly, groebner, det, quotient_norm)

def verify_chart_and_boundary_fibre(D,m):
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    # Diagram 2 chart: derive its printed rows by substitution in the original
    # equations, then check the pair-block solution and boundary identities.
    t2=V('t2'); Z2c={j:V('D2Z'+str(j)) for j in m}; W2c={j:V('D2W'+str(j)) for j in m}
    D2c={e:V('D2D%d%d'%e) for e in D}; pair2={1,2}; tri2={3,4,5}
    z2raw={j:(t2**-1*Z2c[j] if j in pair2 else t2*Z2c[j]) for j in m}
    w2raw={j:(t2*W2c[j] if j in pair2 else t2**-1*W2c[j]) for j in m}
    d2raw={e:(D2c[e] if e==(1,2) or e[0] in tri2 else t2*D2c[e]) for e in D}
    Ez2,Ew2,T2raw=equations(m,z2raw,w2raw,d2raw)
    FZ2={j:(t2 if j in pair2 else t2**-1)*Ez2[j] for j in m}
    FW2={j:(t2**-1 if j in pair2 else t2)*Ew2[j] for j in m}
    for ai in pair2:
        aj=3-ai
        expected_z=(Z2c[ai]-m[aj]*D2c[(1,2)]**3*(Z2c[ai]-Z2c[aj])-
                    t2**3*sum((m[bj]*D2c[edge(ai,bj)]**3*(Z2c[ai]-t2**2*Z2c[bj]) for bj in tri2),P()))
        expected_w=(W2c[ai]-m[aj]*D2c[(1,2)]**3*(W2c[ai]-W2c[aj])+
                    t2*sum((m[bj]*D2c[edge(ai,bj)]**3*(W2c[bj]-t2**2*W2c[ai]) for bj in tri2),P()))
        equal(FZ2[ai],expected_z,'D2 chart pair z row '+str(ai))
        equal(FW2[ai],expected_w,'D2 chart pair w row '+str(ai))
    for bj in tri2:
        expected_z=(Z2c[bj]-sum((m[other]*D2c[edge(bj,other)]**3*(Z2c[bj]-Z2c[other]) for other in tri2 if other!=bj),P())+
                    t2*sum((m[ai]*D2c[edge(ai,bj)]**3*(Z2c[ai]-t2**2*Z2c[bj]) for ai in pair2),P()))
        expected_w=(W2c[bj]-sum((m[other]*D2c[edge(bj,other)]**3*(W2c[bj]-W2c[other]) for other in tri2 if other!=bj),P())-
                    t2**3*sum((m[ai]*D2c[edge(ai,bj)]**3*(W2c[bj]-t2**2*W2c[ai]) for ai in pair2),P()))
        equal(FZ2[bj],expected_z,'D2 chart triangle z row '+str(bj))
        equal(FW2[bj],expected_w,'D2 chart triangle w row '+str(bj))
    for e in D:
        i,j=e
        if e==(1,2):expected=D2c[e]**2*(Z2c[2]-Z2c[1])*(W2c[2]-W2c[1])-1
        elif i in tri2:expected=D2c[e]**2*(Z2c[j]-Z2c[i])*(W2c[j]-W2c[i])-1
        else:expected=D2c[e]**2*(t2**2*Z2c[j]-Z2c[i])*(W2c[j]-t2**2*W2c[i])-1
        equal(T2raw[e],expected,'D2 chart distance row '+str(e))
    M2=m[1]+m[2]; h2=V('D2_M_inverse'); d12_2=D2c[(1,2)]
    pair_eval={'D2Z1':-m[2]*h2,'D2Z2':m[1]*h2,
               'D2W1':-m[2]*d12_2,'D2W2':m[1]*d12_2}
    md3_2=M2*d12_2**3-1
    pair_targets=[m[2]*h2*md3_2,-m[1]*h2*md3_2,
                  m[2]*d12_2*md3_2,-m[1]*d12_2*md3_2,
                  (M2*h2)*(M2*d12_2**3-1)+(M2*h2-1)]
    pair_rows=[FZ2[1].coefficient('t2',0),FZ2[2].coefficient('t2',0),
                FW2[1].coefficient('t2',0),FW2[2].coefficient('t2',0),T2raw[(1,2)].coefficient('t2',0)]
    for i,(row,target) in enumerate(zip(pair_rows,pair_targets)):
        equal(substitute(row,pair_eval),target,'D2 pair boundary factorization '+str(i+1))
    equal(substitute(Z2c[2]-Z2c[1]-1,pair_eval),M2*h2-1,'D2 first gauge')
    equal(substitute(sum((m[j]*FZ2[j].coefficient('t2',0) for j in tri2),P()),{}),
          sum((m[j]*Z2c[j] for j in tri2),P()),'D2 triangle weighted centre z')
    equal(sum((m[j]*FW2[j].coefficient('t2',0) for j in tri2),P()),
          sum((m[j]*W2c[j] for j in tri2),P()),'D2 triangle weighted centre w')
    cross_ties=0
    for ai in pair2:
        for bj in tri2:
            equal(T2raw[edge(ai,bj)].coefficient('t2',0),
                  -D2c[edge(ai,bj)]**2*Z2c[ai]*W2c[bj]-1,
                  'D2 crossed boundary tie '+str((ai,bj)))
            cross_ties+=1
    dw3,dw4,dw5=map(V,('D2shapeW3','D2shapeW4','D2shapeW5'))
    ed45,ed35,ed34=dw5-dw4,dw3-dw5,dw4-dw3
    def Nshape(signs):
        s3,s4,s5=signs
        return (s3*dw3*ed35**2*ed34**2+s4*dw4*ed45**2*ed34**2+
                s5*dw5*ed45**2*ed35**2)
    S4d2=(dw3**3+dw4**3+dw5**3-dw3**2*dw4-dw3*dw4**2-dw3**2*dw5-
          dw3*dw5**2-dw4**2*dw5-dw4*dw5**2+3*dw3*dw4*dw5)
    S5d2=dw3**2+dw4**2+dw5**2-dw3*dw4-dw3*dw5-dw4*dw5
    shape_targets={(1,-1,-1):Nshape((1,-1,-1)),
                   (1,-1,1):Nshape((1,-1,1)),
                   (1,1,-1):Nshape((1,1,-1)),
                   (1,1,1):S4d2*S5d2}
    for signs in product((1,-1),repeat=3):
        representative=signs if signs[0]==1 else tuple(-s for s in signs)
        target=shape_targets[representative]
        equal(Nshape(signs),target if signs[0]==1 else -target,
              'D2 Euler sign shape factor '+str(signs))
    def omega_pow(j):return [(F(1),F(0)),(F(0),F(1)),(F(-1),F(-1))][j%3]
    def omega_add(x,y):return (x[0]+y[0],x[1]+y[1])
    def omega_neg(x):return (-x[0],-x[1])
    def omega_sub(x,y):return omega_add(x,omega_neg(y))
    def omega_mul(x,y):return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]-x[1]*y[1])
    def omega_scale(x,s):return (s*x[0],s*x[1])
    def omega_zero(x):return x==(F(0),F(0))
    lagrange_classes=fold_classes=equal_classes=two_equal_classes=0
    two_equal_patterns={}
    for a34,a35,a45 in product(range(3),repeat=3):
        Aom=omega_pow(a45);Bom=omega_sub(omega_add(omega_pow(a34),omega_pow(a45)),omega_pow(a35));Com=omega_pow(a34)
        disc=omega_sub(omega_mul(Bom,Bom),omega_scale(omega_mul(Aom,Com),4))
        lagrange_classes+=1
        if a34==a35==a45:
            equal_classes+=1
            if omega_zero(disc):raise AssertionError('D2 all-equal Lagrange class unexpectedly folds')
            if Bom!=Aom or Com!=Aom:raise AssertionError('D2 all-equal class does not give kappa^2+kappa+1')
            continue
        if len({a34,a35,a45})==3:
            fold_classes+=1
            if not omega_zero(disc):raise AssertionError(('D2 distinct Lagrange class is not a fold',(a34,a35,a45),disc))
            root=omega_scale(omega_mul(omega_neg(Bom),omega_pow(-a45)),F(1,2))
            root_poly=omega_add(omega_add(omega_mul(root,root),root),(F(1),F(0)))
            if not omega_zero(root_poly):raise AssertionError(('D2 fold root is not on S5',(a34,a35,a45),root_poly))
            continue
        two_equal_classes+=1
        repeated=(a34,a35) if a34==a35 else (a34,a45) if a34==a45 else (a35,a45)
        two_equal_patterns[tuple(sorted(repeated))]=two_equal_patterns.get(tuple(sorted(repeated)),0)+1
        for root in (omega_pow(1),omega_pow(2)):
            value=omega_add(omega_add(omega_mul(Aom,omega_mul(root,root)),omega_mul(Bom,root)),Com)
            if omega_zero(value):raise AssertionError(('D2 two-equal class unexpectedly lies on S5',(a34,a35,a45),root))
    if (lagrange_classes,fold_classes,equal_classes,two_equal_classes)!=(27,6,3,18):
        raise AssertionError('D2 Lagrange class count mismatch')
    if sorted(two_equal_patterns.values())!=[6,6,6]:raise AssertionError(('D2 two-equal pattern counts',two_equal_patterns))
    record('D2.6.2.chart_and_boundary_fibre',69,
           {'chart_coordinate_rows':10,'chart_distance_rows':10,'pair_block_factorizations':5,
            'triangle_centres':2,'crossed_boundary_ties':cross_ties,'Lagrange_classes':lagrange_classes,
            'fold_classes':fold_classes,'two_equal_classes':two_equal_classes,
            'two_equal_pattern_counts':{','.join(map(str,k)):v for k,v in two_equal_patterns.items()},'Euler_sign_shapes':8,
            'scope':'boundary and Lagrange fold/S5 classification checked; mass-resultant tables and interval enclosures remain'})
    return checks, {}

def verify_pair_second_order():
    a,d,h=map(V,('d2pair_a','d2pair_d','d2pair_h'))
    b,z1,z2,w1,w2=map(V,('d2pair_b','d2pair_z1','d2pair_z2','d2pair_w1','d2pair_w2'))
    rows=[z1-b*d**3*(z1-z2),w1-b*d**3*(w1-w2),
          w2-a*d**3*(w2-w1),d**2*(z2-z1)*(w2-w1)-1,z2-z1-1]
    # Eliminate M=a+b and M*d^3=1 without sampling either mass.
    point={'d2pair_b':d**-3-a,'d2pair_z1':-1+a*d**3,
           'd2pair_z2':a*d**3,'d2pair_w1':-(d**-3-a)*d,'d2pair_w2':a*d}
    variables=['d2pair_z1','d2pair_z2','d2pair_w1','d2pair_w2','d2pair_d']
    J=[[substitute(diff_poly(row,var),point) for var in variables] for row in rows]
    equal(det(J),-3*d**-1,'D2 selected pair determinant')
    response=[P(),P(),-h,-h,P()]
    forcing=[P(),-h,-h,P(),P()]
    for row,rhs in zip(J,forcing):
        equal(sum((x*y for x,y in zip(row,response)),P()),rhs,'D2 equal-forcing pair response')
    return [dict(check='D2.second_order.pair_block_and_equal_forcing',status='PASS',items=6,
        detail=dict(rows=['Fz1','Fw1','Fw2','T12','G1'],determinant='-3/D12',
            derivative=dict(Z1=0,Z2=0,W1='-h',W2='-h',D12=0),
            mass_relation='(m1+m2)D12^3=1',
            scope='the two computer-assisted pair-block steps in the second-order coefficient lemma; fold cases are separate'))],{}


if __name__ == "__main__":
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_blowup'))
