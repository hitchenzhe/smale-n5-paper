"""Paper Table 5: ak_diagram8_relation40_derivation.

Reconstruct all 22 ratio-chart rows from the original five-body equations.
Verify every polynomial statement of Lemma 11.1 and the first mass eliminations.
"""
from fractions import Fraction as F
from smale_verify.exact import P, V, zero, equal, equations, substitute, normalized_system


def chart_rows(m):
    a,eta=V('a8'),V('eta8')
    names=('A1','A2','C','B3','B4','D','u','v','R12','R25','P','Q')
    x={name:V('d8_'+name) for name in names}
    A1,A2,C,B3,B4,D,u,v,R12,R25,pp,qq=(x[name] for name in names)
    nu=a**2*eta**2
    z={1:A1,2:A2,3:C+m[5]*nu*u,4:C+m[5]*nu*v,
       5:C-nu*(m[3]*u+m[4]*v)}
    w={1:a**-2*eta**-2*D+m[5]*a**4*pp,
       2:a**-2*eta**-2*D+m[5]*a**4*qq,
       3:a**-2*eta**-2*B3,4:a**-2*eta**-2*B4,
       5:a**-2*eta**-2*D-a**4*(m[1]*pp+m[2]*qq)}
    deltas={(1,2):a**-2*R12,(1,5):a**-2,(2,5):a**-2*R25}
    for pair in ((1,3),(1,4),(2,3),(2,4),(3,4),(3,5),(4,5)):
        name='D%d%d'%pair
        x[name]=V('d8_'+name)
        deltas[pair]=a*eta*x[name] if pair[0] in (1,2) else x[name]
    rz,rw,ties=equations(m,z,w,deltas)
    E={}
    for k in m:
        E[2*k-2]=(a**6 if k in (1,2,5) else 1)*rz[k]
        E[2*k-1]=a**2*eta**2*rw[k]
    for j,pair in enumerate(sorted(deltas),10): E[j]=ties[pair]
    E[20]=A2-A1-1
    E[21]=B4-B3-1
    X,Y=A1-C,A2-C
    f={0:-(m[2]*(X-Y)*R12**3+m[5]*X),
       2:-(m[1]*(Y-X)*R12**3+m[5]*Y*R25**3),
       8:m[1]*X+m[2]*Y*R25**3,
       20:A2-A1-1,10:m[5]*(X-Y)*(pp-qq)*R12**2-1,
       13:X*((m[1]+m[5])*pp+m[2]*qq)-1,
       16:Y*R25**2*(m[1]*pp+(m[2]+m[5])*qq)-1,
       4:C,6:C,1:D,3:D,
       11:X*(D-B3)*x['D13']**2-1,12:X*(D-B4)*x['D14']**2-1,
       14:Y*(D-B3)*x['D23']**2-1,15:Y*(D-B4)*x['D24']**2-1,
       5:B3+m[4]*(B4-B3)*x['D34']**3-m[5]*(B3-D)*x['D35']**3,
       7:B4-m[3]*(B4-B3)*x['D34']**3-m[5]*(B4-D)*x['D45']**3,
       9:D+m[3]*(B3-D)*x['D35']**3+m[4]*(B4-D)*x['D45']**3,
       21:B4-B3-1,17:m[5]*(B3-B4)*(u-v)*x['D34']**2-1,
       18:(B3-D)*((m[3]+m[5])*u+m[4]*v)*x['D35']**2-1,
       19:(B4-D)*(m[3]*u+(m[4]+m[5])*v)*x['D45']**2-1}
    beta1=-m[5]*(m[2]*(pp-qq)*R12**3+(m[1]+m[5])*pp+m[2]*qq)
    beta3=-m[5]*(m[1]*(qq-pp)*R12**3+(m[1]*pp+(m[2]+m[5])*qq)*R25**3)
    tau1=m[3]*(B3-D)*x['D13']**3+m[4]*(B4-D)*x['D14']**3
    tau3=m[3]*(B3-D)*x['D23']**3+m[4]*(B4-D)*x['D24']**3
    return dict(E=E,f=f,x=x,a=a,eta=eta,nu=nu,beta1=beta1,beta3=beta3,tau1=tau1,tau3=tau3,
                z=z,w=w,deltas=deltas)


def verify_polynomial_faces(m):
    data=chart_rows(m)
    E,f,x,a,eta,nu=(data[n] for n in ('E','f','x','a','eta','nu'))
    for j in range(22):
        if min(E[j].exponents('a8'),default=0)<0 or min(E[j].exponents('eta8'),default=0)<0:
            raise AssertionError(('non-polynomial D8 row',j))
        equal(E[j].coefficient('a8',0).coefficient('eta8',0),f[j],'D8 face '+str(j))
    M=m[3]+m[4]+m[5]
    equal(m[1]*E[0]+m[2]*E[2]+m[5]*E[8]+a**6*(m[3]*E[4]+m[4]*E[6]),
          a**6*(m[1]*x['A1']+m[2]*x['A2']+M*x['C']),'D8 mass-centre identity')
    for j,beta,tau in ((1,data['beta1'],data['tau1']),(3,data['beta3'],data['tau3'])):
        remainder=E[j]-x['D']-eta**2*beta-a**3*eta**3*tau
        for mon in remainder.terms:
            ex=dict(mon)
            if ex.get('a8',0)<6 or ex.get('eta8',0)<2:
                raise AssertionError(('D8 w-row remainder exponent',j,mon))
    difference=E[1]-E[3]
    exponents={tuple(dict(mon).get(v,0) for v in ('a8','eta8')) for mon in difference.terms}
    if exponents!={(0,2),(3,3),(6,2),(9,5)}:
        raise AssertionError(('D8 E1-E3 monomial census',exponents))
    correction=nu*(m[3]*x['u']+m[4]*x['v'])
    mapping={'d8_C':x['C']-correction}
    for j in (0,2,8,10,13,16,20):
        remainder=E[j]-substitute(f[j],mapping)
        if any(dict(mon).get('a8',0)<6 for mon in remainder.terms):
            raise AssertionError(('D8 upper invariant remainder',j))
    L3=-m[5]*(m[4]*(x['u']-x['v'])*x['D34']**3+
              ((m[3]+m[5])*x['u']+m[4]*x['v'])*x['D35']**3-x['u'])
    F3=-(m[1]*(x['C']-x['A1'])*x['D13']**3+m[2]*(x['C']-x['A2'])*x['D23']**3)
    G3=-m[5]*x['u']*(m[1]*x['D13']**3+m[2]*x['D23']**3)
    equal(E[4],x['C']+nu*L3+a**3*eta**3*F3+a**5*eta**5*G3,'D8 exact E4 decomposition')
    return [dict(check='D8.11.1.all_ratio_chart_rows',status='PASS',items=34,
                 detail=dict(rows=22,polynomial=True,faces=22,mass_centre=True,
                             w_difference_exponents=sorted(exponents),upper_invariant_rows=7,
                             scope='all five polynomial assertions in Lemma 11.1; higher order ranks and exclusion lemmas remain'))],{}


def verify_first_mass_eliminations(m):
    data=chart_rows(m)
    f,x=data['f'],data['x']
    X,Y=x['A1']-x['C'],x['A2']-x['C']
    R12,R25=x['R12'],x['R25']
    lhs=(X-Y)*R12**2*X*f[16]-(X-Y)*R12**2*Y*R25**2*f[13]+X*Y*R25**2*f[10]
    rhs=(X-Y)*R12**2*(Y*R25**2-X)-X*Y*R25**2
    equal(lhs,rhs,'D8 elimination of P,Q')
    M12=m[1]+m[2]
    force=M12*R12**3+m[5]
    G=M12*R12**2*(m[1]*R25**2+m[2])-m[1]*m[2]*R25**2
    H=M12*(m[1]+m[2]*R25)*R12**2-m[1]*m[2]
    Phi=m[5]**2*M12*(m[1]+m[2]*R25)**3-(m[1]*m[2])**3
    equal(H,R25*G+(1-R25**3)*m[1]*(M12*R12**2-m[2]),'D8 H syzygy')
    equal(Phi,H*(H**2+3*m[1]*m[2]*H+3*m[1]**2*m[2]**2)-
          M12*(m[1]+m[2]*R25)**3*force*(force-2*m[5]),'D8 Phi syzygy')
    omega=V('d8_phase')
    remainder=substitute(Phi,{'d8_R25':omega}).reduce('d8_phase',[1,1,1])
    equal(remainder.coefficient('d8_phase',1),3*m[1]*m[2]*m[5]**2*M12*(m[1]-m[2]),'D8 phase coefficient')
    equal(substitute(remainder.coefficient('d8_phase',0),{'m2':m[1]}),-m[1]**4*(m[1]**2+2*m[5]**2),'D8 positive-mass phase exclusion')
    return [dict(check='D8.first_condition.mass_elimination_syzygies',status='PASS',items=5,
                 detail=dict(scope='P,Q elimination, H/Phi identities, and exact phase exclusion'))],{}


if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram8_relation40_derivation'))
