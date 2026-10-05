"""Table 5: exact lower determinant and second-order constants, Lemma 11.8."""
from fractions import Fraction as F
from smale_verify.exact import P,V,zero,equal,substitute,diff_poly,det,solve_rational
from ak_diagram8_relation40_derivation import chart_rows
from ak_diagram8_kappa_row import VARIABLES,rank_pivots

def equal_data():
    mass=V('cancel_m')
    m={1:P.coerce(1),2:P.coerce(1),3:mass,4:mass,5:P.coerce(F(1,4))}
    data=chart_rows(m)
    point=dict(zip(('d8_A1','d8_A2','d8_C','d8_R12','d8_R25','d8_P','d8_Q'),map(F,('-1/2','1/2','0','-1/2','1','-8','8'))))
    k=list(map(F,('27/272','27/272','9/272','-9/68','-19/17','1')))
    return data,m,point,k

def scalar(poly):
    if poly.names():raise AssertionError(('not a rational scalar',poly.names()))
    return poly.terms.get((),F(0))

def verify_forcing_support():
    data,m,point,k=equal_data();E,f,x=data['E'],data['f'],data['x']
    rows=(1,3,5,7,9,17,18,19,21,11,12,14,15)
    exempt=[]
    for j in rows:
        for mon in (E[j]-f[j]).terms:
            ex=dict(mon);i=ex.get('a8',0);q=ex.get('eta8',0)
            if (i,q)==(0,2) and j in (1,3,9):continue
            if i<2 or i+2*q<6:raise AssertionError(('D8 lower forcing weight',j,mon))
        if j in (1,3,9):
            coefficient=E[j].coefficient('a8',0).coefficient('eta8',2)
            zero(substitute(coefficient,point),'D8 exceptional eta2 coefficient '+str(j))
            exempt.append(j)
    correction=data['nu']*(m[3]*x['u']+m[4]*x['v'])
    mapping={'d8_C':x['C']-correction}
    for j in (0,2,8,10,13,16,20):
        amplitude=x['A1'] if j==0 else x['A2'] if j==2 else x['C'] if j==8 else P()
        remainder=E[j]-substitute(f[j],mapping)-data['a']**6*amplitude
        if any(dict(mon).get('a8',0)<8 for mon in remainder.terms):
            raise AssertionError(('D8 upper row below a8',j))
    for j,beta,tau,pname in ((1,data['beta1'],data['tau1'],'P'),(3,data['beta3'],data['tau3'],'Q')):
        remainder=E[j]-x['D']-data['eta']**2*beta-(data['a']*data['eta'])**3*tau-m[5]*data['a']**6*data['eta']**2*x[pname]
        found={tuple(dict(mon).get(v,0) for v in ('a8','eta8')) for mon in remainder.terms}
        if found!={(9,5)}:raise AssertionError(('D8 exact upper w remainder',j,found))
    return [dict(check='D8.11.8.forcing_and_upper_support',status='PASS',items=25,
                 detail=dict(lower_rows=list(rows),eta2_exemptions=exempt,upper_rows=7,
                             w_remainder_exponents=[9,5]))],{}

def verify_second_order_constants():
    data,m,point,k=equal_data()
    f=data['f'];beta=data['beta1']-data['beta3']
    rows=(2,10,13,16,20)
    def grad(row):return [scalar(substitute(diff_poly(row,v),point)) for v in VARIABLES]
    def quadratic(row):
        return sum((k[i]*k[j]*scalar(substitute(diff_poly(diff_poly(row,v),w),point))
                    for i,v in enumerate(VARIABLES) for j,w in enumerate(VARIABLES)),F(0))/2
    B5=[grad(f[j])[:5] for j in rows]
    determinant=scalar(det(B5))
    if determinant!=F(-17,16):raise AssertionError(('D8 B5 determinant',determinant))
    delta_mu=solve_rational(B5,[F(-1,2),0,0,0,0])+[F(0)]
    delta_2=solve_rational(B5,[-quadratic(f[j]) for j in rows])+[F(0)]
    def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
    constants={'c0':F(-1,2)+dot(grad(f[0]),delta_mu),
               'q2':quadratic(f[0])+dot(grad(f[0]),delta_2),
               'e2':quadratic(beta)+dot(grad(beta),delta_2),
               'c62':dot(grad(beta),delta_mu)+F(1,4)*(-8-8)}
    constants['det_pairing']=constants['c0']*constants['e2']-constants['q2']*constants['c62']
    expected=dict(c0=F(-1),q2=F(243,73984),e2=F(1215,18496),c62=F(-4),det_pairing=F(-243,4624))
    if constants!=expected:raise AssertionError(('D8 second-order constants',constants,expected))
    # Both solved systems are independently replayed, including all five rows.
    for r,row in enumerate(B5):
        if dot(row,delta_mu[:5])!=([F(-1,2),0,0,0,0][r]):raise AssertionError('invalid delta_mu')
        if dot(row,delta_2[:5])!=-quadratic(f[rows[r]]):raise AssertionError('invalid delta_2')
    upper=[grad(f[j]) for j in (0,2,8,10,13,16,20)]
    if rank_pivots(upper)[0]!=5:raise AssertionError('D8 seven upper rows do not have rank five')
    for row in upper:
        if dot(row,k):raise AssertionError('D8 cancelling-curve kernel mismatch')
    if k[0]+k[1]!=F(27,136):raise AssertionError('D8 varphi(k0) mismatch')
    return [dict(check='D8.11.8.second_order_constants',status='PASS',items=19,
                 detail=dict(B5=[[str(q) for q in row] for row in B5],B5_det=str(determinant),
                             delta_mu=list(map(str,delta_mu)),delta_2=list(map(str,delta_2)),
                             constants={key:str(value) for key,value in constants.items()},
                             kernel=list(map(str,k)),varphi_kernel='27/136'))],{}

def verify_lower_block():
    data,m,point,k=equal_data();f=data['f']
    mass=m[3];g,lam,omega,S=map(V,('cancel_g','cancel_lambda','cancel_omega','cancel_S'))
    u=-32*g**-2*S**-1*(mass+F(1,4)+mass*omega)
    v=32*g**-2*S**-1*(mass+(mass+F(1,4))*omega)
    lower_point={'d8_B3':F(-1,2),'d8_B4':F(1,2),'d8_D':0,'d8_u':u,'d8_v':v,
                 'd8_D34':lam*g,'d8_D35':g,'d8_D45':omega*g}
    indices=(1,5,7,9,17,18,19,21)
    variables=('d8_B3','d8_B4','d8_D','d8_u','d8_v','d8_D34','d8_D35','d8_D45')
    def reduce_label(poly):
        poly=poly.clear('cancel_S')
        poly=substitute(poly,{'cancel_S':8*mass+1})
        poly=substitute(poly,{'cancel_m':(4*g**-3-1)/(8*lam**3)})
        for name in ('cancel_g','cancel_lambda','cancel_omega'):poly=poly.clear(name)
        poly=poly.reduce('cancel_lambda',[-(1-omega+omega**2)/4,0,1])
        return poly.reduce('cancel_omega',[-1,0,0,1])
    for j in indices:zero(reduce_label(substitute(f[j],lower_point)),'D8 lower atlas face '+str(j))
    J=[[substitute(diff_poly(f[j],vname),lower_point) for vname in variables] for j in indices]
    actual=det(J)
    expected=-9*g**7*lam*mass*(8*mass+1)*(omega**2+1)*(8*mass*lam**3+1)/(128*omega)
    zero(reduce_label(actual-expected),'D8 lower 8x8 Jacobian determinant')
    equal(((omega**2+1)*(1+omega-omega**2)/2).reduce('cancel_omega',[-1,0,0,1]),1,'D8 omega2+1 invertible')
    return [dict(check='D8.11.8.lower_atlas_and_jacobian',status='PASS',items=10,
                 detail=dict(face_rows=list(indices),variables=list(variables),
                             determinant='-9 g^7 lambda m (8m+1)(omega^2+1)(8m lambda^3+1)/(128 omega)',
                             reduction_relations=['omega^3=1','lambda^2=(1-omega+omega^2)/4','m=(4g^-3-1)/(8lambda^3)'],
                             nonvanishing='g,lambda,omega nonzero; m>0; 8m+1>0; scale factor=4/g^3'))],{}

if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram8_cancelling_curve'))
