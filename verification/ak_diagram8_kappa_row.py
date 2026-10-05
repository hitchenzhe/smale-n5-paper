"""Table 5: exact beta, Jacobian and tau computations of Lemma 11.7."""
from fractions import Fraction as F
from itertools import combinations, product
from smale_verify.exact import P,V,zero,equal,substitute,diff_poly,det
from ak_diagram8_relation40_derivation import chart_rows

UPPER_ROWS=(0,2,10,13,16,20)
VARIABLES=('d8_A1','d8_A2','d8_R12','d8_R25','d8_P','d8_Q')

def clear_substitute(poly,mapping):
    for name in mapping: poly=poly.clear(name)
    return substitute(poly,mapping)

def rank_pivots(matrix):
    a=[list(row) for row in matrix]
    row_indices=list(range(len(a)))
    r=0; columns=[]; rows=[]
    for col in range(len(a[0])):
        pivot=next((j for j in range(r,len(a)) if a[j][col]),None)
        if pivot is None: continue
        a[r],a[pivot]=a[pivot],a[r]
        row_indices[r],row_indices[pivot]=row_indices[pivot],row_indices[r]
        rows.append(row_indices[r]); columns.append(col)
        lead=a[r][col];a[r]=[x/lead for x in a[r]]
        for j in range(r+1,len(a)):
            q=a[j][col];a[j]=[x-q*y for x,y in zip(a[j],a[r])]
        r+=1
        if r==len(a):break
    return r,rows,columns

def verify_beta_identity():
    a,b,c,s,h=map(V,('b8_m1','b8_m2','b8_m5','b8_sum','b8_total'))
    m={1:a,2:b,3:V('b8_m3'),4:V('b8_m4'),5:c}
    data=chart_rows(m)
    diff=a**2*b**2*c**-3*s**-2
    pp=(-s*b**-1-b*diff)*h**-1
    point={'d8_A1':-b/s,'d8_A2':a/s,'d8_C':0,
           'd8_R12':-c*s/(a*b),'d8_R25':1,'d8_P':pp,'d8_Q':pp+diff}
    actual=substitute(data['beta1']-data['beta3'],point)
    expected=(a**3*b**3-c**2*s**4)/(a*b*c*s**2)
    zero(clear_substitute(actual-expected,{'b8_sum':a+b,'b8_total':a+b+c}),'D8 beta at I0 exact mass identity')
    return [dict(check='D8.11.7.beta_mass_identity',status='PASS',items=1,
                 detail=dict(scope='exact identity before imposing the mass relation'))],{}

def general_upper_data():
    alpha,s,g=map(V,('al8','s8','g8'))
    m={1:alpha**2,2:P.coerce(1),3:V('b3'),4:V('b4'),5:alpha**3*s**-2}
    data=chart_rows(m)
    pp=(-s-s**4*alpha**-5)*s**2/g
    point={'d8_A1':-s**-1,'d8_A2':alpha**2/s,'d8_C':0,
           'd8_R12':-alpha/s,'d8_R25':1,'d8_P':pp,'d8_Q':pp+s**4*alpha**-5}
    J=[[substitute(diff_poly(data['f'][j],v),point) for v in VARIABLES] for j in UPPER_ROWS]
    return data,m,point,alpha,s,g,J

def verify_general_jacobian():
    data,m,point,alpha,s,g,J=general_upper_data()
    determinant=det(J)
    positive_quartic=alpha**4-alpha**3+3*alpha**2-alpha+1
    expected=3*alpha**10*(alpha**2-1)*(alpha**2+alpha+1)*positive_quartic*s**-10
    s_poly=alpha**2+1
    mapping={'s8':s_poly,'g8':s_poly**3+alpha**3}
    zero(clear_substitute(determinant-expected,mapping),'D8 general 6x6 Jacobian determinant')
    equal(positive_quartic,(alpha**2-alpha/2+1)**2+F(3,4)*alpha**2,'D8 quartic positivity witness')
    zero(clear_substitute(substitute(data['beta1']-data['beta3'],point),mapping),'D8 beta on positive mass parametrization')
    return [dict(check='D8.11.7.general_jacobian_determinant',status='PASS',items=3,
                 detail=dict(determinant='3 alpha^10 (alpha^2-1) (alpha^2+alpha+1) (alpha^4-alpha^3+3alpha^2-alpha+1)/(alpha^2+1)^10',
                             positive_mass_singularities=[1],
                             quartic_positivity='(alpha^2-alpha/2+1)^2+3 alpha^2/4'))],{}

def verify_equal_pair_kernel():
    m={1:P.coerce(1),2:P.coerce(1),3:V('m3'),4:V('m4'),5:P.coerce(F(1,4))}
    data=chart_rows(m)
    point=dict(zip(('d8_A1','d8_A2','d8_C','d8_R12','d8_R25','d8_P','d8_Q'),map(F,('-1/2','1/2','0','-1/2','1','-8','8'))))
    k=list(map(F,('27/272','27/272','9/272','-9/68','-19/17','1')))
    rows=[data['f'][j] for j in UPPER_ROWS]
    J=[[substitute(diff_poly(row,v),point) for v in VARIABLES] for row in rows]
    numeric=[[entry.terms.get((),F(0)) for entry in row] for row in J]
    if any(entry.names() for row in J for entry in row): raise AssertionError('D8 equal B is not rational')
    rank,ri,ci=rank_pivots(numeric)
    if rank!=5:raise AssertionError(('D8 equal upper rank',rank))
    minor=det([[J[r][c] for c in ci] for r in ri])
    if not minor.terms:raise AssertionError('zero D8 rank-five witness')
    for n,(row,jacrow) in enumerate(zip(rows,J)):
        zero(substitute(row,point),'D8 equal face '+str(n))
        zero(sum((entry*q for entry,q in zip(jacrow,k)),P()),'D8 equal B kernel row '+str(n))
    beta=data['beta1']-data['beta3']
    zero(substitute(beta,point),'D8 equal beta value')
    zero(sum((substitute(diff_poly(beta,v),point)*q for v,q in zip(VARIABLES,k)),P()),'D8 beta gradient on kernel')
    U=('d8_A1','d8_A2','d8_C','d8_R12','d8_R25','d8_P','d8_Q')
    center=data['x']['A1']+data['x']['A2']+(m[3]+m[4]+m[5])*data['x']['C']
    augmented=rows+[data['x']['C'],center]
    AJ=[[substitute(diff_poly(row,v),point) for v in U] for row in augmented]
    witness=None
    for ids in combinations(range(8),7):
        if 6 not in ids or 7 not in ids:continue
        candidate=det([AJ[i] for i in ids])
        if candidate.terms and not candidate.names():witness=(ids,candidate.terms[()]);break
    if witness is None:raise AssertionError('no mass-uniform D8 rank-seven augmented minor')
    return [dict(check='D8.11.7.equal_kernel_and_augmented_rank',status='PASS',items=15,
                 detail=dict(kernel=list(map(str,k)),B_matrix=[[str(q) for q in row] for row in numeric],
                             rank=5,minor_rows=ri,minor_columns=ci,minor_value=str(minor.terms[()]),
                             augmented_minor_rows=list(witness[0]),augmented_minor_value=str(witness[1]),
                             augmented_rank=7,lower_masses='free; chosen minor is constant'))],{}

def verify_tau_sign_identity():
    mus={j:V('mu8_'+str(j)) for j in range(1,5)}
    m={j:mus[j]**2 for j in mus};m[5]=V('m5')
    data=chart_rows(m)
    ii,L,h=V('imag8'),V('L8'),V('M34inv8')
    for signs in product((1,-1),repeat=4):
        s13,s14,s23,s24=signs
        point={'d8_D':0,'d8_B3':-m[4]*h,'d8_B4':m[3]*h,
               'd8_D13':ii*s13*L/(mus[2]*mus[4]),'d8_D14':s14*L/(mus[2]*mus[3]),
               'd8_D23':s23*L/(mus[1]*mus[4]),'d8_D24':ii*s24*L/(mus[1]*mus[3])}
        actual=substitute(data['tau1']-data['tau3'],point).reduce('imag8',[1,0,1])
        Xi=s14*(mus[1]*mus[4])**3+s23*(mus[2]*mus[3])**3+ii*(s13*(mus[1]*mus[3])**3+s24*(mus[2]*mus[4])**3)
        expected=m[3]*m[4]*h*L**3*Xi/(mus[1]*mus[2]*mus[3]*mus[4])**3
        equal(actual,expected,'D8 tau-Xi identity '+str(signs))
    return [dict(check='D8.11.7.tau_Xi_all_signs',status='PASS',items=16,
                 detail=dict(sign_cases=16,field='rational Laurent polynomials in positive mass square roots and i^2=-1'))],{}

if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram8_kappa_row'))
