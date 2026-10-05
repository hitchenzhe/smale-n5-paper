"""Table 5: exact D2 shape resultant and class assignments (Lemma 6.2).

Mass resultants L_j and P_j remain a distinct unfinished task in the manifest.
"""
from fractions import Fraction as F

def verify_mass_resultants():
    from smale_verify.d2_mass_resultants import verify_mass_resultants as replay
    return replay()
from itertools import product
from smale_verify.exact import P,V,zero,equal,substitute,det

def shape_factors():
    w3,w4,w5=map(V,('shape_w3','shape_w4','shape_w5'))
    U,VV=w4-w3,w5-w4
    def numerator(signs):
        s3,s4,s5=signs
        return s3*w3*(w3-w5)**2*(w4-w3)**2+s4*w4*(w5-w4)**2*(w4-w3)**2+s5*w5*(w5-w4)**2*(w3-w5)**2
    S={j:numerator(signs) for j,signs in enumerate(((1,-1,-1),(1,-1,1),(1,1,-1)),1)}
    S[4]=w3**3+w4**3+w5**3-w3**2*w4-w3*w4**2-w3**2*w5-w3*w5**2-w4**2*w5-w4*w5**2+3*w3*w4*w5
    S[5]=U**2+U*VV+VV**2
    S[6]=U**4+2*U**3*VV-U*VV**3+VV**4
    S[7]=U**4+5*U**3*VV+9*U**2*VV**2+5*U*VV**3+VV**4
    S[8]=U**4-U**3*VV+2*U*VV**3+VV**4
    return S,U,VV

def verify_shape_resultant():
    S,U,VV=shape_factors()
    # Sylvester determinant of the two cubics in the definition, in p:
    # U^3-p^3 V^3, U^3-(1+p)^3(U+V)^3.
    first=[-VV**3,P(),P(),U**3]
    cube=(U+VV)**3
    second=[-cube,-3*cube,-3*cube,U**3-cube]
    matrix=[]
    for coeffs in (first,second):
        for shift in range(3):matrix.append([P()]*shift+coeffs+[P()]*(2-shift))
    actual=det(matrix)
    equal(actual,S[5]**3*S[6]*S[7]*S[8],'D2 complete shape resultant factorization')
    for j,expected in ((6,-2),(7,1),(8,-5)):
        coeff=S[j].coefficient('shape_w3',3).coefficient('shape_w4',1).coefficient('shape_w5',0)
        equal(coeff,expected,'D2 S'+str(j)+' identifying coefficient')
    equal(S[5],V('shape_w3')**2+V('shape_w4')**2+V('shape_w5')**2-V('shape_w3')*V('shape_w4')-V('shape_w3')*V('shape_w5')-V('shape_w4')*V('shape_w5'),'D2 S5 translation form')
    return [dict(check='D2.6.2.shape_resultant_factorization',status='PASS',items=5,
                 detail=dict(resultant='S5^3 S6 S7 S8 (exact constant 1)',
                             quartic_labels={'S6_w3^3w4':-2,'S7_w3^3w4':1,'S8_w3^3w4':-5}))],{}

def verify_lagrange_assignments():
    shapes,_,_=shape_factors()
    k,omega=V('shape_kappa'),V('shape_omega')
    assignment_counts={5:0,6:0,7:0,8:0};records=[]
    for a34,a35,a45 in product(range(3),repeat=3):
        c=(omega**((a34-a45)%3)).reduce('shape_omega',[1,1,1])
        b=(omega**((a34-a45)%3)+1-omega**((a35-a45)%3)).reduce('shape_omega',[1,1,1])
        folded=len({a34,a35,a45})==3
        j=5 if len({a34,a35,a45}) in (1,3) else 6 if a34==a35 else 7 if a34==a45 else 8
        shape=substitute(shapes[j],{'shape_w3':0,'shape_w4':k,'shape_w5':k+1})
        if folded:
            discriminant=(b**2-4*c).reduce('shape_omega',[1,1,1])
            zero(discriminant,'D2 fold discriminant '+str((a34,a35,a45)))
            reduced=substitute(shape,{'shape_kappa':-b/2}).reduce('shape_omega',[1,1,1])
        else:
            reduced=shape.reduce('shape_kappa',[c,b,1]).reduce('shape_omega',[1,1,1])
        zero(reduced,'D2 full S'+str(j)+' class assignment '+str((a34,a35,a45)))
        assignment_counts[j]+=1
        records.append(dict(exponents=[a34,a35,a45],factor='S'+str(j),fold=folded))
    if assignment_counts!={5:9,6:6,7:6,8:6}:raise AssertionError(('D2 shape assignments',assignment_counts))
    return [dict(check='D2.6.2.all_Lagrange_shape_assignments',status='PASS',items=27,
                 detail=dict(class_assignments=records,counts=assignment_counts,
                             field='Q(omega), omega^2+omega+1=0'))],{}

def verify_euler_on_S5(m):
    k=V('euler_kappa')
    cases=[]
    for e34,e35,e45 in product((1,-1),repeat=3):
        # Difference of the two Euler force equations after eliminating c^3.
        # On S5, k^3=1 and (k+1)^3=-1; inverse squares are k and k^2.
        A=(m[3]+m[4])*e34*k+m[5]*e35*k**2-m[5]*e45
        B=(m[4]+m[5])*e45+m[3]*e35*k**2-m[3]*e34*k
        eliminant=(k*B-A).reduce('euler_kappa',[1,1,1])
        constant=m[3]*(e35+e34)+m[5]*(e35+e45)
        coefficient=m[4]*(e45-e34)+m[5]*(e45+e35)
        equal(eliminant,constant+coefficient*k,'D2 Euler S5 force eliminant '+str((e34,e35,e45)))
        allowed=e34==e45==-e35
        if allowed:
            zero(eliminant,'D2 Euler S5 surviving signs')
            zero((e34-e45*k**3).reduce('euler_kappa',[1,1,1]),'D2 Euler D34^3=D45^3')
            zero((e35-e45*(k+1)**3).reduce('euler_kappa',[1,1,1]),'D2 Euler D35^3=D45^3')
        else:
            values=list(constant.terms.values())
            if not values or not (all(v>0 for v in values) or all(v<0 for v in values)):
                raise AssertionError('D2 Euler S5 missing strict positive-mass exclusion')
        cases.append(dict(signs=[e34,e35,e45],allowed=allowed,
                          argument='cubes coincide' if allowed else 'nonzero real constant coefficient and nonreal kappa'))
    return [dict(check='D2.6.2.Euler_on_S5_cube_equality',status='PASS',items=8,
                 detail=dict(cases=cases,surviving_signs='e34=e45=-e35',
                             scope='all Euler signs on S5 at positive masses'))],{}

if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_resultants'))
