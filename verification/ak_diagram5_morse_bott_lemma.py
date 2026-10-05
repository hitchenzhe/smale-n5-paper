"""Table 5: exact triangle Hessian, Lemma 8.16."""
from fractions import Fraction as F
from smale_verify.exact import P,V,zero,equal,divide_poly,det


def verify_triangle_hessian():
    b,r,a,rho=map(V,('d5mb_b','d5mb_r','d5mb_a','d5mb_rho'))
    ap=rho**2/a
    mass=[b**2,b**2/r**2,r**-2]
    offset=[P(),r,-r/b]
    inverse={(0,1):r**-1,(0,2):b/r,(1,2):-b}
    relation=r*(b+1)-1
    def check(value,label):
        value=value.clear('d5mb_b').clear('d5mb_r').clear('d5mb_a').clear('d5mb_rho')
        zero(divide_poly(value,[relation],('d5mb_r','d5mb_b','d5mb_a','d5mb_rho')),label)
    U=sum((mass[i]*mass[j]*D/rho for (i,j),D in inverse.items()),P())
    check(U,'D5 U vanishes on the triangle orbit')
    force=[P(),P(),P()]
    L=[[P() for _ in range(3)] for _ in range(3)]
    for (i,j),D in inverse.items():
        check(D**2*(offset[j]-offset[i])**2-1,'D5 base distance tie')
        force[i]+=mass[j]*D**3*(offset[i]-offset[j])
        force[j]+=mass[i]*D**3*(offset[j]-offset[i])
        weight=mass[i]*mass[j]*D**3
        L[i][i]+=weight;L[j][j]+=weight;L[i][j]-=weight;L[j][i]-=weight
    for value in force:check(value,'D5 orbit critical point')
    # Differentiate each signed inverse-distance summand using D^2 dx dy=1.
    # Its xx,xy,yy derivatives are 3D/(4dx^2),D/(4dxdy),3D/(4dy^2).
    H=[[P() for _ in range(6)] for _ in range(6)]
    for i in range(3):
        for j in range(3):
            H[i][j]=3*L[i][j]/(4*rho*a**2)
            H[i+3][j+3]=3*L[i][j]/(4*rho*ap**2)
            H[i][j+3]=H[j+3][i]=L[i][j]/(4*rho*a*ap)
    for vector in ([1,1,1,0,0,0],[0,0,0,1,1,1],
                   offset+[0,0,0],[0,0,0]+offset):
        for row in H:check(sum((x*y for x,y in zip(row,vector)),P()),'D5 Hessian tangent kernel')
    X=b**4/(4*r**6)
    block=[[H[0][0],H[0][3]],[H[3][0],H[3][3]]]
    check(block[0][0]-3*X/(rho*a**2),'D5 complement xx')
    check(block[0][1]-X/(rho*a*ap),'D5 complement xy')
    check(block[1][1]-3*X/(rho*ap**2),'D5 complement yy')
    check(det(block)-8*X**2/rho**6,'D5 complement determinant')
    inverse_xx=3*a**2*rho/(8*X)
    check(det(block)*inverse_xx-block[1][1],'D5 inverse xx')
    # The two directions are a genuine complement: [1,s,e3] has nonzero
    # determinant s5-s4=-1/b, including the central-body slice b=1.
    check(offset[2]-offset[1]+b**-1,'D5 tangent/complement determinant')
    return [dict(check='D5.8.16.orbit_Hessian',status='PASS',items=39,
        detail=dict(masses=['b^2','b^2(b+1)^2','(b+1)^2'],
            X='b^4(b+1)^6/4',rho_square='a a_prime',
            complement_H=['3X/(rho a^2)','X/(rho a a_prime)','3X/(rho a_prime^2)'],
            determinant='8X^2/rho^6',inverse_xx='3a^2rho/(8X)',
            Hessian_rank=2,kernel_dimension=4,
            nonvanishing='b>0 and nonzero a,a_prime,rho; complement determinant -1/b'))],{}


if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram5_morse_bott_lemma'))

