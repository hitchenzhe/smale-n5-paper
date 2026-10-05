"""Global regular S6/S8 exclusion by a weighted-Laplacian cofactor.

All root and mass enclosures use rational outward interval arithmetic.
The four small boxes exhaust the possible cofactor zeros, not a sampling
of the mass cone. No floating arithmetic participates in acceptance.
"""
from fractions import Fraction as Q
from math import isqrt
import json
from pathlib import Path
from .interval import R,Z,GRID,certify_power_root
from .symbolic import sympy

CENTRES={6:[('-1.4735614833535065','-.4447718087620662'),('.4735614833535065','.4447718087620662')],
         8:[('-.6219644269518568','.1877303704569451'),('1.1219644269518568','-1.0537557742413837')]}

def square_root(x):
    if x.lo<0:raise ArithmeticError('Negative square-root interval')
    lower=isqrt((x.lo*GRID**2).__floor__());upper=isqrt((x.hi*GRID**2).__ceil__())+1
    return R(Q(lower,GRID),Q(upper,GRID))

def verify():
    s=sympy();a,b,c,h12,h13,h23=s.symbols('a b c h12 h13 h23')
    masses=[a,b,c];weights={(0,1):h12,(0,2):h13,(1,2):h23}
    C=s.zeros(3)
    for (i,j),h in weights.items():
        C[i,i]+=masses[j]*h;C[i,j]-=masses[j]*h
        C[j,j]+=masses[i]*h;C[j,i]-=masses[i]*h
    tau=a*h12*h13+b*h12*h23+c*h13*h23
    if s.expand(C[:2,:2].det()-c*tau):raise ArithmeticError('Laplacian cofactor identity')
    if C*s.ones(3,1)!=s.zeros(3,1):raise ArithmeticError('Translation kernel identity')
    root3,root_record=certify_power_root(Z(3),2,('1.732050807568877293527446341506','0'),'1e-27')
    omega=(-1+Z(0,1)*root3)/2;records=[];basis_checks=[]
    kap=s.symbols('kap');om=(-1+s.I*s.sqrt(3))/2
    for sheet,coefficient,constant in ((6,s.Integer(1),om**2),(8,om,om)):
        determinant=s.resultant(kap**2+coefficient*kap+constant,kap**2-constant,kap)
        if s.simplify(determinant)==0:raise ArithmeticError('Nonfold 1,W,Z basis determinant')
        basis_checks.append(dict(sheet=sheet,resultant=str(s.simplify(determinant))))
    for j,phases in ((6,(Z(1),Z(1),omega)),(8,(omega,Z(1),Z(1)))):
        A=Z(1);B=Z(1)+phases[2]**2/phases[0]**2-phases[2]**2/phases[1]**2
        D=phases[2]**2/phases[0]**2;disks=[]
        for centre in CENTRES[j]:
            k0=Z(*centre);radius=Q('1e-12')
            dominant=(2*k0+B).modulus_lower()*radius
            remainder=(k0*k0+B*k0+D).upper()+radius**2
            if dominant<=remainder:raise ArithmeticError('Quadratic Rouche inequality')
            k=Z(R(Q(centre[0])-radius,Q(centre[0])+radius),R(Q(centre[1])-radius,Q(centre[1])+radius))
            if not k.im.lower():raise ArithmeticError('Shape is not separated from real axis')
            h=[phases[0]**2*k**2,phases[1]**2*(k+1)**2,phases[2]**2]
            u,v,w=h[0]*h[1],h[0]*h[2],h[1]*h[2]
            den=u.re*v.im-u.im*v.re
            aa=(-w.re*v.im+w.im*v.re)/den
            bb=(-u.re*w.im+u.im*w.re)/den
            if aa.lo<=0 or bb.lo<=0:raise ArithmeticError('Mass box not positive')
            W4=(k*aa-1)/(Z(aa+bb)+k*aa);W3=-(W4*bb+1)/aa
            amplitudes=[aa/square_root(square_root(W3.re.square()+W3.im.square())),bb/square_root(square_root(W4.re.square()+W4.im.square())),R(1)]
            winner=max(range(3),key=lambda i:amplitudes[i].lo)
            margin=amplitudes[winner].lo-sum(v.hi for i,v in enumerate(amplitudes)if i!=winner)
            if margin<=0:raise ArithmeticError('No certified strict amplitude dominance')
            records.append(dict(sheet=j,root_disk=dict(centre=list(centre),radius=str(radius),linear_lower=str(dominant),remainder_upper=str(remainder)),
                real_linear_determinant=den.encode(),mass_box=[aa.encode(),bb.encode()],W_box=[W3.encode(),W4.encode()],
                amplitudes=[v.encode()for v in amplitudes],dominant_body=winner+3,strict_margin=str(margin)))
            disks.append(k)
        if (disks[0].midpoint()-disks[1].midpoint()).modulus_lower()<=Q('2e-12'):raise ArithmeticError('Root disks overlap')
    payload=dict(check='D2.6.17.S6_S8_global_cofactor_exclusion',status='PASS',omega_root=root_record,cases=records,nonfold_basis_resultants=basis_checks,
        cofactor_identity='det C[3:5,3:5]=m5*(m3*h34*h35+m4*h34*h45+m5*h35*h45)',
        proof=['D^3=1/MT. For a nonfold Lagrange shape, the vectors 1,W,Z are independent.',
            'The homogeneous triangle differential forces its edge-Laplacian variation to annihilate 1,W,Z, hence all three D derivatives vanish. The distance equations leave only rotation; W5=1 kills it. Thus the selected nine-row triangle Jacobian is invertible.',
            'Centring and gauge give eta3=-m4/m3*eta4 and eta5=0. A constant forcing translates Z and does not change W. Therefore Psi1=-m4/2*g*(xi4/W4-xi3/W3)^2.',
            'If g=0 then eta=0. Eliminating D derivatives in the W equations gives C*zeta=0, with weights proportional to c_edge^2*(Wdifference)^2. Zeta is not constant, since a constant Zeta would force Dprime=0 and could not solve the unit FZ4 forcing. Hence rank C<=1 and tau=0.',
            'The quadratic shape equations have exactly the two certified disjoint roots each. In each root disk tau=0 is a nonsingular pair of real linear equations for (a,b); the displayed interval solution exhausts all real mass solutions.',
            'At each possible tau zero one amplitude mB/sqrt(abs(WB)) strictly exceeds the sum of the other two. Thus Phi=0 is impossible for ANY unit phases, including every choice of rho and signs.',
            'X=xi4/W4-xi3/W3 is excluded separately by the exact reciprocal-polynomial gcd certificate. Body exchanges, conjugation and common cubic scales cover all S6/S7/S8 labels.'],
        depends_on=['D2.6.17.all_regular_X_unit_circle_exclusions','D2.6.2.all_Lagrange_shape_assignments'])
    path=Path(__file__).resolve().parents[1]/'certificates/d2_S6_S8_global.json';path.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=4,detail=dict(certificate=path.name,complete_positive_mass_domain=True))],{}
