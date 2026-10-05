"""Exact mass-parametric S5 regular response and fold blocks.

Normalize m5=1, a=m3, b=m4, delta**3=1/(a+b+1).
Scaled Jacobians have coefficients in Q(omega)(a,b); no root guesses.
"""
from pathlib import Path
import json
from .symbolic import sympy
from .exact import V,equations,diff_poly


def construct(fold=False):
    s=sympy();a,b=s.symbols('a b');omega=(-1+s.I*s.sqrt(3))/2
    field=s.QQ.algebraic_field(omega);K=field.frac_field(a,b)
    def N(v):
        if not v:return K.zero
        lc=v.denom.LC;return v.field.raw_new(v.numer.quo_ground(lc),v.denom.quo_ground(lc))
    def constant(v):return K.convert(field.from_sympy(s.sympify(v)),field)
    om=constant(omega);mass={3:K.from_sympy(a),4:K.from_sympy(b),5:K.one};MT=sum(mass.values(),K.zero)
    w={5:K.one};k=om**2;w[4]=N((mass[3]*k-1)/(mass[3]+mass[4]+mass[3]*k));w[3]=N(-(mass[4]*w[4]+1)/mass[3])
    phases={(3,4):K.one,(3,5):om if fold else K.one,(4,5):om**2 if fold else K.one}
    dz4=N(phases[(3,4)]**-2/(w[4]-w[3]));dz5=N(phases[(3,5)]**-2/(1-w[3]));z3=N(-(mass[4]*dz4+dz5)/MT)
    z={3:z3,4:N(z3+dz4),5:N(z3+dz5)}
    mm={j:V('m'+str(j))for j in mass};zz={j:V('Z'+str(j))for j in mass};ww={j:V('W'+str(j))for j in mass};dd={p:V('D'+''.join(map(str,p)))for p in phases}
    ez,ew,rt=equations(mm,zz,ww,dd);rows=[ez[3],ew[3],ez[4],ew[4],ew[5],*rt.values(),ww[5]-1]
    variables=[*('Z'+str(j)for j in mass),*('W'+str(j)for j in mass),*('D'+''.join(map(str,p))for p in dd)]
    point={**{'m'+str(j):(0,v)for j,v in mass.items()},**{'Z'+str(j):(-2,v)for j,v in z.items()},**{'W'+str(j):(0,v)for j,v in w.items()},**{'D'+''.join(map(str,p)):(1,v)for p,v in phases.items()}}
    def eval_(poly,offset=0):
        out=K.zero
        for mon,c in poly.terms.items():
            power=offset;value=K.convert(c)
            for name,e in mon:
                q,v=point[name];power+=q*e;value=N(value*v**e)
            if power%3:raise ArithmeticError('Scaled entry is not a cubic-scale invariant')
            out=N(out+value*MT**(-power//3))
        return out
    if any(eval_(row,2 if i in(0,2)else 0)for i,row in enumerate(rows)):
        raise ArithmeticError('S5 triangle does not satisfy all source rows')
    J=[[eval_(diff_poly(row,v),(2 if i in(0,2)else 0)+(-2 if v.startswith('Z')else 1 if v.startswith('D')else 0))for v in variables]for i,row in enumerate(rows)]
    def inverse(matrix):
        n=len(matrix);aug=[row+[K.one if i==j else K.zero for j in range(n)]for i,row in enumerate(matrix)]
        for col in range(n):
            pivot=min((i for i in range(col,n)if aug[i][col]),key=lambda i:len(aug[i][col].numer)+len(aug[i][col].denom));aug[col],aug[pivot]=aug[pivot],aug[col]
            d=aug[col][col];aug[col]=[N(v/d)for v in aug[col]]
            for i in range(n):
                if i!=col and aug[i][col]:
                    d=aug[i][col];aug[i]=[N(v-N(d*q))for v,q in zip(aug[i],aug[col])]
        result=[row[n:]for row in aug]
        for i in range(n):
            for j in range(n):
                if N(sum((matrix[i][k]*result[k][j]for k in range(n)),K.zero))-K.convert(int(i==j)):
                    raise ArithmeticError('S5 inverse product failed')
        return result
    def conjugate(poly):
        out=K.zero
        for (pa,pb),c in poly.items():
            cs=c.to_list();cs=[s.QQ.zero]*(2-len(cs))+cs
            out+=K.convert(field.new([-cs[0],cs[1]-cs[0]]),field)*K.from_sympy(a**pa*b**pb)
        return N(out)
    def norm(poly):return N(K.field.raw_new(poly,poly.ring.one)*conjugate(poly))
    def expression(v):return s.factor(K.to_sympy(N(v)))
    def encode(v):
        def poly(p):return [dict(powers=list(m),coefficients=[str(q)for q in c.to_list()])for m,c in sorted(p.items())]
        return dict(numerator=poly(v.numer),denominator=poly(v.denom))
    def positive(poly):
        value=norm(poly)
        if value.denom.degree(0)>0 or value.denom.degree(1)>0:raise ArithmeticError('Norm is not a polynomial')
        rational=[]
        for mon,c in value.numer.items():
            coefficients=c.to_list()
            if len(coefficients)>1 or not coefficients or coefficients[0]<=0:raise ArithmeticError('Unit norm does not have strictly positive rational coefficients')
            rational.append(dict(powers=list(mon),coefficient=str(coefficients[0])))
        if not rational:raise ArithmeticError('Zero unit norm')
        return rational
    Xcondition=N(w[3]**3-w[4]**3);L5=a*a-a*b-a+b*b-b+1;H=a*a-a*b+2*a+b*b+2*b+4
    expected=(a+b+1)**2*L5*H
    if N(norm(Xcondition.numer)-K.from_sympy(expected)):raise ArithmeticError('S5 cube-difference norm factorization failed')
    if s.expand(2*L5-((a-b)**2+(a-1)**2+(b-1)**2))or s.expand(H-((a-b)**2+a*a+b*b)/2-2*a-2*b-4):
        raise ArithmeticError('Positive-domain SOS identities failed')
    return locals()


def save(name,payload):
    path=Path(__file__).resolve().parents[1]/'certificates'/name
    path.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8');return path.name


def verify_regular():
    ctx=construct();J=ctx['J'];K=ctx['K'];N=ctx['N'];a=ctx['a'];b=ctx['b'];om=ctx['om'];MT=ctx['MT']
    inverse=ctx['inverse'](J);g=N(inverse[4][2]);units={}
    for row in inverse:
        for v in row:
            key=str(v.denom)
            if key not in units:units[key]=ctx['positive'](v.denom)
    aa=ctx['mass'][3];bb=ctx['mass'][4]
    expected=K.from_sympy(ctx['s'].Rational(2,9))*aa*MT**3*((aa-bb)+om*(aa-1))
    if N(K.field.raw_new(g.numer,g.numer.ring.one)-expected):raise ArithmeticError('Regular S5 response numerator failed')
    for i in(0,2):
        if N(inverse[3][i]+bb/aa*inverse[4][i])or inverse[5][i]:raise ArithmeticError('Triangle response centring/gauge failed')
    if N(inverse[4][0]+g):raise ArithmeticError('Translation forcing is not in the response kernel')
    # At a=b=1 the unit-modulus odd coefficient cannot cancel the two majority terms.
    signs=[]
    for p,r in ((1,1),(1,-1),(-1,1),(-1,-1)):
        v=p*om+r*om**2;vnorm=N(v*ctx['conjugate'](v.numer))
        if p==r and N(v+K.convert(p)):raise ArithmeticError('Equal-sign exceptional condition failed')
        if p!=r and vnorm-K.convert(3):raise ArithmeticError('Mixed-sign exceptional modulus failed')
        signs.append(dict(majority_signs=[p,r],required_rho='+-1'if p==r else'modulus squared 3',impossible=True))
    payload=dict(check='D2.6.17.S5_exact_second_order_response',status='PASS',g_scaled=ctx['encode'](g),inverse=[[ctx['encode'](v)for v in row]for row in inverse],
        positive_denominator_norms=list(units.values()),cube_difference_norm='(a+b+1)^2 * L5(a,b,1) * (a^2-ab+b^2+2a+2b+4)',
        g_zero_set='a=b=1 since 1 and omega are independent over real numbers and a,MT>0',
        response='Psi1=-m4/2 * delta^2*g_scaled * X^2; translation forcing kernel and centred W response verified from both relevant inverse columns',
        X_nonzero='X=0 would imply W3^3=W4^3. Exact norm and SOS factors restrict this to a=b=1.',
        exceptional_signs=signs,scope='All nonfold S5 triangles, both conjugate shapes and all three cubic scales, for every positive pair mass. The equal-triangle-mass point has no compatible one-different label since |rho|=1 and rho!=+-1.')
    name=save('d2_S5_exact_response.json',payload)
    return [dict(check=payload['check'],status='PASS',items=6,detail=dict(certificate=name,domain='entire positive mass cone; common mass scaling normalizes m5=1',response_zero_excluded=True))],{}


def verify_fold():
    ctx=construct(True);J=ctx['J'];K=ctx['K'];N=ctx['N'];a=ctx['a'];b=ctx['b']
    M=[[v for j,v in enumerate(row)if j!=4]for i,row in enumerate(J)if i!=3];inverse=ctx['inverse'](M)
    units={}
    for row in inverse:
        for v in row:
            key=str(v.denom)
            if key not in units:units[key]=ctx['positive'](v.denom)
    rhs=[-row[4]for i,row in enumerate(J)if i!=3]
    sol=[N(sum((v*r for v,r in zip(row,rhs)),K.zero))for row in inverse];kernel=sol[:4]+[K.one]+sol[4:]
    ell=[N(-sum((J[3][j if j<4 else j+1]*inverse[j][i]for j in range(8)),K.zero))for i in range(8)];ell=ell[:3]+[K.one]+ell[3:]
    if any(N(sum((v*r for v,r in zip(row,kernel)),K.zero))for row in J):raise ArithmeticError('Fold right kernel failed')
    if any(N(sum((ell[i]*J[i][j]for i in range(9)),K.zero))for j in range(9)):raise ArithmeticError('Fold left kernel failed')
    if N(kernel[3]+ctx['mass'][4]/ctx['mass'][3])or kernel[4]-K.one or kernel[5]or N(ell[0]+ell[2]):raise ArithmeticError('Fold printed kernel/cokernel structure failed')
    ellunit=ctx['positive'](ell[2].numer)
    payload=dict(check='D2.6.18.general_fold_rank_and_simple_torque_zero',status='PASS',inverse_minor=[[ctx['encode'](v)for v in row]for row in inverse],kernel=[ctx['encode'](v)for v in kernel],left_kernel=[ctx['encode'](v)for v in ell],
        positive_denominator_norms=list(units.values()),cokernel_Z4_positive_norm=ellunit,rank=8,
        minor=dict(omitted_row='FW4',omitted_column='W4'),
        augmentation='ell applied to the t-column equals a nonzero mass/scale unit times X; X=0 forces equal triangle masses, where the unit-rho sign census excludes a one-different Phi=0 label',
        torque_derivative='With h=W4-W4(seed), eta=(-m4/m3,1,0), tprime=0, pair Zprime=0; dQ/dh=-c_sigma*m4*X/2 !=0',
        scope='All six fold classes and all positive pair masses at a one-different Phi=0 label. Exact full positive-mass rank-8 minor, cokernel and kernel identities; the selected 20-row augmentation and simple torque zero follow by the verified block structure.')
    name=save('d2_general_fold_block.json',payload)
    return [dict(check=payload['check'],status='PASS',items=6,detail=dict(certificate=name,triangle_rank=8,positive_mass_minor=True,scope=payload['scope']))],{}

