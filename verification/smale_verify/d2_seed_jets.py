"""Safe exact replay of equal-mass D2 ordinary and fold seed certificates."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
import json,hashlib
from .exact import V,normalized_system,equations,diff_poly,substitute
from .d2_fields import s,nu,delta,CF,decode_new
from .d2_laurent import LF,L,q,decode,regular


def build(kind):
    base=normalized_system();m=base['m'];D=base['D'];t=V('t')
    Z={j:V('Z'+str(j))for j in m};W={j:V('W'+str(j))for j in m}
    z={j:t**(-1 if j<3 else 1)*Z[j]for j in m};w={j:t**(1 if j<3 else -1)*W[j]for j in m}
    ds={p:D[p]if p==(1,2)or p[0]>2 else t*D[p]for p in D}
    ez,ew,rt=equations(m,z,w,ds)
    rz={j:t**(1 if j<3 else -1)*ez[j]for j in m};rw={j:t**(-1 if j<3 else 1)*ew[j]for j in m}
    raw=list(rz.values())+list(rw.values())+list(rt.values())+[Z[2]-Z[1]-1,W[5]-1]
    a=4*nu**2/(1+nu**2)**2;b=(1-nu**2)**2/(1+nu**2)**2;omega=(-1+s.I*s.sqrt(3))/2
    constants={f'm{j}':LF.from_sympy(value)for j,value in enumerate((a,b,1,1,1),1)}
    constants['m1']*=q**3;constants['m2']*=q**3
    pt={'Z1':-b,'Z2':a,'W1':-b,'W2':a,'D12':1,'D34':delta,
        'D35':omega*delta if kind=='fold'else delta,'D45':omega**2*delta if kind=='fold'else delta,
        'Z3':-delta*omega if kind=='fold'else delta*omega**2,
        'Z4':-delta*omega**2 if kind=='fold'else delta*omega,
        'Z5':-delta if kind=='fold'else delta,'W3':omega,'W4':omega**2,'W5':1}
    for j,v in zip((3,4,5),(omega,omega**2,1)):
        pt[f'D1{j}']=(1+nu**2)/(1-nu**2)*v;pt[f'D2{j}']=-s.I*(1+nu**2)/(2*nu)*v
    pt={name:LF.from_sympy(s.sympify(value))for name,value in pt.items()}
    pt['W1']*=q**2;pt['W2']*=q**2;pt['D12']*=q**-1
    variables=[f'Z{j}'for j in ((1,2,4,5)if kind=='fold'else range(1,6))]+[f'W{j}'for j in range(1,6)]+[f'D{i}{j}'for i,j in D]
    omitted=(6,8,9)if kind=='fold'else(1,4)
    return dict(raw=raw,constants=constants,point=pt,variables=variables,
                allvars=variables+(['Z3']if kind=='fold'else[]),selected=[r for i,r in enumerate(raw)if i not in omitted],omitted=omitted,omega=omega)


def eval_point(poly,seed):
    out=LF.zero
    for mon,coefficient in poly.terms.items():
        term=LF.convert(coefficient)
        for name,power in mon:
            value=LF.zero if name=='t'else seed['point'].get(name,seed['constants'].get(name))
            term=term*value**power
        out+=term
    return out


def verify_inverse(seed,reference):
    variables=seed['variables'];n=len(variables)
    if reference['variables']!=variables:raise ValueError('Wrong inverse coordinate order')
    inverse=[[decode(v)for v in row]for row in reference['inverse']]
    if len(inverse)!=n or any(len(row)!=n for row in inverse):raise ValueError('Wrong inverse dimensions')
    J=[[eval_point(diff_poly(row.coefficient('t',0),name),seed)for name in variables]for row in seed['selected']]
    for row in inverse:
        for value in row:regular(value)
    for i in range(n):
        for j in range(n):
            if sum((J[i][k]*inverse[k][j]for k in range(n)),LF.zero)-LF.convert(int(i==j)):
                raise ArithmeticError('Seed Jacobian inverse identity failed')
    return n


def replay_series(seed,reference,order,tScale,tDelay,free=None):
    allvars=seed['allvars'];index={name:j for j,name in enumerate(allvars)}
    if reference['variables']!=allvars or reference['order']!=order:raise ValueError('Wrong seed series metadata')
    jets=[]
    if len(reference['jets'])!=len(allvars):raise ValueError('Incomplete jet coordinate list')
    for name,raw in zip(allvars,reference['jets']):
        jet={}
        for term in raw:
            mon=tuple(term['powers'])
            if len(mon)!=2 or any(type(v)is not int or v<0 for v in mon)or sum(mon)>order or mon in jet:raise ValueError('Invalid jet monomial')
            jet[mon]=decode(term['coefficient']);regular(jet[mon])
        if jet.get((0,0),LF.zero)-seed['point'][name]:raise ValueError('Certificate constant is not the source seed')
        jets.append(jet)
    if free is not None:
        expected={(0,0):seed['point']['Z3'],free[0]:free[1]}
        if any(jets[-1].get(p,LF.zero)-expected.get(p,LF.zero)for p in set(jets[-1])|set(expected)):
            raise ValueError('Incorrect free-coordinate path')
    compiled=[]
    for poly in seed['raw']:
        terms=[]
        for mon,coefficient in poly.terms.items():
            factor=LF.convert(coefficient);delay=0;powers=[]
            for name,power in mon:
                if name=='t':delay=power*tDelay;factor*=tScale**power
                elif name in index:powers.append((index[name],power))
                else:factor*=seed['constants'][name]**power
            terms.append((factor,delay,tuple(powers)))
        compiled.append(terms)
    def multiply(a,b):
        out={}
        for (h,t),x in a.items():
            for (k,u),y in b.items():
                if h+t+k+u<=order:out[(h+k,t+u)]=out.get((h+k,t+u),LF.zero)+x*y
        return {p:v for p,v in out.items()if v}
    one={(0,0):LF.one};monomials={():one};powers={};values=[]
    for terms in compiled:
        value={}
        for factor,delay,mon in terms:
            if delay>order:continue
            if mon not in monomials:
                prod=one
                for index_,power in mon:
                    if (index_,power)not in powers:
                        v=one
                        for _ in range(power):v=multiply(v,jets[index_])
                        powers[(index_,power)]=v
                    prod=multiply(prod,powers[(index_,power)])
                monomials[mon]=prod
            for (h,t),coefficient in monomials[mon].items():
                if h+t+delay<=order:value[(h,t+delay)]=value.get((h,t+delay),LF.zero)+factor*coefficient
        values.append({p:v for p,v in value.items()if v})
    if any(values[i]for i in range(22)if i not in seed['omitted']):raise ArithmeticError('A selected source row has a nonzero residual')
    g=[values[i]for i in seed['omitted']]
    if len(reference['g'])!=len(g):raise ValueError('Incorrect omitted residual list')
    for row,encoded in zip(g,reference['g']):
        supplied={tuple(r['powers']):decode(r['coefficient'])for r in encoded}
        if any(row.get(p,LF.zero)-supplied.get(p,LF.zero)for p in set(row)|set(supplied)):
            raise ArithmeticError('Supplied omitted residual differs from exact source replay')
    return g


def verify_ordinary():
    root=Path(__file__).resolve().parents[1];path=root/'certificates'/'d2_ordinary_seed_inputs.json'
    data=json.loads(path.read_text(encoding='utf-8'));seed=build('ordinary')
    verify_inverse(seed,data['inverse'])
    g=replay_series(seed,data['series'],3,LF.one,1)
    if any(value for mon,value in g[1].items()if sum(mon)<3):raise ArithmeticError('Ordinary omitted row has a lower-order obstruction')
    chi=(1-nu**2)/(2*nu);R1=18*nu**2*(1-nu**2)/(1+nu**2)**3;R2=6*nu**2*(1+nu**2)/(1-nu**2)**3
    expected=q**5*LF.from_sympy(R1*(1+s.I*chi))+q**3*LF.from_sympy(3*delta*R2*(1+s.I*chi**5))
    if g[1].get((0,3),LF.zero)-expected:raise ArithmeticError('Printed ordinary coefficient failed in the conjugate-i seed')
    x=s.symbols('x',positive=True)
    if s.expand((1+x**6)-(x+x**5)-(1-x)*(1-x**5)):raise ArithmeticError('Phase-sector inequality identity failed')
    payload=dict(check='D2.6.13.ordinary_general_mass_kappa3',status='PASS',input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        Jacobian_rank=20,orders_replayed=3,coefficient=expected.encode(),
        source_transport='i -> -i and sqrt(3) -> -sqrt(3) fixes omega and delta, and gives the printed 1-i chi coefficient. Independent delta and q cubic embeddings give all nine root choices.',
        nonvanishing='R1,R2,q,delta positive. (1+x^6)-(x+x^5)=(1-x)(1-x^5)>=0 for x>0 gives a phase sector of width pi/2. Its three rotations by 2pi/3 avoid -1, so the two terms cannot cancel.',
        domain='0<nu<1 and q=M^(1/3)>0; all coefficient and inverse denominators stripped into proved nonzero units')
    output=root/'certificates'/'d2_ordinary_seed.json';output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=72,detail=dict(certificate=output.name,scope='ordinary equal-triangle-mass labels, general positive pair mass, all cubic root embeddings'))],{}


def verify_fold():
    root=Path(__file__).resolve().parents[1];path=root/'certificates'/'d2_fold_seed_inputs.json'
    data=json.loads(path.read_text(encoding='utf-8'));seed=build('fold');verify_inverse(seed,data['inverse'])
    quadratic=replay_series(seed,data['quadratic'],2,q**-3,1,((1,0),LF.one))
    gamma11=quadratic[1].get((1,1),LF.zero);gamma02=quadratic[1].get((0,2),LF.zero)
    U=nu**4-2*s.I*nu**3-6*nu**2+2*s.I*nu+1
    printedgamma=LF.from_sympy(9*s.I*delta*(nu-s.I)*(nu+s.I)**3*U/(4*nu*(nu**2-1)*(nu**2+1)**2))
    if gamma11+printedgamma:raise ArithmeticError('Fold gamma11 differs from the reversed source seed')
    n5=quadratic[2].get((1,1),LF.zero)/gamma11;n=[LF.zero,LF.one,n5]
    for i,row in enumerate(quadratic):
        for mon,value in row.items():
            if sum(mon)<2 or mon==(2,0):raise ArithmeticError('Fold quadratic has an unexpected lower or a^2 term')
            if value-n[i]*quadratic[1].get(mon,LF.zero):raise ArithmeticError('Fold quadratic coefficient matrix is not rank one')
    b2=decode(data['boundary_b2']);regular(b2)
    slanted=replay_series(seed,data['slanted'],3,gamma11*q**-3,1,((0,1),-gamma02))
    boundary=replay_series(seed,data['boundary'],4,b2*q**-3,2,((0,1),LF.one))
    # The first omitted row is M^-2/3 psi2 in the source's tilde g.
    slanted[0]={p:q**-2*v for p,v in slanted[0].items()};boundary[0]={p:q**-2*v for p,v in boundary[0].items()}
    if any(v for row in slanted for p,v in row.items()if sum(p)<3):raise ArithmeticError('Slanted fold path has a lower residual')
    if any(v for row in boundary for p,v in row.items()if sum(p)<4):raise ArithmeticError('Boundary fold compatibility failed')
    mu3=-gamma11**2*slanted[0].get((0,3),LF.zero);mu4=-gamma11*boundary[0].get((0,4),LF.zero)
    expected3=q**-9*LF.from_sympy(-3**10*(nu-s.I)**4*(nu+s.I)**16*U**5/(2**12*delta*nu**6*(nu**2-1)**6*(nu**2+1)**8))
    expected4=q**-5*LF.from_sympy(-27*s.I*delta*(nu-s.I)**3*(nu+s.I)**5*U/(64*nu**3*(nu**2-1)**3))
    if mu3-expected3 or mu4+expected4:raise ArithmeticError('General-mass fold minor formula failed')
    # Exact row covariance explains gamma/minor signs in the source seed.
    replacements={'t':-V('t'),**{'D'+str(i)+str(j):-V('D'+str(i)+str(j))for i in (1,2)for j in(3,4,5)}}
    for row in seed['raw']:
        if (substitute(row,replacements)-row).terms:raise ArithmeticError('Reversal is not a source chart symmetry')
    if s.expand(U-(3*(nu-s.I)**4+(nu+s.I)**4)/4):raise ArithmeticError('Upsilon nonvanishing identity failed')
    payload=dict(check='D2.6.13.fold_general_mass_quadratic_and_two_paths',status='PASS',input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),Jacobian_rank=19,
        kernel_variables=['a=Z3-Z3(seed)','T=M*t'],orders=dict(quadratic=2,slanted=3,boundary=4),gamma11=gamma11.encode(),gamma02=gamma02.encode(),boundary_b2=b2.encode(),mu3=mu3.encode(),mu4=mu4.encode(),
        source_transport='The exact t/crossed-distance reversal gives printed gamma11 and mu4; mu3 is unchanged. The source normalized psi2=M^-2/3*psi2 is applied before forming minors.',
        nonvanishing='Upsilon=(3(nu-i)^4+(nu+i)^4)/4; the first summand has three times the modulus of the second for real nu. All other numerator/denominator factors are nonzero for 0<nu<1, q>0.',
        scope='General positive pair mass: exact quadratic rank one and both required strict-transform obstruction minors; all 22 source rows replayed on each path, 19 selected rows vanish to the required order.')
    output=root/'certificates'/'d2_fold_seed.json';output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=72,detail=dict(certificate=output.name,rank=19,source_minor_scaling=['M^-3 mu3','M^-5/3 mu4'],scope='rank-19 fold seed; finite label orbit transport is a separate check'))],{}
