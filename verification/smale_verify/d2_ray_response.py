"""All 72 same-sign ray triangles, with certified second-order responses.

Certificate guesses supply power-root disks and rational preconditioners only.
The source models, finite label list and all bounds are rebuilt by the checker.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product,permutations,combinations
import hashlib,json
from .symbolic import sympy
from .interval import R,Z,decode,certify_power_root,verified_linear_solve
from .interval_refine import refine_real_root
from .exact import V,equations,diff_poly


def symbolic_identities():
    s=sympy();a,b,x=s.symbols('a b x');A=-(b*x+1)/a;B=x
    cube=-A/(-b/(B-A)**2+1/(1-A)**2)
    from ak_resultants import shape_factors
    shapes,_,_=shape_factors()
    def expression(poly):
        return sum(s.Rational(c)*s.prod({'shape_w3':A,'shape_w4':B,'shape_w5':s.Integer(1)}[name]**power for name,power in mon)for mon,c in poly.terms.items())
    shape=s.Poly(s.cancel(expression(shapes[1])).as_numer_denom()[0],x)
    remaining=s.cancel(B+cube*(a/(B-A)**2+1/(1-B)**2)).as_numer_denom()[0]
    if s.rem(s.Poly(remaining,x),shape):raise ArithmeticError('Euler force compatibility does not reduce to S1')
    # S8 is the product of the two class quadratics over Q(omega).
    k,w=s.symbols('k w');quad1=k**2+w*k+w;quad2=k**2+w**2*k+w**2
    if s.rem(s.Poly(s.expand(quad1*quad2-(k**4-k**3+2*k+1)),w),s.Poly(w**2+w+1,w)):
        raise ArithmeticError('S8 quadratic class factorization failed')
    # Rational signed square roots, with xi5=1, follow from Lambda=0.
    A,B,c,d,e=s.symbols('A B c d e');f3=c**2/A;f4=d**2/B;f5=e**2
    Lambda=f3**2+f4**2+f5**2-2*f3*f4-2*f3*f5-2*f4*f5
    xi3=(f4-f3-f5)/(2*c*e);xi4=(-e-c*xi3)/d
    if s.cancel(xi3**2-1/A-Lambda/(4*c**2*e**2)) or s.cancel(xi4**2-1/B-Lambda/(4*d**2*e**2)):
        raise ArithmeticError('Signed-root Lambda identities failed')


def replay(d2_mass_resultants,d2_mass_rays):
    symbolic_identities()
    root=Path(__file__).resolve().parents[1];path=root/'certificates'/'d2_ray_response_inputs.json'
    inputs=json.loads(path.read_text(encoding='utf-8'))
    if len(inputs['families'])!=2 or {r['j']for r in inputs['families']}!={1,8}:raise ValueError('Incomplete ray families')
    sr=inputs['sqrt3_disk'];sqrt3,proof=certify_power_root(Z(3),2,sr['center'],sr['radius'])
    if Q(sr['center'][1]) or sqrt3.re.lo<=0:raise ValueError('sqrt3 must select the positive real branch')
    sqrt3=Z(sqrt3.re);omega=(-1+Z(0,1)*sqrt3)/2
    masses={j:V('m'+str(j))for j in (3,4,5)};z={j:V('Z'+str(j))for j in masses};w={j:V('W'+str(j))for j in masses}
    ds={p:V('D'+''.join(map(str,p)))for p in combinations(masses,2)}
    ez,ew,rt=equations(masses,z,w,ds)
    rows=[ez[3],ew[3],ez[4],ew[4],ew[5],*rt.values(),w[5]-1]
    variables=[*('Z'+str(j)for j in masses),*('W'+str(j)for j in masses),*('D'+''.join(map(str,p))for p in ds)]
    derivatives=[[diff_poly(row,name)for name in variables]for row in rows]
    def evaluate(poly,values):
        out=Z(0)
        for mon,c in poly.terms.items():
            term=Z(c)
            for name,power in mon:term=term*values[name]**power
            out=out+term
        return out
    references=json.loads((root/'certificates'/'d2_mass_ray_root_inputs.json').read_text(encoding='utf-8'))
    records=[]
    for ref in inputs['families']:
        j=ref['j'];ray=next(r for r in d2_mass_rays['root_records']if r['j']==j)
        case=next(r['krawczyk_case']for r in references if r['j']==j)
        a,b,u,v=[z.re for z in refine_real_root(case,ray['krawczyk']['krawczyk_image'])]
        mass={3:Z(a),4:Z(b),5:Z(1)};W={4:Z(u,v),5:Z(1)};W[3]=-(mass[4]*W[4]+1)/mass[3]
        MT=sum(mass.values(),Z(0));rootdisks=ref['power_root_disks']
        if len(rootdisks)!=1:raise ValueError('Exactly one cubic-scale disk is required per family')
        # Every other L factor is nonzero throughout this certified mass box.
        exclusions=[]
        L5=((Z(a)-Z(b))**2+(Z(b)-1)**2+(1-Z(a))**2)/2
        if L5.re.lo<=0:raise ArithmeticError('L5 is not excluded on the ray box')
        exclusions.append(dict(factor='L5',range=L5.encode()))
        for index,data in d2_mass_resultants['resultants'].items():
            if int(index)==j:continue
            terms=data['L_factors'][-1]['polynomial'];value=Z(0)
            for term in terms:
                value+=Q(term['coefficient'])*Z(a)**term['powers'].get('m3',0)*Z(b)**term['powers'].get('m4',0)
            if not value.lower():raise ArithmeticError('Other L factor is not excluded on the ray box')
            exclusions.append(dict(factor='L5prime'if index=='5'else'L'+index,range=value.encode()))
        if j==1:
            signs={(3,4):-1,(3,5):1,(4,5):1}
            c3=-W[3]/sum((mass[k]*signs[(3,k)]/(W[k]-W[3])**2 for k in (4,5)),Z(0))
            c,scaleproof=certify_power_root(c3,3,rootdisks[0]['center'],rootdisks[0]['radius'])
            D={p:signs[p]*c/(W[p[1]]-W[p[0]])for p in ds};Zs={k:W[k]/c**2 for k in masses}
        else:
            delta,scaleproof=certify_power_root(1/MT,3,rootdisks[0]['center'],rootdisks[0]['radius'])
            kappa=(W[4]-W[3])/(1-W[4]);options=[]
            for e in (1,2):
                if not (kappa**2+omega**e*kappa+omega**e).lower():options.append(e)
            if len(options)!=1 or ref['lagrange_exponent']!=options[0]:raise ArithmeticError('S8 class selection failed')
            D={(3,4):delta*omega**options[0],(3,5):delta,(4,5):delta}
            dz4=1/(D[(3,4)]**2*(W[4]-W[3]));dz5=1/(D[(3,5)]**2*(1-W[3]))
            z3=-(mass[4]*dz4+dz5)/MT;Zs={3:z3,4:z3+dz4,5:z3+dz5}
        expected=set(product((False,True),range(3),permutations(masses)))
        supplied=[(r['conjugated'],r['scale'],tuple(r['permutation']))for r in ref['cases']]
        if len(supplied)!=36 or len(set(supplied))!=36 or set(supplied)!=expected:raise ValueError('Ray label transport list is incomplete')
        label_records=[]
        for supplied_case in ref['cases']:
            conjugated=supplied_case['conjugated'];scale=supplied_case['scale'];permutation=supplied_case['permutation']
            massp={k:mass[old]for k,old in zip(masses,permutation)};factor=W[permutation[2]]
            Wp={k:W[old]/factor for k,old in zip(masses,permutation)}
            Zp={k:Zs[old]*factor*omega**scale for k,old in zip(masses,permutation)}
            Dp={p:D[tuple(sorted((permutation[p[0]-3],permutation[p[1]-3])))]*omega**scale for p in ds}
            if conjugated:
                Wp={k:q.conjugate()for k,q in Wp.items()};Zp={k:q.conjugate()for k,q in Zp.items()};Dp={k:q.conjugate()for k,q in Dp.items()}
            xi={5:Z(1)};xi[3]=(massp[4]**2/Wp[4]-massp[3]**2/Wp[3]-massp[5]**2)/(2*massp[3]*massp[5]);xi[4]=(-massp[5]-massp[3]*xi[3])/massp[4]
            values={**{'m'+str(k):q for k,q in massp.items()},**{'Z'+str(k):q for k,q in Zp.items()},**{'W'+str(k):q for k,q in Wp.items()},**{'D'+''.join(map(str,p)):q for p,q in Dp.items()}}
            J=[[evaluate(poly,values)for poly in row]for row in derivatives]
            rhs=[xi[3]/Wp[3],Z(0),xi[4]/Wp[4],*[Z(0)]*6]
            inverse=[[decode(v)for v in row]for row in supplied_case['preconditioner']]
            sol,linear=verified_linear_solve(J,rhs,inverse)
            psi=-sum((massp[k]*xi[k]*sol[k]/Wp[k]for k in masses),Z(0))/2
            margin=Q(3,10)if j==1 else Q(36,100)
            lower=psi.modulus_lower()
            if lower<=margin:raise ArithmeticError('Printed ray response margin not certified')
            label_records.append(dict(conjugated=conjugated,scale=scale,permutation=permutation,response=psi.encode(),lower_bound=str(lower),linear=linear))
        records.append(dict(j=j,scale_root=scaleproof,other_factor_exclusions=exclusions,labels=label_records,minimum_lower_bound=str(min(Q(r['lower_bound'])for r in label_records))))
        print('D2 L'+str(j)+': all 36 triangle responses certified',flush=True)
    payload=dict(check='D2.6.15.all_72_ray_triangle_responses',status='PASS',sqrt3=proof,families=records,
        input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),arithmetic='exact fractions with outward dyadic rounding to 120 bits; rational Rouche and verified Neumann inverses',
        completeness='All other L factors excluded; each common-shape pair has three cubic-scale choices; all six body permutations and conjugation explicitly checked; signed-root vector unique up to global sign.')
    output=root/'certificates'/'d2_ray_response.json';output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=72,detail=dict(certificate=output.name,Euler_margin='>0.3',Lagrange_margin='>0.36',triangle_linear_system='all nine source rows reconstructed and interval inverse certified'))],{}
