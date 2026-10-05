"""D6 Lemma 9.9: exact leading jets at all eight corank-one orbits.

Rebuilds every jet from chart equations; rational arithmetic in Q(i,lambda0)
and the scale algebra d^3=1/[c(2A+lambda^3)], e^3=1/[c(2B+kappa^3)].
Certificates contain polynomial coefficients, never executable expressions.
"""
import json,time
from pathlib import Path
from smale_verify.symbolic import sympy
from smale_verify.exact import V,diff_poly
from smale_verify.d6_states import census


def units(K,lam,kap):
 s=K.s;A,B,c=K.A,K.B,K.c
 result=[A,B,c,2*A+lam**3,2*B+kap**3]
 for mass,root in ((A,lam),(B,kap)):
  if not root.is_Integer:
   H=2*mass+root**3
   result.append(18*(1-root**2)+9*root*H-s.Rational(9,4)*H**2)
 return result


def strip_units(K,poly,lam,kap):
 factors=[]
 for unit in units(K,lam,kap):
  divisor=K.base.from_sympy(unit).numer;count=0
  while True:
   common=poly.gcd(divisor)
   if not any(any(monomial) for monomial in common):break
   quotient,remainder=poly.div(common)
   if remainder:raise ArithmeticError('Nonexact rank-unit factor division')
   poly=quotient;count+=1
  if count:factors.append(dict(unit=str(unit),power=count))
 if not poly or any(any(monomial) for monomial in poly):
  raise ArithmeticError('Non-unit factor remains in leading coefficient or jet denominator: '+str(poly.as_expr()))
 return dict(factors=factors,nonzero_constant=str(poly.as_expr()))


def verify_leading_unit(K,lam,kap,leading):
 s=K.s;A,B,c=K.A,K.B,K.c
 if lam.is_Integer and kap.is_Integer:
  if lam==2 and kap==2:expected=-(2*A+2*B+1)/(60*c*(A+4)*(B+4))
  elif lam==2 and kap==-2:expected=(2*A+2*B+1)/(4*c*(A+4)*(B-4))
  elif lam==-2 and kap==-2:expected=(2*A+2*B+1)/(68*c*(A-4)*(B-4))
  else:raise ValueError('Unexpected Gaussian representative')
 elif lam.is_Integer:expected=A*kap**3-4*A*kap+B*kap**3-B*kap**2-4*B*kap+B
 elif not s.simplify(lam-kap):expected=1-2*kap**3
 else:expected=2*B/(2*B+kap**3)
 ratio=leading/K.from_sympy(expected)
 if len(ratio.terms)!=1:raise ArithmeticError('Leading ratio is not a scale monomial times a mass unit')
 scale,coefficient=next(iter(ratio.terms.items()))
 return dict(expected=str(expected),scale_powers=list(scale),
  numerator=strip_units(K,coefficient.numer,lam,kap),
  denominator=strip_units(K,coefficient.denom,lam,kap))


def build_case(rz,rw,rt,lam_name,kap_name,signs,jet_data=None):
 s=sympy()
 c,d,e=s.symbols('c d e')
 primitive=(s.sqrt(6)+s.I*s.sqrt(2))/2
 choices={'2':s.Integer(2),'-2':s.Integer(-2),'p':primitive,'-p':-primitive,'q':s.conjugate(primitive),'-q':-s.conjugate(primitive)}
 lam=choices[lam_name]
 kap=choices[kap_name]
 from smale_verify.cube_scales import ScaleAlgebra
 K=ScaleAlgebra(s,lam,kap)
 A,B=K.A,K.B
 def N(v):return v
 def total(values):
  result=K.zero
  for v in values:result=result+v
  return result
 signs=list(signs)
 degree=4 if lam.is_Integer and kap.is_Integer else 3 if lam.is_Integer or kap.is_Integer else 2
 aa=c*A;bb=c*B
 constant={f'm{j}':N(K.from_sympy(v)) for j,v in enumerate((aa,aa,bb,bb,c),1)}
 pt={'Z1':-s.Rational(1,2),'Z2':s.Rational(1,2),
  'W1':-d*(aa+2*c*lam),'W2':d*(aa+2*c*(lam**2/2-1)*lam),'W5':aa*d*(4/lam**2-1),
  'D12':d,'D15':lam*d,'D25':(lam**2/2-1)*lam*d,
  'W3':-s.Rational(1,2),'W4':s.Rational(1,2),
  'Z3':-e*(bb+2*c*kap),'Z4':e*(bb+2*c*(kap**2/2-1)*kap),'Z5':bb*e*(4/kap**2-1),
  'D34':e,'D35':kap*e,'D45':(kap**2/2-1)*kap*e,
  'D13':2*s.I*signs[0],'D14':2*signs[1],'D23':2*signs[2],'D24':2*s.I*signs[3]}
 pt={name:N(K.from_sympy(s.sympify(v))) for name,v in pt.items()}
 zrows=[rz[1],rz[2],rz[5],rw[1],rt[(1,2)],rt[(1,5)],rt[(2,5)],V('Z2')-V('Z1')-1]
 wrows=[rw[3],rw[5],rz[3],rz[4],rt[(3,4)],rt[(3,5)],rt[(4,5)],V('W4')-V('W3')-1]
 rows=zrows+wrows+[rt[p] for p in ((1,3),(1,4),(2,3),(2,4))]
 raw=list(rz.values())+list(rw.values())+list(rt.values())+[V('Z2')-V('Z1')-1,V('W4')-V('W3')-1]
 variables=['Z1','Z2','W1','W2','W5','D12','D15','D25','W3','W4','Z3','Z4','Z5','D34','D35','D45','D13','D14','D23','D24']
 def evalpoint(poly):
  result=K.zero
  for mon,coef in poly.terms.items():
   term=K.convert(coef)
   for name,power in mon:
    if name=='t':term*=K.zero**power
    else:term=N(term*(pt[name] if name in pt else constant[name])**power)
   result=N(result+term)
  return result
 J=[[evalpoint(diff_poly(row.coefficient('t',0),v)) for v in variables] for row in rows]
 def inverse(mat):
  n=len(mat);aug=[row+[K.one if i==j else K.zero for j in range(n)] for i,row in enumerate(mat)]
  for col in range(n):
   pivot=min((k for k in range(col,n) if aug[k][col]),key=lambda k:aug[k][col].cost());aug[col],aug[pivot]=aug[pivot],aug[col]
   q=aug[col][col];aug[col]=[N(x/q) for x in aug[col]]
   for j in range(n):
    if j!=col and aug[j][col]:
     q=aug[j][col];aug[j]=[N(x-N(q*y)) for x,y in zip(aug[j],aug[col])]
  return [row[n:] for row in aug]
 
 Ji=inverse(J) if jet_data is None else None
 if Ji is not None:
  for i in range(20):
   for j in range(20):
    value=total(N(J[i][k]*Ji[k][j]) for k in range(20))
    if value-K.convert(int(i==j)):raise ArithmeticError('Jacobian inverse identity failed')
 
 index={v:i for i,v in enumerate(variables)}
 def compile(poly):
  terms=[]
  for mon,coef in poly.terms.items():
   value=K.convert(coef);powers=[];delay=0
   for name,power in mon:
    if name=='t':delay=power
    elif name in index:powers.append((index[name],power))
    else:value=N(value*constant[name]**power)
   terms.append((value,delay,tuple(powers)))
  return terms
 compiled=[compile(row) for row in rows+raw]
 def mul(a,b,n):return [total(N(a[j]*b[k-j]) for j in range(k+1)) for k in range(n+1)]
 def evaluate(jets,n):
  one=[K.one]+[K.zero]*n;monomials={():one};powers={};out=[]
  for terms in compiled:
   value=[K.zero]*(n+1)
   for coef,delay,mon in terms:
    if delay>n:continue
    if mon not in monomials:
     product=one
     for idx,power in mon:
      if (idx,power) not in powers:
       v=one
       for _ in range(power):v=mul(v,jets[idx][:n+1],n)
       powers[(idx,power)]=v
      product=mul(product,powers[(idx,power)],n)
     monomials[mon]=product
    v=monomials[mon]
    for k in range(delay,n+1):value[k]=N(value[k]+N(coef*v[k-delay]))
   out.append(value)
  return out
 if jet_data is None:jets=[[pt[v]]+[K.zero]*degree for v in variables]
 else:
  if len(jet_data)!=20 or any(len(row)!=degree+1 for row in jet_data):raise ValueError('Wrong jet dimensions')
  jets=[[K.decode(value) for value in row] for row in jet_data]
  if any(jets[i][0]-pt[v] for i,v in enumerate(variables)):raise ArithmeticError('Jet constant is not the claimed atlas point')
 if any(v[0] for v in evaluate(jets,0)):raise ArithmeticError('Boundary point does not satisfy all 22 rows')
 for order in ([degree] if jet_data is not None else range(1,degree+1)):
  start=time.time()
  if Ji is not None:
   forcing=[-v[order] for v in evaluate(jets,order)[:20]]
   sol=[total(N(a*b) for a,b in zip(row,forcing)) for row in Ji]
   for j in range(20):jets[j][order]=sol[j]
  values=evaluate(jets,order)
  if any(any(v) for v in values[:20]):raise ArithmeticError('Selected jet rows do not vanish')
  print(f'D6 {lam_name},{kap_name} signs={signs}: exact order {order} replayed',flush=True)
 
 leading=values[26][degree]
 if any(values[26][:degree]) or not leading:raise ArithmeticError('Wrong leading order of Rw2')
 if lam.is_Integer and kap.is_Integer:
  if any(values[28][:degree]) or not values[28][degree]:raise ArithmeticError('Wrong Gaussian leading order of Rw4')
 unit=verify_leading_unit(K,lam,kap,leading)
 payload=dict(lam=lam_name,kap=kap_name,signs=signs,order=degree,variables=variables,
  jets=[[K.encode(v) for v in row] for row in jets],
  raw_residual=[[K.encode(v) for v in row] for row in values[20:]],leading_unit=unit)
 # Re-read the safe polynomial encoding and compare every coefficient exactly.
 denominator_cache=set()
 for original,encoded in zip(jets,payload['jets']):
  for value,data in zip(original,encoded):
   if jet_data is None and K.decode(data)-value:raise ArithmeticError('Jet certificate encoding changed a coefficient')
   for coefficient in value.terms.values():
    key=tuple(sorted(coefficient.denom.items()))
    if key not in denominator_cache:
     strip_units(K,coefficient.denom,lam,kap);denominator_cache.add(key)
 return payload


def verify_leading_jets(rz,rw,rt):
 orbits=[orbit for orbit in census() if orbit['corank_one']]
 names=['2','-2','p','-p','q','-q']
 seeds=[tuple(o['seed']) for o in orbits]
 expected=[(0,0,-1,-1,1,1),(0,1,-1,-1,1,1),(0,2,-1,-1,1,1),
           (1,1,-1,-1,1,1),(1,2,-1,-1,1,1),(2,2,-1,-1,1,1),
           (2,3,-1,-1,1,1),(2,3,-1,1,-1,1)]
 if seeds!=expected:raise ArithmeticError('Corank orbit representatives changed')
 cases=[]
 path=Path(__file__).resolve().parents[1]/'certificates'/'d6_leading_jets.json'
 supplied=json.loads(path.read_text(encoding='utf-8'))['cases'] if path.is_file() else None
 if supplied is not None and len(supplied)!=8:raise ValueError('Certificate needs exactly eight corank orbit cases')
 for orbit,seed in zip(orbits,seeds):
  ref=None if supplied is None else supplied[len(cases)]
  if ref is not None and (ref['state']!=list(seed) or ref['lam']!=names[seed[0]] or ref['kap']!=names[seed[1]] or ref['signs']!=list(seed[2:])):
   raise ValueError('Certificate orbit labels mismatch')
  case=build_case(rz,rw,rt,names[seed[0]],names[seed[1]],seed[2:],None if ref is None else ref['jets'])
  case['state']=list(seed);case['orbit_size']=orbit['size']
  cases.append(case)
 certificate=dict(check='D6.9.9.leading_jets_all_eight_orbits',status='PASS',
  arithmetic='exact Q(i,lambda0)(A,B,c) with two cubic scale relations',
  algebraic_generator='theta=i+lambda0; theta^8+14theta^4-48theta^2+49=0',
  coefficient_encoding='descending rational coefficients in theta; sparse mass powers A,B,c; scale powers d,e in {0,1,2}',
  cases=cases)
 path.write_text(json.dumps(certificate,indent=2)+'\n',encoding='utf-8')
 counts={4:0,3:0,2:0}
 for case in cases:counts[case['order']]+=9*case['orbit_size']
 if counts!={4:144,3:576,2:576}:raise ArithmeticError('D6 leading-jet label census mismatch')
 return [dict(check=certificate['check'],status='PASS',items=8,
  detail=dict(certificate=path.name,orders_and_generic_labels=counts,
   selected_rows=20,substitution_rows=22,jet_denominators='only nonzero finite-atlas scale and rank units',
   leading='each exact quotient by the printed factor is a verified nonzero mass/scale unit',
   coverage='eight corank orbits; nine cubic scale choices per state; exceptional mass faces delete absent Gaussian scale classes'))],{}



