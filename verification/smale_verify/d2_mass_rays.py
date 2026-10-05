"""Global D2 L1/L8 ray classification using Bernstein, Rouche and Krawczyk.

Every certificate polynomial is reconstructed from the manuscript's exact
mass resultant inputs. Numerical root guesses are not used by this checker.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import permutations
import hashlib,json
from .symbolic import sympy
from .bernstein import coefficients,corner_sectors,positive_cover
from .exact import P,V,substitute
from verify_interval_certificates import C,I,q,point_complex,box_complex,verify_uniform_rouche,verify_krawczyk


def path_box(path):
 lo=[Q(0),Q(0)];hi=[Q(1),Q(1)]
 for i,bit in enumerate(path):
  axis=i%2;mid=(lo[axis]+hi[axis])/2
  if bit=='0':hi[axis]=mid
  elif bit=='1':lo[axis]=mid
  else:raise ValueError('Invalid Bernstein binary path')
 return list(zip(lo,hi))


def classify(d2_mass_resultants):
 root=Path(__file__).resolve().parents[1]
 inputs=root/'certificates'/'d2_mass_ray_root_inputs.json'
 references=json.loads(inputs.read_text(encoding='utf-8'))
 if len(references)!=2 or {r['j'] for r in references}!={1,8}:raise ValueError('Need exactly L1 and L8 ray inputs')
 s=sympy();a,b,u,v,x=s.symbols('a b u v x',real=True);records=[]
 def decode(terms):return sum(s.Rational(t['coefficient'])*a**t['powers'].get('m3',0)*b**t['powers'].get('m4',0)*x**t['powers'].get('x',0) for t in terms)
 def encode(poly):return [dict(coefficient=[str(c),0],powers=list(p)) for p,c in s.Poly(poly,a,b,u,v).terms()]
 def evaluate(terms,mass_box,xbox):
  result=C(0)
  for (pa,pb,px),coef in terms:
   result=result+C(I(*mass_box[0]),0)**pa*C(I(*mass_box[1]),0)**pb*xbox**px*int(coef)
  return result
 global_covers=[]
 for j in (1,8):
  ref=next(r for r in references if r['j']==j)
  massdata=d2_mass_resultants['resultants'][str(j)]
  covers=[]
  for sector in corner_sectors(massdata['L_factors'][-1]['polynomial']):
   power=sector.pop('power_terms');cover=positive_cover(coefficients(power),52)
   if cover['sign']!=1:raise ArithmeticError('L polynomial sign changed')
   covers.append(dict(sector=sector,cover=cover))
  holes=[r for r in covers if r['cover']['unresolved']]
  if len(holes)!=2:raise ArithmeticError('Expected exactly two symmetric positive-mass ray neighborhoods')
  target={'m3':'1','m5':'r','m4':'r*s'} if j==1 else {'m5':'1','m4':'r','m3':'r*s'}
  chosen=next(r for r in holes if r['sector']['substitution']==target)
  if holes[0]['cover']['unresolved']!=holes[1]['cover']['unresolved']:
   # Minimum bounds may differ by source normalization; domains must coincide.
   if [r['path'] for r in holes[0]['cover']['unresolved']]!=[r['path'] for r in holes[1]['cover']['unresolved']]:
    raise ArithmeticError('Symmetric ray holes do not have matching domains')
  boxes=[path_box(r['path']) for r in chosen['cover']['unresolved']]
  rlo=min(z[0][0] for z in boxes);rhi=max(z[0][1] for z in boxes)
  slo=min(z[1][0] for z in boxes);shi=max(z[1][1] for z in boxes)
  mass_box=[(1/rhi,1/rlo),(slo,shi)] if j==1 else [(rlo*slo,rhi*shi),(rlo,rhi)]
  f=decode(massdata['shape_flat']);g=decode(d2_mass_resultants['lambda_flat'])
  shape_terms=s.Poly(f,a,b,x).terms();lambda_terms=s.Poly(g,a,b,x).terms()
  disks=[dict(center=r['center'],radius=r['radius']) for r in ref['shape_discs']]
  rouche=dict(parameters=2,kind='uniform_rouche',
   parameter_box=[dict(re=[str(lo),str(hi)],im=[0,0]) for lo,hi in mass_box],
   polynomial=[dict(coefficient=[str(c),0],powers=list(p)) for p,c in shape_terms],disks=disks)
  rouche_result=verify_uniform_rouche(rouche)
  candidates=[];disk_results=[]
  for disk in disks:
   center=point_complex(disk['center'],'shape disk');radius=q(disk['radius'])
   xbox=C(I(center.re.lo-radius,center.re.lo+radius),I(center.im.lo-radius,center.im.lo+radius))
   value=evaluate(lambda_terms,mass_box,xbox)
   excludes=not(value.re.lo<=0<=value.re.hi and value.im.lo<=0<=value.im.hi)
   disk_results.append(dict(disk=disk,Lambda_excluded=excludes,
       Lambda_range=dict(re=[str(value.re.lo),str(value.re.hi)],im=[str(value.im.lo),str(value.im.hi)])))
   if not excludes:candidates.append((disk,xbox))
  if len(candidates)!=2:raise ArithmeticError('Common shapes are not confined to one conjugate pair')
  plus=next((r for r in candidates if r[1].im.lo>0),None)
  minus=next((r for r in candidates if r[1].im.hi<0),None)
  if plus is None or minus is None:raise ArithmeticError('Common-shape disks are not nonreal')
  if plus[0]['radius']!=minus[0]['radius'] or q(plus[0]['center'][0])!=q(minus[0]['center'][0]) or q(plus[0]['center'][1])!=-q(minus[0]['center'][1]):
   raise ArithmeticError('Common disks are not exact conjugates')
  F4=[s.re(s.expand_complex(f.subs(x,u+s.I*v))),s.im(s.expand_complex(f.subs(x,u+s.I*v))),
      s.re(s.expand_complex(g.subs(x,u+s.I*v))),s.im(s.expand_complex(g.subs(x,u+s.I*v)))]
  case=ref['krawczyk_case']
  if case['variables']!=4 or case['polynomials']!=[encode(poly) for poly in F4]:raise ValueError('Krawczyk polynomial model differs from S/Lambda source')
  casebox=[box_complex(item,'Krawczyk box') for item in case['box']]
  required=mass_box+[(plus[1].re.lo,plus[1].re.hi),(plus[1].im.lo,plus[1].im.hi)]
  for box,(lo,hi),center in zip(casebox,required,case['center']):
   if box.re.lo>lo or box.re.hi<hi:raise ValueError('Krawczyk box misses part of the global mass/shape hole')
   if box.im.lo!=-box.im.hi or q(center[1])!=0:raise ValueError('Reality requires conjugation-symmetric boxes and real centers')
  if min(casebox[0].re.lo,casebox[1].re.lo,casebox[3].re.lo)<=0:raise ValueError('Mass/shape sign not certified')
  krawczyk=verify_krawczyk(case)
  for k,poly in enumerate(F4):
   if s.Poly(poly.subs(v,-v)-(-1 if k%2 else 1)*poly,a,b,u,v):raise ArithmeticError('Conjugate-shape transport failed')
  # Reality follows from uniqueness and complex conjugation of all four
  # polynomial variables. The lower shape is the v -> -v image of this root.
  records.append(dict(j=j,mass_hole=[[str(z) for z in interval] for interval in mass_box],
      uniform_rouche_case=rouche,uniform_rouche_result=rouche_result,
      shape_disk_results=disk_results,krawczyk=krawczyk,
      real_solution='uniqueness in a conjugation-symmetric complex box',
      shape_pair='exactly one non-real conjugate pair',
      mass_ranges=krawczyk['krawczyk_image'][:2]))
  global_covers.append(dict(j=j,sectors=covers))
  print(f'D2 L{j}: complete mass cover, all shape roots, and unique mass/shape root certified',flush=True)
 # Check the polynomial mass symmetries underlying both second rays and
 # the six permutations in each three-factor family.
 def sparse(terms):return P({tuple(sorted(t['powers'].items())):Q(t['coefficient']) for t in terms})
 core={int(j):sparse(r['L_factors'][-1]['polynomial']) for j,r in d2_mass_resultants['resultants'].items()}
 symmetry=[]
 for seed,family in ((1,(1,2,3)),(8,(6,7,8))):
  for perm in permutations(('m3','m4','m5')):
   transformed=substitute(core[seed],{name:V(target) for name,target in zip(('m3','m4','m5'),perm)})
   matches=[j for j in family if not (transformed-core[j]).terms or not (transformed+core[j]).terms]
   if len(matches)!=1:raise ArithmeticError('Mass-factor permutation does not match its printed family')
   symmetry.append(dict(seed=seed,permutation=list(perm),factor=matches[0]))
 for seed,pair in ((1,('m4','m5')),(8,('m3','m4'))):
  transformed=substitute(core[seed],{pair[0]:V(pair[1]),pair[1]:V(pair[0])})
  if (transformed-core[seed]).terms:raise ArithmeticError('Second-ray mass reflection failed')
 payload=dict(check='D2.6.14.global_mass_rays_and_complex_shapes',status='PASS',
     input_sha256=hashlib.sha256(inputs.read_bytes()).hexdigest(),
     global_mass_covers=global_covers,root_records=records,mass_symmetry=symmetry,
     proof='Exact Bernstein exclusion outside four mass holes; uniform Rouche disks exhaust shape roots; Lambda excludes all but a conjugate pair; Krawczyk proves one real mass/shape root per hole orbit; exact mass permutations transport it.')
 output=root/'certificates'/'d2_mass_rays.json';output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
 return [dict(check=payload['check'],status='PASS',items=12,
  detail=dict(certificate=output.name,positive_rays_L1=2,positive_rays_L8=2,
      shape_pair_L1='exactly one nonreal conjugate pair',mass_family_permutations=12,
      scope='Lemma 6.14 mass and common-shape classification; second-order response is a separate obligation'))],{'d2_mass_rays':payload}
