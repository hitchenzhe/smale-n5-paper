"""Table 5 census_model: complete macro branch index, not an LP solver.

The indices are independently enumerated. LP constraints and certificates
must still be added; none of these branch counts asserts LP feasibility.
"""
from itertools import combinations,combinations_with_replacement,product
from pathlib import Path
import json

FAR=('monopole','dipole-z','dipole-w','dipole-equal-F','dipole-equal-G','dipole-equal')
NEAR=('dominant-coincidence','coincident','coincident-Phi','coincident-Psi','coincident-both','inner','dominant-central')
REGIMES=FAR+tuple('z-'+v for v in NEAR)+tuple('w-'+v for v in NEAR)


def mirror(regime):
 if regime.startswith('z-'):return 'w-'+regime[2:]
 if regime.startswith('w-'):return 'z-'+regime[2:]
 return {'dipole-z':'dipole-w','dipole-w':'dipole-z','dipole-equal-F':'dipole-equal-G',
         'dipole-equal-G':'dipole-equal-F'}.get(regime,regime)


def key(pair):return tuple(sorted(pair,key=REGIMES.index))
def canonical(pair):return min(key(pair),key(tuple(map(mirror,pair))),key=lambda p:tuple(REGIMES.index(v) for v in p))


def canonical_case(pair,pattern):
    swapped=pattern[0]+pattern[2]+pattern[1]
    mirrored=tuple(map(mirror,pair))
    target=canonical(pair)
    # The census lists all 27 ordered patterns even for an equal-regime pair.
    # Swap pattern slots only when choosing an orientation of a mixed pair.
    if tuple(pair)==target or mirrored==target:return target,pattern
    return target,swapped


def indices():
 remaining={key(pair) for pair in combinations_with_replacement(REGIMES,2)};orbits=[]
 while remaining:
  seed=min(remaining);images={seed,key(tuple(map(mirror,seed)))};remaining-=images
  orbits.append(dict(regimes=list(canonical(seed)),body_unordered_color_images=[list(p) for p in sorted(images)]))
 if len(REGIMES)!=20 or len(orbits)!=111 or sum(len(o['body_unordered_color_images'])==1 for o in orbits)!=12:
  raise ArithmeticError('D5 regime orbit index mismatch')
 patterns=[''.join(p) for p in product('CRE',repeat=3)]
 survivors={}
 def add(pair,pattern):survivors[canonical_case(pair,pattern)]=None
 far=[('monopole','monopole'),*[('monopole',x) for x in ('dipole-z','dipole-equal','dipole-equal-F')],
      ('dipole-z','dipole-z'),*[('dipole-equal',x) for x in ('dipole-equal','dipole-equal-F')],
      *[('dipole-equal-F',x) for x in ('dipole-equal-F','dipole-equal-G')]]
 for pair in far:add(pair,'RRR')
 coincident=['z-coincident','z-coincident-Phi','z-coincident-Psi']
 for pair in combinations_with_replacement(coincident,2):add(pair,'RRR')
 add(('z-inner','z-inner'),'RRR')
 for x in coincident:add(('z-inner',x),'RCR')
 for pattern in ('CEE','CRR'):add(('z-dominant-coincidence','z-dominant-coincidence'),pattern)
 for pattern in ('RRR','RRC','RCR','CEE','CRR'):add(('z-dominant-central','z-dominant-central'),pattern)
 if len(survivors)!=26:raise ArithmeticError('Restored D5 table does not give 26 pair-pattern entries')
 cases=[]
 for pair,pattern in sorted(survivors):
  for omitted in combinations(range(1,6),2):
   complement=[j for j in range(1,6) if j not in omitted]
   edges=[omitted,*combinations(complement,2)]
   cases.append(dict(regimes=list(pair),pattern=pattern,four_product=list(omitted),
      edges=[list(edge) for edge in edges],status='LP_NOT_YET_CERTIFIED'))
 if len(cases)!=260 or any(len(case['edges'])!=4 for case in cases):raise ArithmeticError('D5 four-product index mismatch')
 return dict(regimes=list(REGIMES),regime_orbits=orbits,patterns=patterns,
             pair_pattern_count=111*27,survivor_table=[dict(regimes=list(pair),pattern=pat) for pair,pat in sorted(survivors)],
             four_product_cases=cases,scope='indices only; table survivors are manuscript declarations, not independently solved LP results')


def verify_branch_index():
 data=indices();path=Path(__file__).resolve().parent/'certificates'/'d5_branch_index.json'
 path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
 return [dict(check='D5.branch_index_111_27_260',status='PASS',items=2997,
   detail=dict(regime_orbits=111,patterns=27,declared_surviving_pair_patterns=26,four_product_cases=260,
    certificate=path.name,source_correction='two inner rows restored from latex screenshots d5/8.png',
    scope='complete macro index; no LP feasibility or branch-completeness proof yet'))],{}


if __name__=='__main__':
 from smale_verify.runner import program_cli
 raise SystemExit(program_cli('ak_diagram5_census_model'))
