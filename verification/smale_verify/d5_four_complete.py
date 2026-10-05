"""Bind all 260 indexed four-product decisions and full inner witnesses."""
from fractions import Fraction as Q
from pathlib import Path
import json
from .d5_valuation import NAMES,V,frame,add_virial,negative_product,eq
from .d5_census import check_result
ROOT=Path(__file__).resolve().parents[1]

def indices():
    from ak_diagram5_census_model import indices as all_indices
    return all_indices()['four_product_cases']

def key(case):return (tuple(case['regimes']),case['pattern'],tuple(case['four_product']))

def inner_system(case):
    rows,groups,force,near=frame(tuple(case['regimes']),case['pattern']);rows,groups=add_virial(rows,groups,force,near)
    rows.append(negative_product(tuple(case['regimes']),case['edges']))
    # The single opposite-sign offset also makes this inner coefficient
    # nonzero; enforce its exact order in these positive witnesses.
    from .d5_valuation import affine
    rows+=eq(V['v21'],affine(Q(25,8)),'inner V2 exact nonzero coefficient')+eq(V['v22'],affine(Q(25,8)),'inner V2 exact nonzero coefficient')
    return rows,groups

def generate():
    d={name:Q(0)for name in NAMES};d.update(a12=Q(7,4),b12=Q(-7,4),v=Q(25,8))
    for j in(1,2):
        for name,value in [('a',Q(7,4)),('b',Q(-7,4)),('f',Q(17,8)),('q',Q(-5,8)),('v1',Q(3,8)),('v2',Q(25,8)),('B',Q(3,8)),('t',Q(3,8))]:d[name+str(j)]=value
    x=[d[name]for name in NAMES];records=[]
    for case in indices():
        if tuple(case['regimes'])!=('z-inner',)*2 or case['pattern']!='RRR' or case['four_product'][0]>=3 or case['four_product']==[1,2]:continue
        rows,groups=inner_system(case);path=[]
        for choices in groups:
            possible=[]
            for i,extra in enumerate(choices):
                if all((sum(q*v for q,v in zip(row[:-1],x))+row[-1]<0)if flag else(sum(q*v for q,v in zip(row[:-1],x))+row[-1]<=0)for row,flag,label in extra):possible.append(i)
            if not possible:raise ArithmeticError('Explicit inner witness misses a valuation disjunction')
            path.append(possible[0])
        result=dict(feasible=True,path=path,x=list(map(str,x)));check_result(result,rows,groups)
        records.append(dict(case=case,result=result,negative_order='-1/4'))
    if len(records)!=6:raise ArithmeticError('Inner exception incidence count')
    path=ROOT/'certificates/d5_inner_witnesses.json';path.write_text(json.dumps(dict(variables=list(NAMES),cases=records),indent=2)+'\n',encoding='utf-8');return path

def verify():
    elementary=json.loads((ROOT/'certificates/d5_elementary_four.json').read_text());dominant=json.loads((ROOT/'certificates/d5_dominant_inputs.json').read_text());inner=json.loads((ROOT/'certificates/d5_inner_witnesses.json').read_text())
    covered={};counts={'infeasible':0,'feasible_inner':0}
    for record in elementary['cases']:
        if record['status']=='EXCLUDED_BY_SUBMODEL':covered[key(record['index'])]='infeasible';counts['infeasible']+=1
    for record in dominant['cases']:
        k=key(record['case'])
        if k in covered:raise ValueError('Duplicate four-product exclusion')
        covered[k]='infeasible';counts['infeasible']+=1
    if inner['variables']!=list(NAMES):raise ValueError('Inner variable order mismatch')
    canonical={key(case):case for case in indices()}
    for record in inner['cases']:
        case=record['case'];k=key(case)
        if k in covered or canonical.get(k)!=case:raise ValueError('Unexpected or duplicate inner exception')
        if not record['result'].get('feasible'):raise ValueError('Inner exception requires a feasible witness')
        rows,groups=inner_system(case);check_result(record['result'],rows,groups)
        covered[k]='feasible_inner';counts['feasible_inner']+=1
    if set(covered)!={key(c)for c in indices()}or counts!={'infeasible':254,'feasible_inner':6}:raise ArithmeticError('Incomplete 260-case four-product cover')
    payload=dict(check='D5.8.21.complete_260_four_product_decisions',status='PASS',counts=counts,
        inner_witness='a12=a1=a2=7/4, b12=b1=b2=-7/4; all six exceptional products have exact order -1/4; every necessary full R/C/T/A/B/V/M row and strict branch inequality is checked',
        depends_on=['D5.8.21.elementary_194_four_product_exclusions','D5.8.21.all_60_dominant_four_product_exclusions','D5.near_coefficients_empty_regimes_R_planes_and_deep_exclusions'],
        scope='Every indexed case has an exclusion proof or a full necessary-model witness. The six inner LP survivors are excluded as germs by the separately checked leading-coefficient argument, not by LP infeasibility.')
    path=ROOT/'certificates/d5_complete_four_products.json';path.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=260,detail=payload)],{}
