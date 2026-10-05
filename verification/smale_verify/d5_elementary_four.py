"""D5 four-product exclusions that already follow from frame/C/M/pattern.

This sufficient submodel proves 194 of the 260 source cases. Six inner cases
remain feasible in this submodel; 60 dominant cases require the remaining
R/T/B/V/M rows. No full-census PASS is inferred from a relaxed LP witness.
"""
from fractions import Fraction as Q
from itertools import combinations,product
from pathlib import Path
import json
from .rational_lp import solve
from verify_lp_certificates import check_case

def model(pair,pattern,edges,cb):
    # affine vectors [a12,b12,a1,b1,a2,b2,constant]
    n=6;basis=[tuple(Q(int(i==j))for i in range(7))for j in range(6)]
    one=(Q(0),)*6+(Q(1),)
    def add(*vs):return tuple(sum(t)for t in zip(*vs))
    def neg(v):return tuple(-q for q in v)
    def constant(q):return tuple(q*v for v in one)
    A=[];b=[];strict=[];labels=[]
    def le(v,w,label,flag=False):
        expr=add(v,neg(w));A.append(list(expr[:-1]));b.append(-expr[-1]);labels.append(label)
        if flag:strict.append(len(A)-1)
    def eq(v,w,label):le(v,w,label+' <=');le(w,v,label+' >=')
    a12,b12,a1,b1,a2,b2=basis
    for name,var in (('a12',a12),('b12',b12)):
        le(constant(-2),var,'frame '+name+' >=-2');le(var,constant(2),'frame '+name+' <=2')
    for index,(regime,aa,bb)in enumerate(zip(pair,(a1,a2),(b1,b2)),1):
        le(constant(-2),aa,'no circle a'+str(index),True);le(constant(-2),bb,'no circle b'+str(index),True)
        if regime.startswith('z-'):
            le(bb,one,'near beta<1',True)
            if regime=='z-inner':le(one,aa,'inner alpha>1',True)
            else:eq(aa,one,'coincident alpha=1')
        else:
            le(aa,one,'far alpha<1',True);le(bb,one,'far beta<1',True)
    if pair==('z-dominant-coincidence','z-dominant-coincidence'):
        le(one,a12,'M same special body a12>1',True)
    for vars_,chosen in zip(((a12,a1,a2),(b12,b1,b2)),cb):
        i,j=chosen;k=3-i-j
        eq(vars_[i],vars_[j],'C minimal pair '+str(chosen));le(vars_[i],vars_[k],'C pair <= third')
    for letter,aa,bb in zip(pattern,(a12,a1,a2),(b12,b1,b2)):
        total=add(aa,bb)
        if letter=='R':eq(total,constant(0),'pattern R')
        elif letter=='C':le(constant(0),total,'pattern C',True)
        elif letter=='E':le(total,constant(0),'pattern E',True)
        else:raise ValueError('Unknown pattern')
    edge_orders=[]
    for i,j in edges:
        if i>=3:order=constant(2)
        elif j==2:order=add(a12,b12)
        else:
            regime=pair[i-1];aa=(a1,a2)[i-1];bb=(b1,b2)[i-1]
            order=add(one,bb)if regime.startswith('z-')else add(aa,bb)
        edge_orders.append(order)
    le(add(*edge_orders),constant(0),'negative four-product (or certified lower bound)',True)
    return dict(variables=n,A=[[str(v)for v in row]for row in A],b=[str(v)for v in b],strict=strict),labels

def verify():
    from ak_diagram5_census_model import indices
    indexed=indices();certificates=[];counts={'excluded':0,'relaxed_feasible':0,'requires_full_model':0}
    for case in indexed['four_product_cases']:
        pair=tuple(case['regimes']);pattern=case['pattern']
        if any('dominant-central'in r for r in pair)or(pair==('z-dominant-coincidence','z-dominant-coincidence')and pattern=='CEE'):
            counts['requires_full_model']+=1;continue
        decisions=[];feasible=False
        for cb in product(combinations(range(3),2),repeat=2):
            system,labels=model(pair,pattern,case['edges'],cb)
            result=solve(system['A'],system['b'],system['strict'])
            check_case(result)
            if result['kind']=='feasible':feasible=True
            decisions.append(dict(C_minimum_pairs=[list(p)for p in cb],row_labels=labels,certificate=result))
        if feasible:
            if pair!=('z-inner','z-inner')or case['four_product'][0]>=3 or case['four_product']==[1,2]:
                raise ArithmeticError('Unexpected relaxed four-product survivor')
            counts['relaxed_feasible']+=1
        else:counts['excluded']+=1
        certificates.append(dict(index=case,status='EXCLUDED_BY_SUBMODEL'if not feasible else'RELAXED_FEASIBLE_NOT_FULL_CERTIFIED',branches=decisions))
    if counts!={'excluded':194,'relaxed_feasible':6,'requires_full_model':60}:raise ArithmeticError('Elementary four-product coverage changed')
    payload=dict(check='D5.8.21.elementary_194_four_product_exclusions',status='PASS',counts=counts,cases=certificates,
        coverage='All nine tropical C minimum-pair disjuncts checked per case. Negative-order exclusions follow from a necessary subset of the manuscript rows. Feasible relaxed witnesses do not prove full-system feasibility.',
        covered_by_other_stages='d5_census checks the full 111x27 census; d5_dominant excludes the 60 remaining dominant products; d5_four_complete checks six full inner witnesses, and d5_source_audit supplies their leading-coefficient exclusion.')
    path=Path(__file__).resolve().parents[1]/'certificates'/'d5_elementary_four.json'
    path.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=194,detail=dict(certificate=path.name,LP_systems=len(certificates)*9,complete_C_disjunction=True,scope='194 sufficient subset exclusions; the separate census/dominant/full-inner stages complete the remaining coverage'))],{}
