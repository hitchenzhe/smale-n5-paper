"""Complete 60 dominant four-product exclusions with replayable LP trees.

Only a necessary subset (frame, R, C, A and uniform-forcing consequences)
is used. This proves exclusions without certifying the full census model.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
from .d5_valuation import frame,negative_product,lp,N,NAMES

ROOT=Path(__file__).resolve().parents[1]

def reconstruct(case):
    pair=tuple(case['regimes']);rows,groups,force,near=frame(pair,case['pattern'])
    rows.append(negative_product(pair,case['edges']))
    return rows,groups

def build_tree(rows,groups,depth=0):
    result=lp(rows)
    if result['kind']=='infeasible':return dict(motzkin=[[i,q]for i,q in enumerate(result['y'])if Q(q)])
    if depth==len(groups):raise ArithmeticError('Dominant four-product remains feasible; no proof claimed')
    return dict(children=[build_tree(rows+extra,groups,depth+1)for extra in groups[depth]])

def replay(tree,rows,groups,depth=0):
    if set(tree)=={'motzkin'}:
        total=[Q(0)]*(N+1);strict=False;seen=set()
        for index,value in tree['motzkin']:
            if type(index)is not int or index in seen or not 0<=index<len(rows):raise ValueError('Invalid Motzkin row index')
            seen.add(index);q=Q(value)
            if q<=0:raise ValueError('Motzkin coefficients must be strictly positive when present')
            row,flag,label=rows[index];total=[v+q*w for v,w in zip(total,row)];strict|=flag
        # Rows are affine <=0, so a strictly positive constant is a
        # contradiction; a zero constant needs a strict row multiplier.
        if any(total[:-1])or total[-1]<0 or(total[-1]==0 and not strict):raise ArithmeticError('Motzkin contradiction failed')
        return 1
    if set(tree)!={'children'}or depth>=len(groups)or len(tree['children'])!=len(groups[depth]):raise ValueError('LP branch tree is incomplete')
    return sum(replay(child,rows+extra,groups,depth+1)for child,extra in zip(tree['children'],groups[depth]))

def cases():
    from ak_diagram5_census_model import indices
    return [r for r in indices()['four_product_cases']if any('dominant-central'in s for s in r['regimes'])or(tuple(r['regimes'])==('z-dominant-coincidence',)*2 and r['pattern']=='CEE')]

def generate():
    records=[]
    for case in cases():
        rows,groups=reconstruct(case);tree=build_tree(rows,groups);leaves=replay(tree,rows,groups)
        records.append(dict(case=case,tree=tree,leaves=leaves));print('D5 dominant',case['pattern'],case['four_product'],'Motzkin leaves',leaves,flush=True)
    path=ROOT/'certificates/d5_dominant_inputs.json';path.write_text(json.dumps(dict(variables=list(NAMES),cases=records),indent=2)+'\n',encoding='utf-8')
    return path

def verify():
    path=ROOT/'certificates/d5_dominant_inputs.json';data=json.loads(path.read_text(encoding='utf-8'));expected=cases()
    if data['variables']!=list(NAMES)or len(data['cases'])!=60 or len(expected)!=60:raise ValueError('Dominant case index mismatch')
    total=0
    for record,case in zip(data['cases'],expected):
        if record['case']!=case:raise ValueError('Dominant certificate case differs from incidence table')
        rows,groups=reconstruct(case);count=replay(record['tree'],rows,groups)
        if count!=record['leaves']:raise ValueError('Dominant proof leaf count mismatch')
        total+=count
    payload=dict(check='D5.8.21.all_60_dominant_four_product_exclusions',status='PASS',cases=60,motzkin_leaves=total,
        certificate=path.name,scope='Complete exclusions for all 60 previously outstanding dominant four-product cases, using a necessary subset. This does not assert full regime-census feasibility.',
        model_proof=['Exact single-special-edge forcing orders f=-(1+e)/2-3beta/2 and q=-3(1+e)/2-beta/2; 0<e<=1 follows from the frame.',
            'Uniform forcing and the single-source nonuniformity lemma require the two raw leading crossed source orders to tie.',
            'For a common central body, shared g>1, g>3+q, both centre/special tropical relations, the special pair relation, and the excluded deep conjunction are precisely the proved consequences of Lemma 8.19.',
            'R has all nonzero coefficients, rank two on every coordinate plane and common one-dimensional kernel. Sixteen input/output valuation alternatives exhaust unique, double and triple input minima; C and A are separately imposed.',
            'Each proof-tree child represents every alternative of its deterministic disjunction. Exact nonnegative Motzkin combinations close every leaf, with a positive strict-row multiplier when the combined constant is zero.',
            'For dominant coincidence the special body is normalized to 3. Permutations of the triangle bodies permute the complete list of ten four-products, so the ten checked exclusions cover every special placement.'])
    output=ROOT/'certificates/d5_dominant.json';output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=60,detail=payload)],{}
