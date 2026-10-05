"""Reconstruct and certify the complete macro index with rational LP trees.

The model deliberately keeps conservative lower bounds where a coefficient
has not yet been proved nonzero. Extra feasible cases remain audit failures,
not silently asserted source survivors.
"""
from pathlib import Path
import json
from fractions import Fraction as Q
from .d5_valuation import frame,add_virial,lp,NAMES
from .d5_dominant import replay
from verify_lp_certificates import check_case
ROOT=Path(__file__).resolve().parents[1]

def search(rows,groups,depth=0,path=()):
    result=lp(rows)
    if result['kind']=='infeasible':return dict(motzkin=[[i,q]for i,q in enumerate(result['y'])if Q(q)])
    if depth==len(groups):return dict(feasible=True,path=list(path),x=result['x'])
    children=[]
    for index,extra in enumerate(groups[depth]):
        child=search(rows+extra,groups,depth+1,path+(index,))
        if child.get('feasible'):return child
        children.append(child)
    return dict(children=children)

def model(pair,pattern,full):
    rows,groups,force,near=frame(pair,pattern)
    if full:
        rows,groups=add_virial(rows,groups,force,near)
        if full=='TB':
            excluded=('V1 sum','V2 sum','central V1 share','central V2 share','V Morse-Bott lower','V reaches self order','V raised cancelling forcing')
            groups=[choices for choices in groups if not any(tag in choices[0][0][2]for tag in excluded)]
    return rows,groups

def check_result(result,rows,groups):
    if result.get('feasible'):
        path=result['path']
        if len(path)!=len(groups):raise ValueError('Incomplete feasibility path')
        for index,choices in zip(path,groups):
            if type(index)is not int or not 0<=index<len(choices):raise ValueError('Invalid feasibility branch index')
            rows=rows+choices[index]
        return check_case(dict(variables=len(NAMES),A=[[str(q)for q in row[:-1]]for row,flag,label in rows],b=[str(-row[-1])for row,flag,label in rows],strict=[i for i,(row,flag,label)in enumerate(rows)if flag],kind='feasible',x=result['x']))
    return replay(result,rows,groups)

def generate():
    from ak_diagram5_census_model import indices
    import time,os
    (ROOT/'certificates/d5_census_generator.pid').write_text(str(os.getpid()))
    indexed=indices();declared={(tuple(r['regimes']),r['pattern'])for r in indexed['survivor_table']};records=[];start=time.time();completed=0
    path=ROOT/'certificates/d5_census_inputs.json'
    if path.is_file():
        old=json.loads(path.read_text())
        if old['variables']==list(NAMES):
            records=old['cases'];completed=old['completed_orbits']
            for record in records:
                if record['result'].get('feasible'):
                    rows,groups=model(tuple(record['pair']),record['pattern'],record['full_model']);record['result']=search(rows,groups)
                rows,groups=model(tuple(record['pair']),record['pattern'],record['full_model']);check_result(record['result'],rows,groups)
            print('Resumed',completed,'completed D5 regime orbits; all stored proofs replayed',flush=True)
    for i,orbit in enumerate(indexed['regime_orbits']):
        if i<completed:continue
        pair=tuple(orbit['regimes']);survivors=[]
        for pattern in indexed['patterns']:
            full=False;rows,groups=model(pair,pattern,False);result=search(rows,groups)
            if result.get('feasible'):
                full='TB';rows,groups=model(pair,pattern,full);result=search(rows,groups)
                if result.get('feasible'):
                    full=True;rows,groups=model(pair,pattern,full);result=search(rows,groups)
            check_result(result,rows,groups)
            empty=any('coincident-both'in r for r in pair)
            wanted=(pair,pattern)in declared
            if result.get('feasible'):survivors.append(pattern+('!'if not wanted else''))
            if wanted and not result.get('feasible'):print('AUDIT source-declared survivor is infeasible:',pair,pattern,flush=True)
            records.append(dict(pair=list(pair),pattern=pattern,empty_regime=empty,declared_survivor=wanted,full_model=full,result=result))
        print('D5 macro',i+1,'/111',pair,'survivors',','.join(survivors),'elapsed',round(time.time()-start,1),flush=True)
        path.write_text(json.dumps(dict(variables=list(NAMES),completed_orbits=i+1,cases=records),separators=(',',':'))+'\n',encoding='utf-8')
    return path

def verify():
    from ak_diagram5_census_model import indices
    indexed=indices();path=ROOT/'certificates/d5_census_inputs.json';data=json.loads(path.read_text())
    if data['variables']!=list(NAMES)or data['completed_orbits']!=111 or len(data['cases'])!=2997:raise ValueError('Incomplete 111x27 census certificate')
    declared={(tuple(r['regimes']),r['pattern'])for r in indexed['survivor_table']};feasible=[];extra=[];missing=[];leaves=0
    for record,(orbit,pattern)in zip(data['cases'],((o,p)for o in indexed['regime_orbits']for p in indexed['patterns'])):
        pair=tuple(orbit['regimes'])
        if record['pair']!=list(pair)or record['pattern']!=pattern:raise ValueError('Census index mismatch')
        level=record['full_model']
        if not(type(level)is bool or level=='TB'):
            raise ValueError('Unknown census model level')
        if record['result'].get('feasible') and level is not True:
            raise ValueError('A census survivor requires a full-model witness')
        rows,groups=model(pair,pattern,record['full_model']);value=check_result(record['result'],rows,groups)
        if record['result'].get('feasible'):
            feasible.append([list(pair),pattern])
            if (pair,pattern)not in declared:extra.append([list(pair),pattern])
        else:
            leaves+=value
            if (pair,pattern)in declared:missing.append([list(pair),pattern])
    status='PASS'if not(extra or missing)else'INCOMPLETE'
    payload=dict(check='D5.8.20.full_111_by_27_regime_census',status=status,macro_cases=2997,motzkin_leaves=leaves,feasible=feasible,extra_feasible=extra,missing_survivors=missing,certificate=path.name,
        scope='Every generated row and every exhaustive infeasible branch is replayed. A conservative coefficient bound may leave extra feasible cases; these prevent a source-census PASS.')
    output=ROOT/'certificates/d5_census.json';output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status=status,items=2997,detail=payload)],{}
