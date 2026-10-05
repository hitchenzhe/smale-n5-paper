"""Necessary valuation constraints reconstructed from D5 rows.

This development model is intentionally not registered as a full census:
the exact near-source coefficient audit and central V rows are separate.
"""
from fractions import Fraction as Q
from itertools import combinations
from .rational_simplex import solve

NAMES=('a12','b12','a1','b1','a2','b2','f1','q1','f2','q2','e1','e2','g','t1','t2','B1','B2','v11','v12','v21','v22','v','h','specialB1','specialB2')
N=len(NAMES)
def affine(name):
    if isinstance(name,str):return tuple(Q(int(j==NAMES.index(name)))for j in range(N))+(Q(0),)
    return (Q(0),)*N+(Q(name),)
def add(*vs):return tuple(sum(t)for t in zip(*vs))
def scale(v,k):return tuple(Q(k)*q for q in v)
def sub(u,v):return add(u,scale(v,-1))
V={name:affine(name)for name in NAMES};ONE=affine(1);ZERO=affine(0)
def le(u,v,label,strict=False):return (sub(u,v),strict,label)
def eq(u,v,label):return [le(u,v,label+' <= '),le(v,u,label+' >= ')]
def lp(rows):
    return solve([r[0][:-1]for r in rows],[-r[0][-1]for r in rows],[i for i,r in enumerate(rows)if r[1]])

def tropical(values,label):
    for i,j in combinations(range(len(values)),2):
        # Put the all-tied case in the first minimum pair only. Later
        # pairs use a strict comparison; this is an exhaustive partition
        # for the three-term identities used by this model.
        yield eq(values[i],values[j],label+' tied')+[le(values[i],v,label+' minimal',len(values)==3 and(i,j)!=(0,1))for k,v in enumerate(values)if k not in(i,j)]

def sum_order_zero(values,label):
    for i,value in enumerate(values):
        yield eq(value,ZERO,label+' zero minimum')+[le(ZERO,v,label+' others nonnegative')for k,v in enumerate(values)if k!=i]
    for i,j in combinations(range(len(values)),2):
        yield eq(values[i],values[j],label+' negative tied')+[le(values[i],ZERO,label+' negative minimum',True)]+[le(values[i],v,label+' others above')for k,v in enumerate(values)if k not in(i,j)]

def row_plane(inputs,outputs,label):
    for i,minimum in enumerate(inputs):
        yield [le(minimum,v,label+' unique input minimum',True)for k,v in enumerate(inputs)if k!=i]+sum((eq(out,minimum,label+' no cancellation')for out in outputs),[])
    for i,j in combinations(range(3),2):
        minimum=inputs[i];k=3-i-j;base=eq(inputs[i],inputs[j],label+' two input minimum')+[le(minimum,inputs[k],label+' third higher',True)]
        for cancelled in range(3):
            yield base+sum((eq(out,minimum,label+' noncancelled')if h!=cancelled else[le(minimum,out,label+' possibly cancelled',cancelled!=0)]for h,out in enumerate(outputs)),[])
    base=eq(inputs[0],inputs[1],label+' three input minimum')+eq(inputs[0],inputs[2],label+' three input minimum')
    for cancelled in range(3):
        yield base+sum((eq(out,inputs[0],label+' noncancelled')if h!=cancelled else[le(inputs[0],out,label+' possibly cancelled',cancelled!=0)]for h,out in enumerate(outputs)),[])
    yield base+[le(inputs[0],out,label+' all cancelled',True)for out in outputs]

def frame(pair,pattern):
    rows=[];force=[];near=[]
    for name in ('a12','b12'):
        rows+=[le(affine(-2),V[name],'frame lower'),le(V[name],affine(2),'frame upper')]
    rows.append(le(add(V['a12'],V['b12']),affine(2),'12 is not a zw edge',True))
    for j,regime in enumerate(pair,1):
        aa,bb,ff,gg,ee=(V[name+str(j)]for name in ('a','b','f','q','e'))
        rows+=[le(affine(-2),aa,'no circle alpha',True),le(affine(-2),bb,'no circle beta',True)]
        mirror=regime.startswith('w-')
        if mirror:aa,bb,ff,gg=bb,aa,gg,ff
        regime=regime.replace('w-','z-')
        if regime.startswith('z-'):
            rows+=[le(bb,ONE,'near beta<1',True)]
            rows+=eq(aa,ONE,'near alpha=1')if not regime.endswith(('inner','central'))else[le(ONE,aa,'near alpha>1',True)]
            special='dominant'in regime;A=add(ONE,ee)if special else ONE
            if special:rows+=[le(ZERO,ee,'special e>0',True),le(ee,ONE,'frame special z<=2')]
            f=scale(add(A,scale(bb,3)),Q(-1,2));q=scale(add(scale(A,3),bb),Q(-1,2))
            rows+=[le(affine(-2),f,'cross Z force not maximal',True),le(affine(-2),q,'cross W force not maximal',True),le(add(A,bb),affine(2),'cross is not zw edge',True)]
            fzero=regime=='z-coincident-Phi';qzero=regime=='z-coincident-Psi'
            rows+=[le(f,ff,'Phi zero raises F',True)]if fzero else eq(ff,f,'near F order')
            rows+=[le(q,gg,'Psi zero raises G',True)]if qzero else eq(gg,q,'near G order')
            force.append((add(affine(2),q),add(affine(2),f)))
            if regime=='z-coincident-both':rows.append(le(ONE,ZERO,'empty regime'))
            if regime=='z-dominant-central':
                near.append((j,'central',mirror,A,q))
                # Central V/B lower bounds are constructed separately.
                continue
            near.append((j,regime,mirror,A,q))
            v1=add(ONE,q)if special else scale(add(ONE,bb),Q(-1,2))
            v2=add(ONE,f)if special else scale(add(ONE,scale(bb,3)),Q(-1,2))
            B=v1;exact= special or regime in ('z-inner','z-coincident-Phi','z-coincident-Psi')
            vv1,vv2=V['v1'+str(j)],V['v2'+str(j)]
            if mirror:vv1,vv2=vv2,vv1
            rows+=eq(vv1,v1,'near V1 exact')if exact else[le(v1,vv1,'near V1 lower')]
            rows+=[le(v2,vv2,'near V2 lower')]
            rows+=eq(V['t'+str(j)],v1,'near torque exact')if exact else[le(v1,V['t'+str(j)],'near torque lower')]
            Bexact=special or regime in ('z-inner','z-coincident-Phi','z-coincident-Psi')
            rows+=eq(V['B'+str(j)],B,'near B exact')if Bexact else[le(B,V['B'+str(j)],'near B lower')]
        else:
            rows+=[le(aa,ONE,'far alpha<1',True),le(bb,ONE,'far beta<1',True)]
            f=scale(add(aa,scale(bb,3)),Q(-1,2));q=scale(add(scale(aa,3),bb),Q(-1,2))
            rawf,rawq=f,q
            rows+=[le(affine(-2),rawf,'cross Z force not maximal',True),le(affine(-2),rawq,'cross W force not maximal',True)]
            v1=add(ONE,q);v2=add(ONE,f);baseB=scale(add(aa,bb),Q(-1,2))
            if regime=='monopole':
                rows+=eq(ff,f,'monopole F')+eq(gg,q,'monopole G')+eq(V['B'+str(j)],baseB,'monopole B')
                rows+=[le(v1,V['v1'+str(j)],'monopole V1 lower'),le(v2,V['v2'+str(j)],'monopole V2 lower')]
                # T lower bound is min(v1,v2), deferred as a disjunction.
            else:
                if regime=='dipole-z':rows.append(le(bb,aa,'dipole-z',True));f=add(f,sub(ONE,aa));q=add(q,sub(ONE,aa));border=add(baseB,sub(ONE,aa))
                elif regime=='dipole-w':rows.append(le(aa,bb,'dipole-w',True));f=add(f,sub(ONE,bb));q=add(q,sub(ONE,bb));border=add(baseB,sub(ONE,bb))
                else:rows+=eq(aa,bb,'equal dipole');f=add(f,sub(ONE,aa));q=add(q,sub(ONE,aa));border=add(baseB,sub(ONE,aa))
                rows+=[le(f,ff,'dipole F zero',True)]if regime=='dipole-equal-F'else eq(ff,f,'dipole F')
                rows+=[le(q,gg,'dipole G zero',True)]if regime=='dipole-equal-G'else eq(gg,q,'dipole G')
                rows+=eq(V['B'+str(j)],border,'dipole B')if regime!='dipole-equal'else[le(border,V['B'+str(j)],'dipole B lower')]
                rows+=eq(V['v1'+str(j)],v1,'dipole V1')+eq(V['v2'+str(j)],v2,'dipole V2')
            force.append((add(affine(2),rawq),add(affine(2),rawf)))
    for letter,aa,bb in zip(pattern,(V['a12'],V['a1'],V['a2']),(V['b12'],V['b1'],V['b2'])):
        r=add(aa,bb)
        rows+=eq(r,ZERO,'pattern R')if letter=='R'else[le(ZERO,r,'pattern C',True)]if letter=='C'else[le(r,ZERO,'pattern E',True)]
    if len(near)==2 and near[0][2]==near[1][2]:
        kinds=[rec[1]for rec in near]
        rows+=eq(near[0][4],near[1][4],'uniform forcing requires tied leading sources')
        if kinds==['z-dominant-coincidence']*2:
            rows.append(le(ONE,V['b12'if near[0][2]else'a12'],'M same special difference',True))
            rows+=eq(near[0][4],near[1][4],'M forcing tie')
        if 'z-dominant-coincidence'in kinds and any(k in('z-coincident-Phi','z-coincident-Psi')for k in kinds):rows.append(le(ONE,ZERO,'M Phi/Psi excluded'))
        if 'central'in kinds and 'z-inner'in kinds:rows.append(le(ONE,ZERO,'slice inner empty'))
        if kinds==['central']*2:
            rows+=eq(near[0][4],near[1][4],'central forcing tie')+[le(add(affine(3),near[0][4]),V['g'],'central g>3+phi',True)]
    if any(rec[1]=='central'for rec in near):
        if any(regime.startswith('dipole')or regime.endswith('inner')for regime in pair):rows.append(le(ONE,ZERO,'central slice excludes dipoles and inner bodies'))
    r2=[add(V['a12'],V['b12']),add(V['a1'],V['b1']),add(V['a2'],V['b2'])]
    negative=[i for i,ch in enumerate(pattern)if ch=='E']
    Agroups=[]
    if len(negative)==1 or pattern=='CCC':rows.append(le(ONE,ZERO,'A impossible unique negative or all positive'))
    elif len(negative)==2:rows+=eq(r2[negative[0]],r2[negative[1]],'A two negative leading orders tie')
    elif len(negative)==3:Agroups=[list(tropical(r2,'A negative minimum'))]
    groups=[list(row_plane([scale(add(V['a12'],scale(V['b12'],3)),Q(-1,2)),V['f1'],V['f2']],[V['a12'],V['a1'],V['a2']],'Rz')),
            list(row_plane([scale(add(scale(V['a12'],3),V['b12']),Q(-1,2)),V['q1'],V['q2']],[V['b12'],V['b1'],V['b2']],'Rw')),
            list(tropical([V['a12'],V['a1'],V['a2']],'Cz')),list(tropical([V['b12'],V['b1'],V['b2']],'Cw'))]+Agroups
    if len(near)==2 and all(rec[1]=='central'for rec in near)and near[0][2]==near[1][2]:
        rows.append(le(ONE,V['g'],'central offset g>1',True))
        aa=[V['b'+str(j)]if near[0][2]else V['a'+str(j)]for j in (1,2)]
        bb=[V['a'+str(j)]if near[0][2]else V['b'+str(j)]for j in (1,2)]
        A=[rec[3]for rec in near]
        groups=[list(tropical([aa[0],A[0],V['g']],'central source 1')),
                list(tropical([aa[1],A[1],V['g']],'central source 2')),
                list(tropical([A[0],A[1],V['b12'if near[0][2]else'a12']],'central special C'))]+groups
        groups.append([[le(aa[0],V['g'],'M not deep 1')],[le(aa[1],V['g'],'M not deep 2')],
            [le(A[0],V['g'],'M A1<g',True)],[le(V['g'],A[0],'M g<A1',True)],
            [le(A[1],V['g'],'M A2<g',True)],[le(V['g'],A[1],'M g<A2',True)],
            [le(scale(V['g'],-1),bb[0],'M beta1>=-g')],[le(scale(V['g'],-1),bb[1],'M beta2>=-g')]])
    return rows,groups,force,near

def explore(rows,groups,limit=100000,stop_at_first=False):
    leaves=[];visited=0
    def visit(rows,depth,path):
        nonlocal visited
        visited+=1
        if visited>limit:raise RuntimeError('Valuation tree budget exhausted; no completeness result')
        result=lp(rows)
        if result['kind']=='infeasible':return False
        if depth==len(groups):leaves.append(dict(path=path,certificate=result));return True
        any_feasible=False
        for index,extra in enumerate(groups[depth]):
            any_feasible=visit(rows+extra,depth+1,path+[index])or any_feasible
            if any_feasible and stop_at_first:return True
        return any_feasible
    visit(rows,0,[]);return leaves,visited

def add_virial(rows,groups,force,near):
    rows=rows+eq(V['t1'],V['t2'],'T equal torque orders')
    central=[rec for rec in near if rec[1]=='central']
    if central and len({rec[2]for rec in central})==1:
        rows=rows+[le(ONE,V['g'],'central z offset>1',True),le(ONE,V['h'],'central w offset>1',True)]
        for j,kind,mirror,A,q in central:
            bb=V['a'+str(j)]if mirror else V['b'+str(j)]
            f=scale(add(A,scale(bb,3)),Q(-1,2));generic1=scale(add(ONE,bb),Q(-1,2));generic2=scale(add(ONE,scale(bb,3)),Q(-1,2))
            vv1,vv2=V['v1'+str(j)],V['v2'+str(j)]
            if mirror:vv1,vv2=vv2,vv1
            groups=groups+[list(tropical([vv1,add(V['g'],q),generic1],'central V1 share')),
                           list(tropical([vv2,add(V['h'],f),generic2],'central V2 share'))]
            special=V['specialB'+str(j)];d=scale(add(A,bb),Q(-1,2))
            groups=groups+[[[le(V['g'],A,'central B g<A',True)]+eq(special,add(V['g'],q),'central B special dipole'),
                             [le(A,V['g'],'central B A<g',True)]+eq(special,d,'central B special monopole'),
                             eq(A,V['g'],'central B g=A')+[le(d,special,'central B coincident lower')]],
                           list(tropical([V['B'+str(j)],special,generic1],'central B share'))]
    groups=groups+[list(sum_order_zero([V['B1'],V['B2'],scale(add(V['a12'],V['b12']),Q(-1,2))],'B')),
                  list(tropical([V['v11'],V['v21'],V['t1']],'source 1 torque difference')),
                  list(tropical([V['v12'],V['v22'],V['t2']],'source 2 torque difference')),
                  list(tropical([V['v11'],V['v12'],V['v']],'V1 sum')),
                  list(tropical([V['v21'],V['v22'],V['v']],'V2 sum'))]
    lower=[affine(2)]+[sub(scale(f,2),ONE)for pair in force for f in pair]
    bound=[[le(form,V['v'],'V Morse-Bott lower')]for form in lower]
    if len(near)==2 and all(rec[1]=='z-dominant-coincidence'for rec in near):
        # Both leading crossed forcings tie and cancel on their common
        # special body; all other crossed and uncrossed terms have strictly
        # higher order. Hence kappa>2+phi rather than merely >=2+phi.
        bound=[[le(affine(2),V['v'],'V reaches self order')],
               [le(add(affine(3),scale(near[0][4],2)),V['v'],'V raised cancelling forcing',True)]]
    return rows,groups+[bound]

def negative_product(pair,edges,special=3):
    values=[]
    for i,j in edges:
        if i>=3:values.append(affine(2))
        elif j==2:values.append(add(V['a12'],V['b12']))
        else:
            regime=pair[i-1];aa=V['a'+str(i)];bb=V['b'+str(i)]
            if regime.startswith('w-'):aa,bb=bb,aa
            if regime.startswith(('z-','w-')):aa=add(ONE,V['e'+str(i)])if 'dominant'in regime and j==special else ONE
            values.append(add(aa,bb))
    return le(add(*values),ZERO,'negative four-product',True)
