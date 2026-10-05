"""Exact conditional ray orders for Remark 11.10, reconstructed from the chart.

This computes orders on a hypothetical solution ray with finite chart
coordinates and a nondegenerate lower wing. It does not construct a germ;
the main Diagram-8 theorem excludes such germs separately.
"""
from fractions import Fraction as Q
from itertools import combinations
from .exact import P,V,equal,zero,substitute
from ak_diagram8_relation40_derivation import chart_rows

p=V('ray_p')
UPPER=((1,2),(1,5),(2,5))
LOWER=((3,4),(3,5),(4,5))
CROSSED=((1,3),(1,4),(2,3),(2,4))
PAIRS=tuple(combinations(range(1,6),2))

def affine_coefficients(poly):
    poly=P.coerce(poly)
    if poly.names()-{'ray_p'} or any(k not in (0,1)for k in poly.exponents('ray_p')):
        raise ValueError('An order must be affine in p')
    return poly.terms.get((),Q(0)),poly.terms.get((('ray_p',1),),Q(0))

def encode(poly):
    a,b=affine_coefficients(poly)
    return dict(constant=str(a),p_coefficient=str(b))

def evaluate(poly,point):
    a,b=affine_coefficients(poly)
    return a+b*Q(point)

def open_positive(poly,lo=-4,hi=2):
    """A nonzero affine function nonnegative at both ends is positive inside."""
    left,right=evaluate(poly,lo),evaluate(poly,hi)
    if left<0 or right<0 or(left==0 and right==0):
        raise ArithmeticError('Strict affine inequality fails on the open interval')
    return dict(polynomial=encode(poly),domain=[str(lo),str(hi)],
                endpoint_values=[str(left),str(right)],endpoints='excluded')

def halfline_positive(poly,endpoint,direction):
    a,b=affine_coefficients(poly);v=evaluate(poly,endpoint)
    if v<0 or b*direction<0 or(v==0 and b==0):
        raise ArithmeticError('Strict affine inequality fails on the half-line')
    return dict(polynomial=encode(poly),endpoint=str(endpoint),endpoint_value=str(v),
                direction=direction,endpoint_included=False)

def support(poly):
    return sorted({(dict(mon).get('a8',0),dict(mon).get('eta8',0))for mon in poly.terms})

def corner_weight(exponents):
    i,j=exponents
    return i*(1-p/2)+j*(2+p/2)

def reconstruction():
    masses={i:V('m'+str(i))for i in range(1,6)}
    data=chart_rows(masses);x=data['x'];z,w,delta=(data[k]for k in('z','w','deltas'))
    if set(delta)!=set(PAIRS):raise ArithmeticError('Incomplete distance index')
    for j,e in enumerate(PAIRS,10):
        i,k=e
        equal(data['E'][j],delta[e]**2*(z[k]-z[i])*(w[k]-w[i])-1,
              'D8 source distance tie '+str(e))
    equal(z[2]-z[1],x['A2']-x['A1'],'D8 exact upper separation')
    equal(w[4]-w[3],data['nu']**-1*(x['B4']-x['B3']),
          'D8 exact lower separation')
    # Positive-mass leading coordinates, using the two gauges and centres.
    M12,M34=V('ray_M12'),V('ray_M34')
    point={'d8_A1':-masses[2]/M12,'d8_A2':masses[1]/M12,'d8_C':0,
           'd8_B3':-masses[4]/M34,'d8_B4':masses[3]/M34,'d8_D':0}
    identities=[masses[1]*x['A1']+masses[2]*x['A2'],
                masses[3]*x['B3']+masses[4]*x['B4'],
                x['A2']-x['A1']-1,x['B4']-x['B3']-1]
    for polynomial in identities:
        polynomial=substitute(polynomial,point).clear('ray_M12').clear('ray_M34')
        zero(substitute(polynomial,{'ray_M12':masses[1]+masses[2],
                                    'ray_M34':masses[3]+masses[4]}),
             'D8 positive leading centres and gauges')
    C_support=support(data['E'][4]-x['C'])
    D_support=support(data['E'][1]-x['D'])
    if C_support!=[(2,2),(3,3),(5,5)]:
        raise ArithmeticError('D8 C-order forcing support')
    if D_support!=[(0,2),(3,3),(6,2),(9,5)]:
        raise ArithmeticError('D8 D-order forcing support')
    C_weights=[corner_weight(e)for e in C_support]
    D_weights=[corner_weight(e)for e in D_support]
    if [encode(v)for v in C_weights]!=[encode(P.coerce(v))for v in(6,9,15)]:
        raise ArithmeticError('D8 C order is not bounded below by six')
    D_gaps=[open_positive(v-(4+p))for v in D_weights[1:]]
    # The endpoint p=2 is also covered by nonnegative D-weight gaps.
    if any(evaluate(v-(4+p),2)<0 for v in D_weights):
        raise ArithmeticError('D8 D-order bound misses the upper endpoint')
    distances={};separations={};components={}
    for e in PAIRS:
        powers=support(delta[e])
        if len(powers)!=1:raise ArithmeticError('Distance is not a single corner monomial')
        inverse_order=corner_weight(powers[0])
        expected=p-2 if e in UPPER else P.coerce(3 if e in CROSSED else 0)
        equal(inverse_order,expected,'D8 inverse-distance order '+str(e))
        squared=-2*inverse_order
        if e in UPPER:
            zo=P();wo=squared
        elif e in LOWER:
            wo=P.coerce(-6);zo=squared-wo
        else:
            zo=P();wo=P.coerce(-6)
        equal(zo+wo,squared,'D8 tie-order sum '+str(e))
        distances[e]=squared;separations[e]=(zo,wo)
        components[e]=(zo+3*inverse_order,wo+3*inverse_order)
    four_products={};product_edges={}
    for e in PAIRS:
        other=[i for i in range(1,6)if i not in e]
        edges=(e,*combinations(other,2))
        product_edges[e]=edges
        four_products[e]=sum((distances[q]for q in edges),P())
    least=-8-2*p;minimum_gaps=[]
    for e,value in four_products.items():
        wanted=4-2*p if e==(1,2)else 12-6*p if e==(3,4)else least
        equal(value,wanted,'D8 four-product order '+str(e))
        if value.terms!=least.terms:minimum_gaps.append(open_positive(value-least))
    negative_minimum=open_positive(-least)
    provenance=dict(C_forcing_support=[list(e)for e in C_support],
        C_forcing_orders=[encode(v)for v in C_weights],C_lower_bound='6',
        D_forcing_support=[list(e)for e in D_support],
        D_forcing_orders=[encode(v)for v in D_weights],D_lower_bound=encode(4+p),
        D_weight_gap_certificates=D_gaps,
        leading_coordinates=dict(A1='-m2/(m1+m2)',A2='m1/(m1+m2)',
            B3='-m4/(m3+m4)',B4='m3/(m3+m4)',C='0',D='0'),
        positive_units=['m1','m2','m3','m4','m5','m1+m2','m3+m4'],
        nondegeneracy='All inverse-distance chart coordinates have nonzero finite limits; both gauge separations are one. Upper z and lower w wing coefficients, and crossed coefficients, are nonzero.',
        assumptions='Hypothetical ratio-chart solution ray with finite coordinates and finite nondegenerate lower wing; -4<p<2, or p=2 where indicated. No existence assertion.')
    return dict(data=data,distances=distances,separations=separations,components=components,
                four_products=four_products,product_edges=product_edges,
                minimum_gaps=minimum_gaps,negative_minimum=negative_minimum,
                provenance=provenance)

def diagram8_gaps(result):
    c=result['components'];proofs=[]
    zmin=3*p-6;wmin=P.coerce(-6)
    for e in PAIRS:
        zz,ww=c[e]
        if e in UPPER:equal(zz,zmin,'D8 upper z stroke')
        else:proofs.append(dict(kind='z edge gap',edge=list(e),proof=open_positive(zz-zmin)))
        if e in LOWER:equal(ww,wmin,'D8 lower w stroke')
        else:proofs.append(dict(kind='w edge gap',edge=list(e),proof=open_positive(ww-wmin)))
    # Lower bounds suffice when a position has a cancelling coefficient.
    for name,gap in [('upper z positions',-zmin),('lower z positions',6-zmin),
                     ('upper w D contribution',p-2-wmin),
                     ('upper w P/Q contribution',4-2*p-wmin)]:
        proofs.append(dict(kind=name,proof=open_positive(gap)))
    return proofs
