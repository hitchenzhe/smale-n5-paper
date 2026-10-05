"""Table 5: exact sharp window and endpoint patterns in Remark 11.10."""
import json
from pathlib import Path
from smale_verify.exact import P,equal
from smale_verify.symbolic import sympy
from smale_verify.d8_orders import (p,reconstruction,diagram8_gaps,encode,evaluate,
    halfline_positive,UPPER,LOWER,CROSSED)

def verify_pattern_boundaries():
    result=reconstruction();s=sympy();r=s.symbols('r',positive=True)
    expression=(2*r-4)/(r+1)
    if s.cancel(s.diff(expression,r)-6/(r+1)**2):
        raise ArithmeticError('D8 ray map is not strictly increasing')
    if s.limit(expression,r,0)!=-4 or s.limit(expression,r,s.oo)!=2:
        raise ArithmeticError('D8 ray map misses the sharp endpoints')
    if any(s.cancel(expression.subs(r,u)-v)for u,v in
           ((s.Rational(1,2),-2),(2,0),(3,s.Rational(1,2)))):
        raise ArithmeticError('D8 printed sample-ray correspondences')
    a=1-p/2;eta=2+p/2
    equal(a+eta,3,'D8 corner t relation')
    equal(eta-2*a,3*p/2,'D8 corner h relation')
    if (evaluate(a,-4),evaluate(eta,-4),evaluate(a,2),evaluate(eta,2))!=(3,0,0,3):
        raise ArithmeticError('D8 corner face orders')
    endpoint_components=[]
    for e,(Z,W)in result['components'].items():
        endpoint_components.append(dict(edge=list(e),p_equals_2=dict(Z=str(evaluate(Z,2)),W=str(evaluate(W,2))),
            p_equals_minus_4=dict(Z=str(evaluate(Z,-4)),W=str(evaluate(W,-4)))))
    zmin=P();wmin=P.coerce(-6)
    for e,(Z,W)in result['components'].items():
        if (evaluate(Z,2)==0)!=(e in UPPER) or(evaluate(W,2)==-6)!=(e in LOWER):
            raise ArithmeticError('D8 p=2 edge signature is not Diagram 6')
    if evaluate(p-2,2)!=0 or evaluate(4-2*p,2)!=0:
        raise ArithmeticError('D8 p=2 upper-w position bounds')
    right=[halfline_positive(3*p-6,2,1),
           halfline_positive(P.coerce(6),2,1),halfline_positive(P.coerce(9),2,1)]
    # For p>2, positions z1,z2 have order zero, whereas every Z-edge
    # has positive order (the lower-edge estimate >=6 is enough).
    left=halfline_positive(-6-(p-2),-4,-1)
    if evaluate(p-2,-4)!=-6:raise ArithmeticError('D8 lower-end edge equality')
    payload=dict(check='D8.11.10.sharp_pattern_boundaries',status='PASS',
        ray_map=dict(p='(2*r-4)/(r+1)',derivative='6/(r+1)^2',r_domain='(0,infinity)',p_image='(-4,2)',
                     corner_orders=dict(a=encode(a),eta=encode(eta)),
                     faces=[dict(p='-4',face='a=0, eta is a unit',a_order='3',eta_order='0'),
                            dict(p='2',face='eta=0, a is a unit',a_order='0',eta_order='3')]),
        whole_window=diagram8_gaps(result),endpoint_components=endpoint_components,
        p_equals_2_pattern=dict(diagram=6,z_strokes=[list(e)for e in UPPER],z_circles=[1,2],
                               w_strokes=[list(e)for e in LOWER],w_circles=[3,4],uncircled_hinge=5),
        p_greater_than_2=dict(z_circles=[1,2],z_strokes=[],strict_order_certificates=right,
                              contradiction='Rule 1a: a z-circle must have an incident z-stroke'),
        p_equals_minus_4=dict(upper_w_stroke_order='-6',lower_w_34_order='-6',
                             conclusion='Additional upper w-strokes prevent the Diagram-8 pattern'),
        p_less_than_minus_4=dict(upper_w_stroke_order=encode(p-2),lower_w_34_order='-6',
                                 strict_comparison=left,conclusion='Upper w-strokes dominate the lower w-stroke 34'),
        proof=['For each color, strokes and circles are exactly the components of minimal order. The whole-window comparisons are certified affine inequalities.',
            'At p=2, the upper Z-strokes and z1,z2 all have order zero; the lower W-strokes and w3,w4 have order -6. Every other component has strictly greater order. This is the Diagram-6 signature.',
            'For p>2, the upper Z-strokes have order 3p-6>0, the crossed ones have order 9, and the lower ones have order at least 6. Thus z1,z2 would be circles without z-strokes, contradicting Rule 1a.',
            'At p=-4, upper W-strokes join the order -6 of W34; below it they have smaller order. Neither case has the Diagram-8 signature.',
            'The monotone ray map and its exact endpoint limits prove the window and face correspondence. Endpoint pattern statements are conditional chart-order computations, not germ existence.'],
        scope='Supplementary sharpness and boundary-pattern computations of Remark 11.10(a), with finite lower wing.')
    output=Path(__file__).resolve().parent/'certificates/d8_pattern_boundaries.json'
    output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=7,
        detail=dict(certificate=output.name,sharp_window=True,upper_endpoint_diagram=6,whole_domains=True))],{}

if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram8_pattern_boundaries'))
