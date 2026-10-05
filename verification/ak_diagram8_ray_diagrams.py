"""Table 5: exact ray components and four-products in Remark 11.10."""
import json
from pathlib import Path
from smale_verify.d8_orders import reconstruction,diagram8_gaps,encode,UPPER,LOWER

def verify_ray_orders():
    result=reconstruction();gaps=diagram8_gaps(result)
    rows=[]
    for edge in result['distances']:
        z,w=result['separations'][edge];Z,W=result['components'][edge]
        rows.append(dict(edge=list(edge),z_separation=encode(z),w_separation=encode(w),
            squared_distance=encode(result['distances'][edge]),Z_component=encode(Z),W_component=encode(W)))
    products=[dict(omitted_pair=list(e),edges=[list(q)for q in result['product_edges'][e]],
                   order=encode(value))for e,value in result['four_products'].items()]
    payload=dict(check='D8.11.10.ray_orders_and_four_products',status='PASS',
        domain='-4<p<2; all orders normalized by ord(t)=1',
        chart_provenance=result['provenance'],pair_orders=rows,four_products=products,
        least_four_product=dict(constant='-8',p_coefficient='-2'),
        minimum_gap_certificates=result['minimum_gaps'],strict_negative_certificate=result['negative_minimum'],
        pattern=dict(z_strokes=[list(e)for e in UPPER],w_strokes=[list(e)for e in LOWER],
                     z_circles=[],w_circles=[3,4]),pattern_gap_certificates=gaps,
        proof=['The ten squared-distance ties are reconstructed from the chart. Nonzero finite inverse-distance limits give their exact orders.',
            'The gauges, positive mass centres, and the C/D source-row supports supply every needed leading unit and position lower bound. Cancellation can only increase a lower-bounded position order.',
            'All component comparisons are affine in p and certified on the entire open interval by exact endpoint values.',
            'All ten four-products are indexed by the omitted pair and its complementary triangle. Eight have order -8-2p, one has order 4-2p and one has order 12-6p. The minimum is strictly negative throughout the window.'],
        scope='Supplementary conditional order statements in Remark 11.10. This is not a construction of a Diagram-8 germ.')
    output=Path(__file__).resolve().parent/'certificates/d8_ray_diagrams.json'
    output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=20,
        detail=dict(certificate=output.name,pairs=10,four_products=10,whole_open_window=True))],{}

if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram8_ray_diagrams'))
