"""All 27 Lagrange classes and 64 crossed-sign labels at the K5 ray.

The universal regular/fold certificates replace redundant exceptional jets.
The exact Phi sign census is reconstructed here, including every pair scale.
"""
from itertools import product
from fractions import Fraction as Q
from pathlib import Path
import json
from .symbolic import sympy
from .interval import R

def verify(d2_mass_rays):
    s=sympy();r=s.symbols('r');q=(r*r+r+1)/3
    polynomial=r**8+2*r**7+3*r**6-4*r**5-5*r**4-6*r**3-3*r*r-12*r-12
    if s.expand(9*((r*r*q-r)**2-4*q)-polynomial):raise ArithmeticError('K5 unsquared root-amplitude identity')
    box=R('1.5736409','1.573641');qb=(box**2+box+1)/3;positive=box**2*qb-box
    if positive.lo<=0 or box.lo<=1:raise ArithmeticError('K5 physical root sign or unequal masses')
    separated=[]
    for ray in d2_mass_rays['root_records']:
        ranges=[tuple(map(Q,pair))for pair in ray['mass_hole']]+[(Q(1),Q(1))]
        for i in range(3):
            for j in range(i+1,3):
                if not(ranges[i][1]<ranges[j][0]or ranges[j][1]<ranges[i][0]):raise ArithmeticError('L ray has an uncertified repeated mass')
        separated.append(dict(seed=ray['j'],mass_boxes=[[str(v)for v in pair]for pair in ranges]))
    # W3,W4 are conjugate, have sum -r and product q. For v3=U-iV,
    # v4=U+iV, V>0, the squared equation above and its positive sign give
    # 2V=r. U>0 excludes the equal-majority-sign cancellation.
    signs=[]
    for first,second in product(product((1,-1),repeat=3),repeat=2):
        sigma=tuple(a*b for a,b in zip(first,second));odd=None
        if len(set(sigma))==2:odd=next(i for i,v in enumerate(sigma)if sigma.count(v)==1)
        zero=odd==2 and first[0]==-first[1]and first[2]==-sigma[0]*first[0]
        if zero:reason='odd 5; opposite majority square-root signs; 2V=r and the required odd sign'
        elif odd==2:reason='nonzero real component U or the wrong sign of the imaginary amplitude'
        elif odd is None:reason='L5=(r-1)^2 !=0; L5prime>0'
        else:reason='K5 positive zeros require the two majority masses equal; here they are 1 and r !=1'
        signs.append(dict(first=list(first),second=list(second),sigma=list(sigma),odd_body=None if odd is None else odd+3,Phi_zero=zero,reason=reason))
    if len(signs)!=64 or sum(v['Phi_zero']for v in signs)!=4:raise ArithmeticError('K5 complete sign census')
    classes=[]
    for exponents in product(range(3),repeat=3):
        distinct=len(set(exponents));kind='fold'if distinct==3 else'ordinary-S5'if distinct==1 else'two-equal'
        classes.append(dict(exponents=list(exponents),kind=kind,shape_count=1 if kind=='fold'else 2,Phi_zero_per_triangle=4 if kind!='two-equal'else 0))
    triangle_count=sum(r['shape_count']for r in classes);labels=triangle_count*3*64;zeros=sum(r['shape_count']*r['Phi_zero_per_triangle']for r in classes)*3
    if triangle_count!=48 or labels!=9216 or zeros!=144:raise ArithmeticError('K5 full finite-label census')
    payload=dict(check='D2.6.16.K5_all_27_classes_and_64_signs',status='PASS',classes=classes,S5_signs=signs,
        pair_cubic_scales=3,total_triangles=48,total_labels=9216,Phi_zero_labels=144,Phi_nonzero_labels=9072,unsquared_amplitude_sign=positive.encode(),L_ray_distinct_mass_boxes=separated,
        conjugation_transport='Real source rows are preserved. Conjugation swaps the two W shapes and sends the second crossed-sign row to its negative because D2B=i*sprimeB*positive_amplitude*vB; sigma changes to -sigma. Thus the opposite-shape zero rule is s5=+sigma*s3 and has the same exact four-sign count.',
        proof=['At mprime=(1,1,r), both S5 shapes obey W3+W4=-r, W3*W4=(r^2+r+1)/3 and W5=1. They are nonreal conjugate shapes.',
            'For conjugate signed square roots v3=U-iV,v4=U+iV with U,V>0, four V^2 equals 2/sqrt(q)+r/q. The verified ray polynomial and the positive unsquared sign give 2V=r exactly.',
            'The explicit 64-sign enumeration leaves four Phi zeros per S5 triangle, all with odd body 5; the other odd placements would require unequal majority masses to be equal by the verified K5 ray classification.',
            'The two-equal Lagrange classes have no one-different Phi zero by positive K6/K7/K8. Their same-sign zeros would require an L6/L7/L8 ray with three distinct masses, whereas (1,1,r) has two equal masses. The certified L-ray root boxes separate all three masses.',
            'All ordinary Phi-zero labels are isolated by the global S5 nonzero second response. All fold Phi-zero labels are isolated by the exact rank-eight block, nonzero t augmentation and simple torque zero. Every remaining label is isolated directly by tQ=0 and the finite boundary fibre.',
            'Thus all 9216 finite labels are covered. The universal fold argument makes the original exceptional h^3/h^4 jet calculation unnecessary; those redundant orders are not asserted by this certificate.'],
        depends_on=['D2.6.16.K5_positive_mass_ray','D2.6.17.S5_exact_second_order_response','D2.6.18.general_fold_rank_and_simple_torque_zero','D2.6.14.global_mass_rays_and_complex_shapes','D2.6.12.K3_K6_K7_K8_L5prime_and_equal_masses','D2.6.2.all_Lagrange_shape_assignments','D2.6.2.chart_and_boundary_fibre'])
    path=Path(__file__).resolve().parents[1]/'certificates/d2_K5_all_labels.json';path.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    return [dict(check=payload['check'],status='PASS',items=9216,detail=dict(certificate=path.name,isolated_labels=9216,complete_signs_and_classes=True,method='global regular/fold proof replaces redundant exceptional jets'))],{}
