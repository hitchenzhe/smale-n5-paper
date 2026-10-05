"""Table 5: the equal upper-pair specialization of Lemma 11.6."""
from fractions import Fraction as F
from smale_verify.exact import P,V,zero,equal,substitute,diff_poly
from ak_diagram8_relation40_derivation import chart_rows

def verify_equal_order_six_pairing():
    m={1:P.coerce(1),2:P.coerce(1),3:V('m3'),4:V('m4'),5:P.coerce(F(1,4))}
    data=chart_rows(m);x=data['x']
    point=dict(zip(('d8_A1','d8_A2','d8_C','d8_R12','d8_R25','d8_P','d8_Q'),map(F,('-1/2','1/2','0','-1/2','1','-8','8'))))
    U=('d8_A1','d8_A2','d8_C','d8_R12','d8_R25','d8_P','d8_Q')
    indices=(0,2,4,10,13,16,20)
    center=x['A1']+x['A2']+(m[3]+m[4]+m[5])*x['C']
    rows=[data['f'][j] for j in indices]+[center]
    ell=list(map(F,('1','-1','0','3/8','-3/16','-3/16','0','0')))
    for v in U:
        zero(sum((q*substitute(diff_poly(row,v),point) for row,q in zip(rows,ell)),P()),'D8 equal order-six left kernel '+v)
    S6=[substitute(data['E'][j].coefficient('a8',6).coefficient('eta8',0),point) for j in indices]+[P()]
    equal(sum((q*value for q,value in zip(ell,S6)),P()),-1,'D8 equal order-six pairing')
    return [dict(check='D8.11.6.equal_pairing',status='PASS',items=8,
                 detail=dict(left_kernel=list(map(str,ell)),pairing='-1',
                             masses='(1,1,m3,m4,1/4); all positive scales follow by mass scaling'))],{}

if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram8_equal_pair_germs'))
