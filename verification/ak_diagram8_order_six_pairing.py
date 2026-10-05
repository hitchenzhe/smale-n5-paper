"""Table 5: support, rank and order-six left pairing in Lemma 11.6."""
from smale_verify.exact import P,zero,equal,substitute,diff_poly,det
from ak_diagram8_relation40_derivation import chart_rows
from ak_diagram8_kappa_row import general_upper_data,clear_substitute

def verify_row_support(m):
    data=chart_rows(m);E,f,x=data['E'],data['f'],data['x']
    supports={0:{(6,0),(2,2),(9,3),(11,5)},2:{(6,0),(2,2),(9,3),(11,5)},
              4:{(2,2),(3,3),(5,5)},13:{(2,2)},16:{(2,2)},10:set(),20:set()}
    for j,allowed in supports.items():
        found={tuple(dict(mon).get(v,0) for v in ('a8','eta8')) for mon in (E[j]-f[j]).terms}
        if found!=allowed:raise AssertionError(('D8 order-six row support',j,found,allowed))
    return [dict(check='D8.11.6.row_support',status='PASS',items=7,
                 detail=dict(supports={str(k):sorted(v) for k,v in supports.items()}))],{}

def verify_general_left_pairing():
    data,m,point,alpha,s,g,B=general_upper_data()
    x=data['x'];E=data['E']
    U=('d8_A1','d8_A2','d8_C','d8_R12','d8_R25','d8_P','d8_Q')
    indices=(0,2,4,10,13,16,20)
    center=m[1]*x['A1']+m[2]*x['A2']+(m[3]+m[4]+m[5])*x['C']
    faces=[data['f'][j] for j in indices]+[center]
    mapping={'s8':alpha**2+1,'g8':(alpha**2+1)**3+alpha**3}
    for j,row in enumerate(faces):
        zero(clear_substitute(substitute(row,point),mapping),'D8 general S face '+str(j))
    J=[[substitute(diff_poly(row,v),point) for v in U] for row in faces]
    h=[(-1)**i*det([row for j,row in enumerate(J) if j!=i]) for i in range(8)]
    for col in range(7):
        zero(sum((h[j]*J[j][col] for j in range(8)),P()),'D8 cofactor left kernel column '+str(col))
    # This cofactor is nonzero for alpha>0, alpha!=1: it is the previously
    # certified det(B) times the positive total mass.  At alpha=1 the separate
    # equal-pair program exhibits a constant rank-seven minor instead.
    zero(clear_substitute(h[2]-(s+m[3]+m[4]+m[5])*det(B),mapping),'D8 nonzero normalizing cofactor')
    S6=[substitute(E[j].coefficient('a8',6).coefficient('eta8',0),point) for j in indices]+[P()]
    for j,value in enumerate(S6):
        target=point['d8_A1'] if j==0 else point['d8_A2'] if j==1 else P()
        equal(value,target,'D8 mu coefficient '+str(j))
    Delta=(m[3]+m[4])*s**2+g
    pairing=sum((h[j]*S6[j] for j in range(8)),P())
    zero(clear_substitute((alpha**2-1)*pairing*alpha*Delta+2*s**4*h[2],mapping),'D8 normalized order-six pairing')
    return [dict(check='D8.11.6.general_rank_and_pairing',status='PASS',items=25,
                 detail=dict(normalization='ell_E4=alpha^2-1',
                             pairing='-2 (alpha^2+1)^4/(alpha Delta)',
                             Delta='(m3+m4)(alpha^2+1)^2+alpha^6+3alpha^4+alpha^3+3alpha^2+1',
                             rank_certificate='det(B) and C row for alpha!=1; constant augmented minor at alpha=1',
                             left_kernel='eight signed maximal minors, checked in all seven columns'))],{}

if __name__=='__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram8_order_six_pairing'))
