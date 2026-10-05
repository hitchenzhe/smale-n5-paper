"""Table 5 first_limit_relevance: exact far-source expansions.

Reconstructs first and second coefficients, including signed-moment zeros.
Near-source coefficient identities and nonvanishing are independently
reconstructed by verify_near_source_audit.
"""
from fractions import Fraction as Q
from itertools import product
from smale_verify.exact import P,V,equal,zero

def verify_near_source_audit():
 from smale_verify.d5_source_audit import verify
 return verify()


def verify_far_source_expansions():
 h,z,w=map(V,('d5far_h','d5far_z','d5far_w'))
 masses=[V('d5far_m'+str(i)) for i in range(3)]
 offsets=[V('d5far_s'+str(i)) for i in range(3)]
 records=[]
 for signs in product((1,-1),repeat=3):
  moments=[sum((m*e*s**n for m,e,s in zip(masses,signs,offsets)),P()) for n in range(3)]
  mu,C1,C2=moments;rows={name:P() for name in ('F','G','V1','V2','B')}
  for m,e,s in zip(masses,signs,offsets):
   d1=1+h*s*(z**-1+w**-1)/2+h**2*s**2*(Q(3,8)*(z**-2+w**-2)+Q(1,4)*z**-1*w**-1)
   d3=1+Q(3,2)*h*s*(z**-1+w**-1)+h**2*s**2*(Q(15,8)*(z**-2+w**-2)+Q(9,4)*z**-1*w**-1)
   for n in range(3):
    zero((d1**2*(1-h*s/z)*(1-h*s/w)-1).coefficient('d5far_h',n),'D5 distance-ratio series')
    zero((d1**3-d3).coefficient('d5far_h',n),'D5 inverse-cube series')
   rows['F']+=m*e*d3*(z-h*s);rows['G']+=m*e*d3*(w-h*s)
   rows['V1']+=m*e*d3*h*s*(h*s-w);rows['V2']+=m*e*d3*h*s*(h*s-z)
   bracket=1+h*s*(z**-1+w**-1)/2+h**2*s**2*(z**-2+w**-2)/2
   rows['B']+=m*e*d1*bracket
  expected=dict(F=[mu*z,C1*(w+3*z)/(2*w)],G=[mu*w,C1*(z+3*w)/(2*z)],
    V1=[P(),-C1*w,-C2*(1+3*w/z)/2],V2=[P(),-C1*z,-C2*(1+3*z/w)/2],
    B=[mu,C1*(z**-1+w**-1),C2*(Q(9,8)*(z**-2+w**-2)+Q(3,4)*z**-1*w**-1)])
  for name,coefficients in expected.items():
   for n,value in enumerate(coefficients):equal(rows[name].coefficient('d5far_h',n),value,'D5 far '+name+' coefficient '+str(n))
  minority=None
  if len(set(signs))==2:
   minority=next(i for i in range(3) if signs.count(signs[i])==1)
   majority=-signs[minority]
   equal(C1,majority*sum((m*s for m,s in zip(masses,offsets)),P())-2*majority*masses[minority]*offsets[minority],'D5 minority signed moment')
  else:
   equal(C1,signs[0]*sum((m*s for m,s in zip(masses,offsets)),P()),'D5 equal-sign signed moment')
  records.append(dict(signs=list(signs),minority_index=minority,
    first_signed_moment='zero by mass centering' if minority is None else 'minus twice the minority signed mass times its centered offset',
    slice_exception='first signed moment zero when the minority body is central' if minority is not None else None))
 return [dict(check='D5.far_source_expansions_and_signed_moment_cases',status='PASS',items=8,
  detail=dict(sign_cases=records,coefficients='inverse-distance and inverse-cube ratios through order two; F/G, V1/V2 and B checked by exact expansion',
    scope='formal coefficients only; complete valuation branch audit and LP certificates remain'))],{}


if __name__=='__main__':
 from smale_verify.runner import program_cli
 raise SystemExit(program_cli('ak_diagram5_first_limit_relevance'))

