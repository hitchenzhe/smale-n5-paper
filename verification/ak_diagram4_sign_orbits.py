"""Table 5: reconstruct all 64 signs from the 16 normalized cases."""
from itertools import product
def positive_variants(signs,R_sign):
    _,s14,s15,_,s24,s25=signs
    return [(f1,f1*s14,f1*s15,f2,f2*s24,f2*s25)
            for f1,f2 in product((1,-1),repeat=2) if f1*f2*R_sign>0]
