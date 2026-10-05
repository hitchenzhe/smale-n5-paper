"""Table 5: exact denominator clearing and signs at Sturm-isolated roots."""
from fractions import Fraction as F
from smale_verify.qpoly import *
def positive_denominator_numerator(r):
    # Every denominator is a product of n, b, c, b+c, n^2+1 and positive constants.
    return primitive(r.n)

def sign_at_root(poly,root,seq,a,b):
    common=gcd(root,poly)
    if len(common)>1 and count_roots(sturm(common),a,b):return 0
    for _ in range(256):
        lo,hi=eval_interval(poly,(a,b))
        if lo>0:return 1
        if hi<0:return -1
        a,b=refine(root,seq,a,b)
    raise ArithmeticError("could not certify sign at isolated root")

