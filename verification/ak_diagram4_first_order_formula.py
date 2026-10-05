"""Table 5: exact complex first-order vectors X,Y in Diagram 4."""
from fractions import Fraction as F
from smale_verify.qpoly import *
def cadd(a,b):return (a[0]+b[0],a[1]+b[1])
def cneg(a):return (-a[0],-a[1])
def csub(a,b):return cadd(a,cneg(b))
def cmul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def cpow(a,n):
    y=(Rat(1),Rat(0))
    for _ in range(n):y=cmul(y,a)
    return y

def vectors(signs,m,g4,g5,power):
    s13,s14,s15,s23,s24,s25=signs
    im_sign=1 if power==1 else -1
    X=(s13*m[0]+s15*m[2]*g5**power,im_sign*s14*m[1]*g4**power)
    Y=(s23*m[0]+s25*m[2]*g5**power,im_sign*s24*m[1]*g4**power)
    return X,Y

