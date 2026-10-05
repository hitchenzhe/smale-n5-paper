"""Exact integer Kronecker encoding of bounded multivariate resultants.

The coefficient bound is n! ||f||_1^deg(g) ||g||_1^deg(f).
With base > 2*bound, balanced digits recover every coefficient uniquely.
Homogeneous mass degree lets m5=1 recover its exponent from m3,m4.
No modular probability or floating arithmetic is used.
"""
from math import factorial
from fractions import Fraction


def integer_det(matrix):
    a=[list(row) for row in matrix];n=len(a);previous=1;sign=1
    for k in range(n-1):
        if not a[k][k]:
            row=next((j for j in range(k+1,n) if a[j][k]),None)
            if row is None:return 0
            a[k],a[row]=a[row],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value=pivot*a[i][j]-a[i][k]*a[k][j]
                quotient,remainder=divmod(value,previous)
                if remainder:raise ArithmeticError('Nonexact Bareiss division')
                a[i][j]=quotient
            a[i][k]=0
        previous=pivot
    return sign*a[-1][-1]


def parameter_resultant(symbolic,first,second,variable,parameter):
    """Two integer bivariate polynomials; recover the exact parameter resultant."""
    s=symbolic;f=s.Poly(first,variable,parameter);g=s.Poly(second,variable,parameter)
    n,m=f.degree(variable),g.degree(variable)
    degree_bound=f.degree(parameter)*m+g.degree(parameter)*n
    bound=factorial(n+m)*sum(abs(int(v)) for _,v in f.terms())**m*sum(abs(int(v)) for _,v in g.terms())**n
    bits=(2*bound).bit_length()+1;base=1<<bits
    def packed(poly):
        return [int(v.subs(parameter,base)) for v in s.Poly(poly.as_expr(),variable).all_coeffs()]
    a,b=packed(f),packed(g);matrix=[]
    for coefficients,shifts in ((a,m),(b,n)):
        for shift in range(shifts):matrix.append([0]*shift+coefficients+[0]*(shifts-1-shift))
    value=integer_det(matrix);sign=-1 if value<0 else 1;value=abs(value);coefficients=[]
    while value:
        digit=value%base;value//=base
        if digit>=base//2:digit-=base;value+=1
        if abs(digit)>bound:raise ArithmeticError('Parameter resultant coefficient bound failed')
        coefficients.append(sign*digit)
    if len(coefficients)>degree_bound+1:raise ArithmeticError('Parameter resultant degree bound failed')
    return coefficients,dict(coefficient_bound=str(bound),degree_bound=degree_bound,radix_bits=bits)


def resultant(f,g,shape_degree,divide_m3=True):
    d=max(t['powers'].get('x',0) for t in f)
    e=max(t['powers'].get('x',0) for t in g)
    raw_mass_degree=e*shape_degree+6*d
    mass_degree=raw_mass_degree-(2*shape_degree if divide_m3 else 0)
    y_degree=d*max(t['powers'].get('y',0) for t in g)
    bound=factorial(d+e)*sum(abs(int(t['coefficient'])) for t in f)**e*sum(abs(int(t['coefficient'])) for t in g)**d
    byte_width=max(1,((2*bound).bit_length()+7)//8)
    bits=8*byte_width;radix=raw_mass_degree+1
    def evaluate(poly,degree):
        coefficients=[0]*(degree+1)
        for term in poly:
            p=term['powers'];value=Fraction(term['coefficient'])
            if value.denominator!=1:raise ArithmeticError('Packed resultant needs integral inputs')
            exponent=p.get('m3',0)+radix*p.get('m4',0)+radix**2*p.get('y',0)
            coefficients[p.get('x',0)]+=value.numerator << (bits*exponent)
        return list(reversed(coefficients))
    a,b=evaluate(f,d),evaluate(g,e)
    matrix=[]
    for coeffs,shift_count in ((a,e),(b,d)):
        for shift in range(shift_count):matrix.append([0]*shift+coeffs+[0]*(shift_count-1-shift))
    value=integer_det(matrix)
    sign=-1 if value<0 else 1;raw=abs(value).to_bytes(max(1,(abs(value).bit_length()+7)//8),'little')
    terms=[];carry=0;limit=radix**2*(y_degree+1)
    for index in range(max((len(raw)+byte_width-1)//byte_width,1)+1):
        digit=int.from_bytes(raw[index*byte_width:(index+1)*byte_width],'little')+carry
        if digit >= (1 << (bits-1)):digit-=1<<bits;carry=1
        else:carry=0
        coefficient=sign*digit
        if not coefficient:continue
        if abs(coefficient)>bound:raise ArithmeticError('Decoded coefficient exceeds determinant bound')
        yp,remainder=divmod(index,radix**2);m4p,m3p=divmod(remainder,radix)
        m5p=raw_mass_degree-m3p-m4p
        if m5p<0 or yp>y_degree or index>=limit:raise ArithmeticError('Decoded monomial exceeds degree bounds')
        if divide_m3:
            m3p-=2*shape_degree
            if m3p<0:raise ArithmeticError('Resultant is not symbolically divisible by claimed m3 power')
        powers={k:v for k,v in dict(y=yp,m3=m3p,m4=m4p,m5=m5p).items() if v}
        terms.append(dict(coefficient=str(coefficient),powers=powers))
    if carry:raise ArithmeticError('Unconsumed balanced digit carry')
    return terms,dict(coefficient_bound=str(bound),radix_bits=bits,mass_degree=mass_degree,raw_mass_degree=raw_mass_degree,y_degree=y_degree,
                     method='exact integer Sylvester determinant; injective balanced Kronecker encoding')
