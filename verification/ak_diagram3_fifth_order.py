"""Table 5: exact fifth-order jets, nine-root norm and Lemma 7.10."""
from itertools import product
from fractions import Fraction as F
from smale_verify.exact import P, V, substitute, equal, zero, quotient_norm, det, diff_poly
from smale_verify.qpoly import divmodp, power, sturm, variations, primitive


def verify_nine_root_norm():
    a, b, y = map(V, ('d35_a', 'd35_b', 'd35_y'))
    f4 = a**4+a**3*b+a**2*b**2+a*b**3+b**4
    f2 = a**2+a*b+b**2
    A, B, C = a**2*f4, -b**4/a**2*f4, -b**2*(a**2-b**2)*f2
    r, s = (a**2+b**2)/(2*a**2), (a**2+b**2)/(2*b**2)
    # Norm first over x^3=r, then over y^3=s. The first equality follows
    # by a separate multiplication-matrix determinant, rather than assuming
    # the product formula used in the second norm.
    x = V('d35_x')
    inner = (B*y+C)**3+A**3*r
    equal(quotient_norm(A*x+B*y+C, 'd35_x', [-r,0,0,1]), inner,
          'D3 first cubic norm')
    norm = quotient_norm(inner, 'd35_y', [-s,0,0,1])
    cleared = norm*a**18
    if any(sum(e for _,e in mon) != 72 for mon in cleared.terms):
        raise AssertionError('norm homogeneity')
    at_one = substitute(cleared, {'d35_a':P.coerce(1)})
    coeffs = [at_one.coefficient('d35_b', n).terms.get((), F(0)) for n in range(73)]
    q, rem = divmodp(coeffs, power([-1,0,1], 3))
    if rem or len(q) != 67:
        raise AssertionError('norm degree-66 factor')
    coefficients = [-8*c for c in q]
    if any(c <= 0 or c.denominator != 1 for c in coefficients):
        raise AssertionError('P66 coefficients not positive integers')
    homogeneous = sum((c*a**(66-k)*b**k for k,c in enumerate(coefficients)), P())
    equal(norm, (a-b)**3*(a+b)**3*homogeneous/(8*a**18),
          'D3 full homogeneous norm factorization')
    minus = [c*(-1)**k for k,c in enumerate(coefficients)]
    sequence = sturm(minus)
    at_zero = variations(sequence, F(0))
    signs = [1 if row[-1] > 0 else -1 for row in sequence]
    at_infinity = sum(x != y for x,y in zip(signs, signs[1:]))
    if at_zero != at_infinity or minus[0] == 0:
        raise AssertionError('P66(1,-rho) has a positive real root')
    return [dict(check='D3.7.10.nine_root_norm_and_Sturm',status='PASS',items=3,
        detail=dict(degree=66, positive_coefficients=67,
            P66_coefficients_ascending=[str(c) for c in coefficients],
            norm='(a-b)^3(a+b)^3 P66(a,b)/(8a^18)',
            minus_positive_roots=0, variations_zero=at_zero,
            variations_positive_infinity=at_infinity,
            sturm_chain=[[str(c) for c in row] for row in sequence],
            scope='norm and positivity; the companion jet function independently proves the G5-to-L identity'))], {}


def verify_fifth_order_jet():
    from ak_diagram3_blowup import compute_second_jet, selected_system
    data,variables,X2=compute_second_jet()
    rows=selected_system(data);m=data['m'];D=data['deltas'];normalize=data['normalize']
    J=[[normalize(substitute(diff_poly(row.coefficient('t3chart',0),v),data['point'])) for v in variables] for row in rows]
    denominator=-5184*m[1]**2*m[2]**3*m[4]**3*(D[(1,2)]*D[(1,4)]*D[(2,3)]*D[(3,4)])**4
    denominator=denominator*D[(1,3)]*D[(2,4)]*D[(1,5)]*D[(2,5)]*D[(3,5)]*D[(4,5)]/(data['s12']*data['s14'])
    rhs=[normalize(-substitute(row.coefficient('t3chart',3),data['point'])) for row in rows]
    X3=[]
    for col in range(19):
        replaced=[[rhs[i] if j==col else J[i][j] for j in range(19)] for i in range(19)]
        X3.append(normalize(det(replaced)/denominator))
    for i in range(19):
        zero(data['reduce_label'](sum((J[i][j]*X3[j] for j in range(19)),P())-rhs[i]),'D3 third jet row '+str(i))
    rz,rw=data['rz'],data['rw']
    left=[rz[1]-rz[2]-m[3]/m[2]*rz[3]+m[4]/m[1]*rz[4],
          rw[1]+m[2]/m[1]*rw[2]-m[2]/m[1]*rw[3]-rw[4]]
    G5=[]
    for ell in left:
        ell=substitute(ell,{'d3c_B':data['c']['A']+1})
        zero(ell.coefficient('t3chart',0),'D3 left row identically zero')
        value=substitute(ell.coefficient('t3chart',5),data['point'])
        for order,jet in ((3,X2),(2,X3)):
            value+=sum((substitute(diff_poly(ell.coefficient('t3chart',order),v),data['point'])*q for v,q in zip(variables,jet)),P())
        G5.append(normalize(value))
    a,b,c,S,H,ii=map(V,('d35f_a','d35f_b','d35f_c','d35f_S','d35f_H','d35f_i'))
    root_relations={c.names().pop():(2,2*S),ii.names().pop():(2,P.coerce(-1)),
        'd3c_D12':(3,(2*a**2)**-1),'d3c_D34':(3,(2*b**2)**-1),'d3c_D14':(3,S**-1)}
    def family_reduction(poly,branch):
        result=P()
        for mon,coefficient in poly.terms.items():
            term=P.coerce(coefficient)
            for name,exponent in mon:
                if name in root_relations:
                    power,value=root_relations[name];q,r=divmod(exponent,power)
                    term*=V(name)**r*value**q
                else:term*=V(name)**exponent
            result+=term
        for name in ('d35f_a','d35f_b','d35f_c','d35f_S','d35f_H'):
            result=result.clear(name)
        return substitute(result,{'d35f_S':a**2+b**2,'d35f_H':a**3+branch*b**3})
    for branch in (1,-1):
        bb=branch*b
        f4=a**4+a**3*bb+a**2*bb**2+a*bb**3+bb**4
        f2=a**2+a*bb+bb**2
        L=a**2*f4*D[(1,2)]-b**4/a**2*f4*D[(3,4)]-b**2*(a**2-b**2)*f2*D[(1,4)]
        for nu in ((1,) if branch==1 else (1,-1)):
            for sigma,eta in product((1,-1),repeat=2):
                mapping={'m31':a**2,'m32':a**2,'m34':b**2,'m35':a**3*b**3/(c*nu*H),
                    's312':2*a**2,'s314':S,'d3c_D13':P.coerce(sigma),'d3c_D24':eta*ii,
                    'd3c_D15':-sigma*nu*c/b,'d3c_D35':-branch*sigma*nu*c/a,
                    'd3c_D25':-eta*nu*ii*c/b,'d3c_D45':-branch*eta*nu*ii*c/a,
                    'd3c_D23':D[(1,4)]}
                expected=3*(sigma-eta*ii)*(a-branch*b)/H*L
                zero(family_reduction(substitute(G5[0],mapping)-expected,branch),'D3 general fifth-order factor')
    # Check the other families by exact covariance of all 21 chart rows.
    t=V('t3chart');coordinates=data['c'];perm={1:1,2:4,3:3,4:2,5:5}
    exchange={'m32':data['md'],'m34':data['mb'],
        'd3c_A':coordinates['C'],'d3c_B':coordinates['D']+t**2*coordinates['v'],
        'd3c_x':coordinates['u'],'d3c_y':-coordinates['v'],'d3c_zhat5':coordinates['what5'],
        'd3c_C':coordinates['A'],'d3c_D':coordinates['B']+t**2*coordinates['y'],
        'd3c_u':coordinates['x'],'d3c_v':-coordinates['y'],'d3c_what5':coordinates['zhat5']}
    for pair in D:
        exchange['d3c_D%d%d'%pair]=D[tuple(sorted(perm[j] for j in pair))]
    for j in m:
        equal(substitute(rz[j],exchange),rw[perm[j]],'D3 z/w exchange z row')
        equal(substitute(rw[j],exchange),rz[perm[j]],'D3 z/w exchange w row')
    for pair,row in data['ties'].items():
        equal(substitute(row,exchange),data['ties'][tuple(sorted(perm[j] for j in pair))],'D3 z/w exchange tie')
    equal(substitute(data['gauge'],exchange),-data['gauge'],'D3 z/w exchange gauge')
    for name in exchange:
        equal(substitute(substitute(V(name),exchange),exchange),V(name),
              'D3 chart exchange is an analytic involution')
    # At equal masses both components vanish over the complete root ring,
    # on each of the two overlapping sets making up the 180-label census.
    for side in ('d3c_D23','d3c_D34'):
        mapping={'m31':P.coerce(F(1,2)),'m32':P.coerce(F(1,2)),'m34':P.coerce(F(1,2)),
            'm35':P.coerce(F(1,8)),'s312':P.coerce(1),'s314':P.coerce(1),
            'd3c_D15':-2*D[(1,3)],'d3c_D35':-2*D[(1,3)],
            'd3c_D25':-2*D[(2,4)],'d3c_D45':-2*D[(2,4)],
            side:D[(1,4)] if side=='d3c_D23' else D[(1,2)]}
        for value in G5:
            value=substitute(value,mapping)
            for name in ('d3c_D12','d3c_D14','d3c_D23','d3c_D34'):
                value=value.reduce(name,[-1,0,0,1])
            value=value.reduce('d3c_D13',[-1,0,1]).reduce('d3c_D24',[1,0,1])
            zero(value,'D3 equal-mass fifth order vanishes')
    witnesses={name:[dict(coefficient=str(coef),powers=dict(mon)) for mon,coef in sorted(value.terms.items())] for name,value in zip(variables,X3)}
    return [dict(check='D3.7.10.fifth_order_jet_and_symmetry',status='PASS',items=56,
        detail=dict(third_jet=witnesses,method='exact Cramer solves and independent replay of all 19 linear rows',
            G5_z='3(sigma-i eta)(a-branch*b)/(a^3+branch*b^3) * D14 * L(D12/D14,D34/D14); branch=+1 or -1, with b replaced by branch*b in F4,F2',
            factor_nonzero='a,b>0; on unequal masses a-b or a+b is nonzero; sigma-i eta has modulus sqrt(2); D14 nonzero',
            symmetry_rows=21,equal_labels_G5_zero=180,
            scope='all four positive mass families and all sign/cubic root classes'))],{}


if __name__ == '__main__':
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram3_fifth_order'))

