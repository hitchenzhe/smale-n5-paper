"""Reconstruction of ak_diagram3_equal_shape from Table 5, paper page 56.

Only the checks reported by this module are implemented; see program_index.json.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from smale_verify.exact import (P, V, zero, equal, edge, equations, substitute,
    diff_poly, divide_poly, groebner, det, quotient_norm)

def verify_twisted_family():
    checks=[]
    def record(name,count,detail=None):
        checks.append(dict(check=name,status="PASS",items=count,detail=detail))
    # Exact check of the explicit twisted Diagram-3 family in Lemma 7.12.
    tau,pp,qq,yy=map(V,('tau','P3','Q3','Y3'))
    alpha,beta,eps,sgn=map(V,('alpha','beta','eps','eta'))
    d3_relations=[(1-pp)*(1+qq)-alpha*tau**2,
                  (1+pp)*(1-qq)-beta*tau**2,
                  4*pp*qq*yy**2+1,alpha**3-1,beta**3-1,eps**2-1,sgn**2-1]
    d3_order=('alpha','beta','eps','eta','tau','P3','Q3','Y3')
    d3_basis=groebner(d3_relations,d3_order)
    d3_z={1:P.coerce(1),2:pp,3:P.coerce(-1),4:-pp,5:P()}
    d3_w={1:P.coerce(1),2:-qq,3:P.coerce(-1),4:qq,5:P()}
    d3_D={(1,2):alpha,(3,4):alpha,(1,4):beta,(2,3):beta,
           (1,3):eps*tau/2,(1,5):-eps*tau,(3,5):-eps*tau,
           (2,4):tau*yy,(2,5):-2*tau*yy,(4,5):-2*tau*yy}
    d3_m={1:F(1,2),2:F(1,2),3:F(1,2),4:F(1,2),5:F(1,8)}
    d3_rows=[]
    for j in range(1,6):
        d3_rows.append(d3_z[j]-sum(d3_m[l]*d3_D[edge(j,l)]**3*(d3_z[j]-d3_z[l]) for l in range(1,6) if l!=j))
        d3_rows.append(d3_w[j]-sum(d3_m[l]*d3_D[edge(j,l)]**3*(d3_w[j]-d3_w[l]) for l in range(1,6) if l!=j))
    for e in combinations(range(1,6),2):
        i,j=e
        d3_rows.append(d3_D[e]**2*(d3_z[j]-d3_z[i])*(d3_w[j]-d3_w[i])-tau**2)
    zero(d3_z[3]-d3_z[1]-d3_w[3]+d3_w[1],'D3 twisted gauge')
    for idx,row in enumerate(d3_rows):
        rem=divide_poly(row,d3_basis,d3_order)
        zero(rem,'D3.7.12 row '+str(idx+1))
    record('D3.7.12.twisted_family',len(d3_rows),{'relations':len(d3_relations),'groebner_basis_size':len(d3_basis),'method':'exact ideal reduction over Q'})
    return checks, {}

def verify_modular_order_six():
    """Replay Lemma 7.11 at every equal-mass label, including all jet data."""
    from ak_diagram3_blowup import boundary_data, selected_system
    from smale_verify.exact import diff_poly
    from smale_verify.modular import (rational_mod, inverse_matrix,
        compile_polynomials, evaluate_series)
    data = boundary_data()
    selected = selected_system(data)
    raw = list(data['rz'].values()) + list(data['rw'].values())
    raw += [data['ties'][pair] for pair in sorted(data['ties'])] + [data['gauge']]
    raw = [substitute(row, {'d3c_B': data['c']['A']+1}) for row in raw]
    variables = ['d3c_'+n for n in ('A','x','y','zhat5','C','D','u','v','what5')]
    variables += ['d3c_D%d%d'%pair for pair in sorted(data['deltas'])]
    certificates = []
    for p in (13, 37):
        half = pow(2, -1, p)
        constants = dict(m31=half, m32=half, m34=half,
                         m35=pow(8, -1, p), s312=1, s314=1)
        cube_roots = [x for x in range(p) if pow(x, 3, p) == 1]
        iota = next(x for x in range(p) if x*x % p == p-1)
        if len(cube_roots) != 3:
            raise AssertionError('missing cubic roots')
        compiled_selected = compile_polynomials(selected, variables, constants, 't3chart', p)
        compiled_raw = compile_polynomials(raw, variables, constants, 't3chart', p)
        jacobian_rows = [diff_poly(row.coefficient('t3chart', 0), name)
                         for row in selected for name in variables]
        compiled_J = compile_polynomials(jacobian_rows, variables, constants, 't3chart', p)
        labels = []
        for cubes in product(cube_roots, repeat=4):
            d12, d14, d23, d34 = cubes
            if d12 != d34 and d14 != d23:
                continue
            for sign13, sign24 in product((1, -1), repeat=2):
                d13 = sign13 % p
                d24 = sign24*iota % p
                deltas = {(1,2):d12, (1,3):d13, (1,4):d14,
                          (1,5):-2*d13 % p, (2,3):d23, (2,4):d24,
                          (2,5):-2*d24 % p, (3,4):d34,
                          (3,5):-2*d13 % p, (4,5):-2*d24 % p}
                values = dict(constants, **{'d3c_D%d%d'%pair: value
                                            for pair, value in deltas.items()})
                point = []
                for name in variables:
                    if name not in data['point']:
                        point.append(values[name])
                    else:
                        value = 0
                        for mon, coefficient in data['point'][name].terms.items():
                            term = rational_mod(coefficient, p)
                            for key, exponent in mon:
                                term = term*pow(values[key], exponent, p) % p
                            value = (value+term) % p
                        point.append(value)
                jets = [[x]+[0]*6 for x in point]
                if any(row[0] for row in evaluate_series(compiled_raw, jets, 0, p)):
                    raise AssertionError('modular boundary label')
                entries = [row[0] for row in evaluate_series(compiled_J, jets, 0, p)]
                J = [entries[19*i:19*(i+1)] for i in range(19)]
                inverse, minor, pivots = inverse_matrix(J, p)
                # Independent specialization of the previously proven determinant.
                expected = (-5184*pow(half, 8, p)
                            *pow(d12*d14*d23*d34, 4, p)
                            *d13*d24*deltas[(1,5)]*deltas[(2,5)]
                            *deltas[(3,5)]*deltas[(4,5)]) % p
                if minor != expected:
                    raise AssertionError('selected unit minor formula')
                for order in range(1, 7):
                    rhs = [-row[order] % p for row in
                           evaluate_series(compiled_selected, jets, order, p)]
                    solution = [sum(x*y for x, y in zip(row, rhs)) % p for row in inverse]
                    for j in range(19):
                        jets[j][order] = solution[j]
                    selected_residual = evaluate_series(compiled_selected, jets, order, p)
                    if any(any(row) for row in selected_residual):
                        raise AssertionError('selected jet row replay')
                residual = evaluate_series(compiled_raw, jets, 6, p)
                if any(any(row[:6]) for row in residual):
                    raise AssertionError('raw residual before order six')
                sixth = [row[6] for row in residual]
                survives = not any(sixth)
                if survives != (d12 == d34 and d14 == d23):
                    raise AssertionError('order-six exclusion differs from Lemma 7.11')
                labels.append(dict(cube_indices=[cube_roots.index(x) for x in cubes],
                    signs=[sign13, sign24], unit_minor=minor, elimination_pivots=pivots,
                    jets=jets, residual_order_six=sixth, survives_order_six=survives))
        survivors = sum(label['survives_order_six'] for label in labels)
        if len(labels) != 180 or survivors != 36:
            raise AssertionError('equal-mass label census')
        certificates.append(dict(prime=p, iota=iota, cube_roots=cube_roots,
            labels=labels, excluded=180-survivors, survivors=survivors))
    for left, right in zip(certificates[0]['labels'], certificates[1]['labels']):
        if (left['cube_indices'], left['signs'], left['survives_order_six']) != (
                right['cube_indices'], right['signs'], right['survives_order_six']):
            raise AssertionError('prime-independent excluded label set')
    return [dict(check='D3.7.11.modular_order_six', status='PASS', items=360,
        detail=dict(variables=variables, order=6, primes=[13,37],
            masses=['1/2','1/2','1/2','1/2','1/8'],
            per_prime_labels=180, per_prime_excluded=144, per_prime_survivors=36,
            method='exact prime-field series; all selected rows replayed; certified unit minors; all 21 raw residuals recorded',
            scope='nonzero modular residual proves a characteristic-zero exclusion; zero residual alone is not an existence proof',
            certificates=certificates))], {}


if __name__ == "__main__":
    from smale_verify.runner import program_cli
    raise SystemExit(program_cli('ak_diagram3_equal_shape'))
