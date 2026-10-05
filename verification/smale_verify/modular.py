"""Exact truncated power series and linear algebra over a prime field."""

def rational_mod(value, p):
    return value.numerator * pow(value.denominator, -1, p) % p


def inverse_matrix(matrix, p):
    n = len(matrix)
    a = [[x % p for x in row] + [int(i == j) for j in range(n)]
         for i, row in enumerate(matrix)]
    determinant = 1
    pivots = []
    for j in range(n):
        k = next((k for k in range(j, n) if a[k][j]), None)
        if k is None:
            raise ArithmeticError('singular matrix modulo ' + str(p))
        if k != j:
            a[j], a[k] = a[k], a[j]
            determinant = -determinant
        pivot = a[j][j]
        pivots.append(pivot)
        determinant = determinant * pivot % p
        a[j] = [x * pow(pivot, -1, p) % p for x in a[j]]
        for i in range(n):
            if i != j:
                factor = a[i][j]
                if factor:
                    a[i] = [(x - factor*y) % p for x, y in zip(a[i], a[j])]
    inverse = [row[n:] for row in a]
    for i in range(n):
        for j in range(n):
            if sum(matrix[i][k]*inverse[k][j] for k in range(n)) % p != int(i == j):
                raise AssertionError('modular inverse replay')
    return inverse, determinant % p, pivots


def convolution(a, b, order, p):
    return [sum(a[j]*b[k-j] for j in range(k+1)) % p for k in range(order+1)]


def compile_polynomials(rows, variables, constants, t_name, p):
    index = {name: i for i, name in enumerate(variables)}
    result = []
    for row in rows:
        terms = []
        for monomial, coefficient in row.terms.items():
            value = rational_mod(coefficient, p)
            delay = 0
            powers = []
            for name, exponent in monomial:
                if name == t_name:
                    delay = exponent
                elif name in index:
                    if exponent < 0:
                        raise ValueError('negative jet variable exponent')
                    powers.append((index[name], exponent))
                else:
                    value = value * pow(constants[name], exponent, p) % p
            if delay < 0:
                raise ValueError('negative chart exponent')
            if value:
                terms.append((value, delay, tuple(powers)))
        result.append(terms)
    return result


def evaluate_series(rows, jets, order, p):
    one = [1] + [0]*order
    powers = {}
    monomials = {(): one}
    result = []
    for terms in rows:
        value = [0]*(order+1)
        for coefficient, delay, mon in terms:
            if delay > order:
                continue
            if mon not in monomials:
                product = one
                for index, exponent in mon:
                    key = (index, exponent)
                    if key not in powers:
                        power = one
                        for _ in range(exponent):
                            power = convolution(power, jets[index][:order+1], order, p)
                        powers[key] = power
                    product = convolution(product, powers[key], order, p)
                monomials[mon] = product
            term = monomials[mon]
            for k in range(delay, order+1):
                value[k] = (value[k] + coefficient*term[k-delay]) % p
        result.append(value)
    return result
