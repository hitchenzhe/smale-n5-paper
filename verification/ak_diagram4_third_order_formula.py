"""Table 5: exact complex second obstruction from first/third-order vectors."""
from ak_diagram4_first_order_formula import vectors, csub, cmul, cpow
def second_obstruction(signs,m,g4,g5):
    X1,Y1=vectors(signs,m,g4,g5,1)
    X3,Y3=vectors(signs,m,g4,g5,3)
    return csub(cmul(Y1,cpow(X3,3)),cmul(X1,cpow(Y3,3)))
