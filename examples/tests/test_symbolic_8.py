# -*- coding: utf-8 -*-

from quantarhei.symbolic.cumulant import Uged, Uegd, Ugde,Uedg, ExpdV
from quantarhei.symbolic.cumulant import gg, g1, g2
from quantarhei.symbolic.cumulant import CumulantExpr
from quantarhei.symbolic.abc import a, b, c, d, e, t, T, tau, x, y
from quantarhei.symbolic.lang import python_code, fortran_code

from sympy import S, Symbol
from sympy import sympify, collect
from sympy import diff
from sympy import exp

t1 = Symbol('t1')
t2 = Symbol('t2')
t3 = Symbol('t3')


def evaluate_cumulant(cum):
    """
    
    """
    A = cum.rewrite(gg)
    expr = CumulantExpr(A)
    expr = expr.evaluate()
    expr = CumulantExpr(expr)._make_positive(t1)
    D = expr._leading_index(a)
    expr = D._getExpr()

    ss = fortran_code(expr.__str__())
    
    return ss

    
def R1g():
    """

    """
    A = Ugde(b,t1)*Uedg(b,t1+t2)*Ugde(a,t1+t2+t3)

    return evaluate_cumulant(A)

def R2g():
    """

    """
    A = Uedg(a,t1+t2)*Ugde(b,t1+t2+t3)*Uedg(b,t1)

    return evaluate_cumulant(A)

               
print("R1g:")
print(R1g())

print("")
print("R2g:")
print(R2g())

print("")
print("R3g:")


