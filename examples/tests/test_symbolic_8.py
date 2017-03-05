# -*- coding: utf-8 -*-
"""
Calculation of cumulant expressions for non-linear response functions
of the third order for a multilevel three band system.



"""
from quantarhei.symbolic.cumulant import  Ugde, Uedg #, Uegd, Uged, ExpdV
from quantarhei.symbolic.cumulant import gg #, g1, g2
from quantarhei.symbolic.cumulant import CumulantExpr
from quantarhei.symbolic.abc import a, b, f, tau #, c, d, e, t, T, tau, x, y
from quantarhei.symbolic.abc import t1, t2, t3
#from quantarhei.symbolic.lang import python_code
from quantarhei.symbolic.lang import fortran_code


def evaluate_cumulant(cum):
    """
    
    """
    A = cum.rewrite(gg)
    expr = CumulantExpr(A)
    expr = expr.evaluate()
    expr = CumulantExpr(expr)._make_positive(t1)
    expr = CumulantExpr(expr)._make_positive(t2)    
    expr = CumulantExpr(expr)._make_positive(t3)    
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

def R3g():
    """
    
    """
    A = Uedg(a,t1)*Ugde(b,t1+t2+t3)*Uedg(b,t1+t2)

    return evaluate_cumulant(A)


def R4g():
    """
    
    """
    A = Ugde(b,t1+t2+t3)*Uedg(b,t1+t2)*Ugde(a,t1)

    return evaluate_cumulant(A)

def R1fs():
    """
    
    """
    A = (Uedg(a,t1+t2+t3)*Ugde(f,t1+t2+t3)*Uedg(f,t1+t2)
        *Ugde(b,t1+t2)*Uedg(b,t1))

    return evaluate_cumulant(A)

def R2fs():
    """
    
    """
    A = (Ugde(b,t1)*Uedg(b,t1+t2+t3)*Ugde(f,t1+t2+t3)
        *Uedg(f,t1+t2)*Ugde(a,t1+t2))

    return evaluate_cumulant(A)


def print_R1gt():
    """
    
    """    
    A = Ugde(b,t3)
    print(evaluate_cumulant(A))
    
    B = Ugde(a,t1)
    print(evaluate_cumulant(B))
    
def print_R2gt():
    """
    
    """    
    A = Ugde(b,t3)
    print(evaluate_cumulant(A))
    
    B = Uedg(a,t1)
    print(evaluate_cumulant(B))    

def print_R1fst():
    """
    
    """    
    A = Uedg(b,t3)*Ugde(f,t3)
    print(evaluate_cumulant(A))
    
    B = Uedg(a,t1)
    print(evaluate_cumulant(B))    
    
 
def print_R2fst():
    """
    
    """    
    A = Uedg(b,t3)*Ugde(f,t3)
    print(evaluate_cumulant(A))
    
    B = Ugde(a,t1)
    print(evaluate_cumulant(B)) 
    
def print_trans_R2g():
    """

    """
    A = (Uedg(a,t1+tau)*Ugde(b,t1+tau)*Uedg(b,t1+t2)*Ugde(b,t1+t2+t3)
        *Uedg(b,t1+tau)*Ugde(a,t1+tau)*Uedg(a,t1))

    print(evaluate_cumulant(A))

       
    
print("R1g:")
print(R1g())

print("")
print("R2g:")
print(R2g())

print("")
print("R3g:")
print(R3g())

print("")
print("R4g:")
print(R4g())

print("")
print("R1fs:")
print(R1fs())

print("")
print("R2fs:")
print(R2fs())

print("")
print("R1gt")
print_R1gt()

print("")
print("R2gt")
print_R2gt()

print("")
print("R1fst")
print_R1fst()

print("")
print("R2fst")
print_R2fst()

print("")
print("Trans_R2g")
print_trans_R2g()
