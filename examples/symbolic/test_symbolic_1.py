# -*- coding: utf-8 -*-
from quantarhei.symbolic.cumulant import Uged, Ugde, Uedg, Uegd
from quantarhei.symbolic.cumulant import gg
from quantarhei.symbolic.cumulant import CumulantExpr
from quantarhei.symbolic.abc import a, b, c, d, t, T


""" 
Test of cumulant expansion method on a first order term of 
modified 

 
    <a|H(t)|b><c|W|d>
    
    = J_ac*<\Psi_d|Dagger(U_a(t))*U_b(t)|\Psi_c>/Norm
    = J_ac*<\Psi_g|U_g(T)Dagger(U_d(T))Dagger(U_a(t))U_g(t)
     x Dagger(U_g(t))U_b(t)U_c(T)Dagger(U_g(T))|\Psi_g>/Norm
    = J_ac*<\Psi_g|[U_g(T)Dagger(U_d(T))][Dagger(U_a(t))U_g(t)]
     x [Dagger(U_g(t))U_b(t)][U_c(T)Dagger(U_g(T))]|\Psi_g>
     x Norm^-1

    => Uged(d,T)*Uedg(a,t)*Ugde(b,t)*Uegd(c,T)*(1/Norm)
    
    Norm = Uged(d,T)*Uegd(c,T)
        
    
"""
A     = Uged(d,T) *Uedg(a,t)*Ugde(b,t)* Uegd(c,T)

Anorm = Uged(d,T)*Uegd(c,T)


verbatim = True
if verbatim:
    print(" ")
    print("Expression to evaluate: ")
    print(" ")
    print("Tr_bath{",A,"W_eq}")
    print(" ")
    print(" ")


A = A.rewrite(gg)
expr = CumulantExpr(A)
""" use option large=T to calculate evaluate in T --> oo """
expr = expr.evaluate(large=T) 
""" use the symetry of lineshape function in the exciton indices """
#D = CumulantExpr(expr)._leading_index(a)
#expr = D._getExpr()

A = Anorm.rewrite(gg)
norm = CumulantExpr(A)
""" use option large=T to calculate evaluate in T --> oo """
norm = norm.evaluate(large=T) 
""" use the symetry of lineshape function in the exciton indices """
#D = CumulantExpr(expr)._leading_index(a)
#expr = D._getExpr()

expr = (expr-norm).simplify()

if verbatim:
    print("Cumulant: ")
    print(" ")
    print(expr)
    print(" ")
    print("In simplified form: ")
    print(" ")

expr = expr.subs(c,b)
expr = expr.subs(d,b).simplify()
print(expr)
