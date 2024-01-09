# -*- coding: utf-8 -*-
from quantarhei.symbolic.cumulant import Uged, Ugde, Uedg, Uegd
from quantarhei.symbolic.cumulant import gg
from quantarhei.symbolic.cumulant import CumulantExpr
from quantarhei.symbolic.abc import a, b, c, d, e, t, T, tau

from sympy import Symbol

""" 
Cummulant expression for the time-non-local Förster theory

    U_a(t-tau)U_c(tau) W_eq Dagger(U_d(tau))Dagger(U_b(t-tau))
    
    =         Dagger(U_d(tau))   Dagger(U_b(t-tau))            U_a(t-tau) U_c(tau) W_eq
    
    =[U_g(tau)Dagger(U_d(tau))] [Dagger(U_b(t-tau))U_g(t-tau)]
     [Dagger(U_g(t-tau))U_a(t-tau)] [U_c(tau)Dagger(U_g(tau))]W_eq
    
    
    = Uged(d,tau) * Uedg(b,t-tau) * Ugde(a,t-tau) * Uegd(c,tau)

"""


A = Uged(d,tau)*Uedg(b,t-tau)*Ugde(a,t-tau)*Uegd(c,tau)
Anorm = Uged(d,tau)*Uegd(c,tau)  

verbatim = True

if verbatim:
    print(" ")
    print("Expression to evaluate: ")
    print(" ")
    print("    Tr_bath{",A,"W_eq}")
    print(" ")
    print("The expression is normalized by:")
    print(" ")
    print("    Tr_bath{",Anorm,"W_eq}")
    print(" ")


A = A.rewrite(gg)
expr = CumulantExpr(A)
""" use option large=T to evaluate in T --> oo """
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
    print("\nCumulant: ")
    print(" ")
    print(expr)
    print(" ")

print("Foerster rate")
#expr = expr.subs(e,d).simplify()

#a1 = Symbol('a1')
#b1 = Symbol('b1')
#expr = expr.subs({d:b1,c:a1,b:a1,a:b1})
#expr = expr.subs({a1:a,b1:b}).simplify()

print(expr)