# -*- coding: utf-8 -*-
from quantarhei.symbolic.cumulant import Uged, Ugde, Uedg, Uegd
from quantarhei.symbolic.cumulant import gg
from quantarhei.symbolic.cumulant import CumulantExpr
from quantarhei.symbolic.abc import a, b, c, d, e, t, T, tau

from sympy import Symbol

""" 
Test of cumulant expansion method on the second order term of non-secular
Foerster equation.

    <a|H(t)|b><c|H(t-tau)|d><d|W|e>
    
    = <\Psi_e|<a|H(t)|b><c|H(t-tau)|d>|\Psi_d>
    = <\Psi_g|Dagger(U_e(T))Dagger(U_a(t))U_b(t)Dagger(U_c(t-tau))U_d(t-tau)
      x U_d(T)|\Psi_g>
      
    = <\Psi_g|[U_g(T)Dagger(U_e(T))][Dagger(U_a(t))U_g(t)][Dagger(U_g(t)U_b(t)]
    x [Dagger(U_c(t-tau))U_g(t-tau)][Dagger(U_g(t-tau)U_d(t-tau)]
    x [U_d(T)Dagger(U_g(T))]|\Psi_g>
    
    = Uged(e,T)*Uedg(a,t)*Ugde(b,t)*Uedg(c,t-tau)*Ugde(d,t-tau)*Uegd(d,T)

"""


A = Uged(e,T)*Uedg(a,t)*Ugde(b,t)*Uedg(c,t-tau)*Ugde(d,t-tau)*Uegd(d,T)
Anorm = Uged(e,T)*Uegd(d,T)  

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
    print("Cumulant: ")
    print(" ")
    print(expr)
    print(" ")

print("Foerster rate")
expr = expr.subs(e,d).simplify()

a1 = Symbol('a1')
b1 = Symbol('b1')
expr = expr.subs({d:b1,c:a1,b:a1,a:b1})
expr = expr.subs({a1:a,b1:b}).simplify()

print(expr)

