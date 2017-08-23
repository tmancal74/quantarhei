# -*- coding: utf-8 -*-

from quantarhei.symbolic.cumulant import Uged, Uegd, Ugde,Uedg, ExpdV
from quantarhei.symbolic.cumulant import gg, g1, g2
from quantarhei.symbolic.cumulant import CumulantExpr
from quantarhei.symbolic.abc import a, b, c, d, e, t, T, tau, x, y

from sympy import S, Symbol
from sympy import sympify, collect
from sympy import diff
from sympy import exp

""" 
Test of cumulant expansion method on the second order term of non-secular
Modified Redfield equation.

    <a|H(t)|b><c|H(t-tau)|d><d|W|e>
    
    = <\Psi_e|<a|H(t)|b><c|H(t-tau)|d>|\Psi_d>
    = <\Psi_g|Dagger(U_e(T))Dagger(U_a(t))<a|dV|b>U_b(t)Dagger(U_c(t-tau))
      <c|dV|d>U_d(t-tau)U_d(T)|\Psi_g>
      
    = <\Psi_g|[U_g(T)Dagger(U_e(T))][Dagger(U_a(t))U_g(t)][Dagger(U_g(t)U_b(t)]
    x [Dagger(U_c(t-tau))U_g(t-tau)][Dagger(U_g(t-tau)U_d(t-tau)]
    x [U_d(T)Dagger(U_g(T))]|\Psi_g>
    
    = Uged(e,T)*Uedg(a,t)*Ugde(b,t)*Uedg(c,t-tau)*Ugde(d,t-tau)*Uegd(d,T)

"""

m = Symbol('m')
n = Symbol('n')
t1 = Symbol('t1')

A = Uged(n,T)*Uedg(n,tau)*ExpdV(a,tau,x)*Ugde(m,tau)*ExpdV(b,0,y)*Uegd(n,T)
A = Uged(b,T)*ExpdV(n,0,x)*Uegd(a,tau)*ExpdV(m,-tau,y)*Uged(b,tau)*Uegd(b,T)

Anorm = Uged(n,T)*Uegd(n,T)  

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
""" use option large=T to  evaluate in T --> oo """
expr = expr.evaluate(large=T)
expr = CumulantExpr(expr)._make_positive(tau)
""" use the symetry of lineshape function in the exciton indices """
D = expr._leading_index(a)
expr = D._getExpr()

A = Anorm.rewrite(gg)
norm = CumulantExpr(A)# -*- coding: utf-8 -*-



""" use option large=T to evaluate in T --> oo """
norm = norm.evaluate(large=T) 
""" use the symetry of lineshape function in the exciton indices """
#D = CumulantExpr(expr)._leading_index(a)
#expr = D._getExpr()
print("Norm: ", norm)

expr = (expr-norm).subs(t1,0).simplify()

info = False

""" 
Test of cumulant expansion method

"""
#A = Uedg(b,t)*ExpdV(a,tau,x)*Ugde(a,tau)*ExpdV(b,t,y)
if info:
    print("Cumulant exponent: ")
    print(" ")
    print(expr)
    print(" ")
    print(" ")
    
    print("Analyzing the cumulant exponent and sorting parameters")
    print(" ")
    terms = collect(expr,[x,y],evaluate=False)
    print("   1:")

    A1 = sympify(terms[S.One])
    print(A1)
    print("   x:")
    Ax1 = sympify(terms[x])
    Ax = Ax1.subs(y,0)
    print(Ax)
    print("   y:")
    Ay1 = sympify(terms[y])
    Ay = Ay1.subs(x,0)
    print(Ay)
    print(" x*y:")
    Axy1 = collect(Ax1,y,evaluate=False)
    Axy = sympify(Axy1[y])
    print(Axy)
    print("x**2:")
    Ax2 = sympify(terms[x**2])
    print(Ax2)
    print(" ")
    print(" ")

if verbatim:     
    print("Cumulant: ")
    print(" ")
    print(expr)
    print(" ")

#expr = expr.subs(c,b)
#expr = expr.subs(d,b).simplify()
#print(expr)

print(" ")
print("Final form: ")
B = diff(diff(exp(expr),x),y).subs({x:0,y:0})
print(" ")
print(B.simplify())
