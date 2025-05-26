# -*- coding: utf-8 -*-
from quantarhei.symbolic.cumulant import Ugde, Uedg, ExpdV
#from quantarhei.symbolic.cumulant import Uegd, Uged
from quantarhei.symbolic.cumulant import gg
from quantarhei.symbolic.cumulant import CumulantExpr
from quantarhei.symbolic.abc import a, b, t, tau, x, y
#from quantarhei import stop

from sympy import diff
from sympy import exp

from sympy import Symbol

m = Symbol('m')
n = Symbol('n')
k = Symbol('k')

""" 
Cummulant expression for 





"""


A = Uedg(n,t)*ExpdV(a,t,x)*Ugde(m,t)*Uedg(m,t-tau)*ExpdV(b,t-tau,y)*Ugde(k,t-tau)

A1 = Uedg(n,t)*ExpdV(a,t,x)*Ugde(m,t)

A_nderiv = 2
A1_nderiv = 1

#Anorm = Uged(d,tau)*Uegd(c,tau)  

verbatim = True

if verbatim:
    print(" ")
    print("Expressions to evaluate: ")
    print(" ")
    print("    Tr_bath{",A,"W_eq}")
    print(" ")
    print(" ")
    print("    Tr_bath{",A1,"W_eq}")
    print(" ")
    #print("The expression is normalized by:")
    #print(" ")
    #print("    Tr_bath{",Anorm,"W_eq}")
    #print(" ")

A = A.rewrite(gg)
expr = CumulantExpr(A)
A1 = A1.rewrite(gg)
expr1 = CumulantExpr(A1)
""" use option large=T to evaluate in T --> oo """
expr = expr.evaluate() #large=T) 
expr1 = expr1.evaluate()
""" use the symetry of lineshape function in the exciton indices """
#D = CumulantExpr(expr)._leading_index(a)
#expr = D._getExpr()


print(" ")
print("Final form: ")
if A_nderiv == 2:
    B = diff(diff(exp(expr),x),y).subs({x:0,y:0})
elif A_nderiv == 1:
    B = diff(exp(expr),x).subs({x:0})
elif A_nderiv == 0:
    B = expr
else:
    raise Exception("nderiv can only be 0, 1 or 2")

print(" ")

final_form = B.simplify()
ss = final_form.__str__()
print(ss)

if A1_nderiv == 2:
    B1 = diff(diff(exp(expr1),x),y).subs({x:0,y:0})
elif A1_nderiv == 1:
    B1 = diff(exp(expr1),x).subs({x:0})
elif A1_nderiv == 0:
    B1 = expr1
else:
    raise Exception("nderiv can only be 0, 1 or 2")

print(" ")

final_form = B1.simplify()
ss = final_form.__str__()
print(ss)
print(" ")


                      
                      
                      
                      
                      
                      
                      
                      