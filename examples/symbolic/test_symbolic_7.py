# -*- coding: utf-8 -*-

from quantarhei.symbolic.cumulant import Uged, Uegd, Ugde,Uedg, ExpdV
from quantarhei.symbolic.cumulant import gg, g1, g2
from quantarhei.symbolic.cumulant import CumulantExpr
from quantarhei.symbolic.lang import python_code, fortran_code
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
t2 = Symbol('t2')
t3 = Symbol('t3')
t4 = Symbol('t4')
f = Symbol('f')
k = Symbol('k')
k1 = Symbol('k1')

#A = Uged(n,T)*Uedg(n,tau)*ExpdV(a,tau,x)*Ugde(m,tau)*ExpdV(b,0,y)*Uegd(n,T)
#A = Uged(b,T)*ExpdV(n,0,x)*Uegd(a,tau)*ExpdV(m,-tau,y)*Uged(b,tau)*Uegd(b,T)

#A = Uedg(b,t)*ExpdV(m,0,x)*Ugde(a,t)*Uedg(a,tau)*ExpdV(n,0,y)*Ugde(b,tau)
#A = Uedg(b,t)*ExpdV(m,t,x)*Ugde(a,t)*Uedg(a,tau)*ExpdV(n,tau,y)*Ugde(b,tau)
#A = Uedg(b,t)*             Ugde(a,t-tau)*             Ugde(b,tau)

#A = Uedg(b,tau)*ExpdV(m,tau,x)*Ugde(a,tau)*Uedg(a,0)*ExpdV(n,0,y)*Ugde(a,0)
#A = Uedg(b,t)*ExpdV(m,t,x)*Ugde(a,t)*Uedg(a,tau)*ExpdV(n,tau,y)*Ugde(a,tau)

#A = Uged(e,T)*Uedg(a,t)*Ugde(b,t)*Uedg(c,t-tau)*Ugde(d,t-tau)*Uegd(d,T)


""" Parts of the relaxation tensor """
part = "test3"

""" M1 K """
if part == "M1_K":
    A = Uged(b,T)*Uedg(b,t)*ExpdV(m,t,x)*Ugde(a,t)*Uedg(a,t-tau)\
    *ExpdV(n,t-tau,y)*Ugde(b,t-tau)*Uegd(b,T)

    use_norm = False
    Anorm = Uged(n,T)*Uegd(n,T)  
    nderiv = 2

    filename = "m1_K.txt"

elif part == "M1_H":
    
    pass


elif part == "K_bk":
    
    A = Uged(a,t)*ExpdV(m,0,x)*Uegd(c,t)
    nderiv = 1

    filename = "K_bk.txt"

elif part == "norm":
    
    A = Uged(a,t)*Uegd(b,t)
    nderiv = 0

    filename = "norm"
    
elif part == "test1":
    nderiv = 2
    filename = "test1"
    A = ExpdV(a,tau,x)*ExpdV(a,0,y)*Uegd(a,t-tau)
    use_norm =True
    Anorm = Ugde(a,t-tau)

elif part == "test2":
    nderiv = 1
    filename = "test2"
    A =Uedg(a,t)*ExpdV(m,t,x)*Ugde(a,t) #m=ac
    use_norm =True
    Anorm = Uedg(a,t)*Ugde(a,t)
    
elif part == "test":
    nderiv = 1
    filename = "test"
    A = ExpdV(a,t,x)*Ugde(a,t)
    use_norm =True
    Anorm = Ugde(a,t)
    
elif part == "Kab,c":
    nderiv = 1
    filename = "Kabc"
    A = Uedg(b,t)*ExpdV(m,t,x)*Ugde(c,t) #m=ac
    use_norm =True
    Anorm = Uedg(b,t)*Ugde(c,t)

elif part == "Kc,ab":
    nderiv = 1
    filename = "Kcab"
    A = Uedg(c,t)*ExpdV(m,t,x)*Ugde(a,t) #m=cb
    use_norm =True
    Anorm = Uedg(c,t)*Ugde(a,t)
    
elif part == "Jeff":
    nderiv = 0
    filename = "jeff"
    use_norm = False
    A = Ugde(e,t1)*Uedg(e,t2)*Ugde(f,t2)*Uedg(f,t3)*Ugde(e,t3)*Uedg(e,t4)
   
    A = Ugde(e,t1)*Uedg(e,t2)*Ugde(e,t3)*Uedg(e,t4)
    
elif part == "test3":
    nderiv = 1
    filename = "jeff"
    use_norm = False
    A = Uedg(a,T-tau)*ExpdV(m,T-tau,x)*Ugde(c,T-tau)
    Anorm = Uedg(a,T)*ExpdV(m,T,x)*Ugde(c,T)

verbatim = False

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
expr = expr.evaluate() #large=T)
#expr = CumulantExpr(expr)._make_positive(t)
""" use the symetry of lineshape function in the exciton indices """
#D = expr._leading_index(b)
#D = D._leading_index(a)
#expr = D._getExpr()



if use_norm:
    A = Anorm.rewrite(gg)
    norm = CumulantExpr(A)

    """ use option large=T to evaluate in T --> oo """
    norm = norm.evaluate(large=T) 
    """ use the symetry of lineshape function in the exciton indices """
    #D = CumulantExpr(expr)._leading_index(a)
    #expr = D._getExpr()
    print("Norm: ", norm)

    expr = (expr-norm).subs(t1,0).simplify()
    
else:

    #expr = expr.subs(t1,0).simplify()    
    pass
    

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
if nderiv == 2:
    B = diff(diff(exp(expr),x),y).subs({x:0,y:0})
elif nderiv == 1:
    B = diff(exp(expr),x).subs({x:0})
elif nderiv == 0:
    B = expr
else:
    raise Exception("nderiv can only be 0, 1 or 2")

print(" ")

final_form = B.simplify()
ss = final_form.__str__()
print(ss)
print(" ")

arrays = ["gg","g1","g2"]
sr = python_code(ss,arrays)

print(" ")
print("Python code:")
print(" ")
print(sr)

f = open(filename,'wt')
f.write(sr)
f.close()

sr = fortran_code(ss,)

print(" ")
print("Fortran code:")
print(" ")
print(sr)