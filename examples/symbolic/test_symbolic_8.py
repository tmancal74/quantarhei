# -*- coding: utf-8 -*-
"""
Calculation of cumulant expressions for non-linear response functions
of the third order for a multilevel three band system.



"""
from quantarhei.symbolic.cumulant import  Ugde, Uedg, Uged, Uegd #, ExpdV
from quantarhei.symbolic.cumulant import gg #, g1, g2
from quantarhei.symbolic.cumulant import CumulantExpr
from quantarhei.symbolic.abc import a, b, f, tau, tau1, tau2, tau3, c, d #, e, t, T, tau, x, y
from quantarhei.symbolic.abc import t1, t2, t3
from quantarhei.symbolic.lang import python_code
from quantarhei.symbolic.lang import fortran_code

import time


def evaluate_cumulant(cum, positive_times = [], leading_index=None,
                      lang = "Python", arrays=None):
    """
    
    """
    t0 = time.time()
    A = cum.rewrite(gg)
    expr = CumulantExpr(A)
    expr = expr.evaluate()
    
    t1 = time.time()
    for tt in positive_times:
        expr = CumulantExpr(expr)._make_positive(tt)    
        
    t2 = time.time()
    #a = leading_index[0]
    if leading_index is not None:
        D = expr._leading_index(leading_index)
        expr = D._getExpr()
        
    t3 = time.time()
    if lang == "Fortran":
        ss = fortran_code(expr.__str__())
    elif lang == "Python":
        ss = python_code(expr.__str__(),arrays=arrays)
    else:
        raise Exception("Unknown language")
    
    print(t1-t0)
    print(t2-t1)
    print(t3-t2)
    
    return ss

    
def R1g():
    """

    """
    A = Ugde(b,t1)*Uedg(b,t1+t2)*Ugde(a,t1+t2+t3)

    return evaluate_cumulant(A, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"])

def R2g():
    """

    """
    A = Uedg(a,t1+t2)*Ugde(b,t1+t2+t3)*Uedg(b,t1)

    return evaluate_cumulant(A, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"])

def R3g():
    """
    
    """
    A = Uedg(a,t1)*Ugde(b,t1+t2+t3)*Uedg(b,t1+t2)

    return evaluate_cumulant(A, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"])


def R4g():
    """
    
    """
    A = Ugde(b,t1+t2+t3)*Uedg(b,t1+t2)*Ugde(a,t1)

    return evaluate_cumulant(A, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"])

def R1fs():
    """
    
    """
    A = (Uedg(a,t1+t2+t3)*Ugde(f,t1+t2+t3)*Uedg(f,t1+t2)
        *Ugde(b,t1+t2)*Uedg(b,t1))

    return evaluate_cumulant(A, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"])

def R2fs():
    """
    
    """
    A = (Ugde(b,t1)*Uedg(b,t1+t2+t3)*Ugde(f,t1+t2+t3)
        *Uedg(f,t1+t2)*Ugde(a,t1+t2))

    return evaluate_cumulant(A, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"])


def print_R1gt():
    """
    
    """    
    A = Ugde(b,t3)
    print(evaluate_cumulant(A, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"]))
    
    B = Ugde(a,t1)
    print(evaluate_cumulant(B, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"]))
    
def print_R2gt():
    """
    
    """    
    A = Ugde(b,t3)
    print(evaluate_cumulant(A, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"]))
    
    B = Uedg(a,t1)
    print(evaluate_cumulant(B, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"]))    

def print_R1fst():
    """
    
    """    
    A = Uedg(b,t3)*Ugde(f,t3)
    print(evaluate_cumulant(A, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"]))
    
    B = Uedg(a,t1)
    print(evaluate_cumulant(B, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"]))    
    
 
def print_R2fst():
    """
    
    """    
    A = Uedg(b,t3)*Ugde(f,t3)
    print(evaluate_cumulant(A, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"]))
    
    B = Ugde(a,t1)
    print(evaluate_cumulant(B, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"])) 
    
def print_trans_R2g():
    """

    """
    A = (Uedg(a,t1+tau)*Ugde(b,t1+tau)*Uedg(b,t1+t2)*Ugde(b,t1+t2+t3)
        *Uedg(b,t1+tau)*Ugde(a,t1+tau)*Uedg(a,t1))

    print(evaluate_cumulant(A, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"]))

    
def print_trans_R2g_alt():
    """

    """
    #A = (Uedg(a,t1+tau)*Ugde(b,t1+tau)*Uedg(b,t1+t2)*Ugde(b,t1+t2+t3)
    #    *Uedg(b,t1+tau)*Ugde(a,t1+tau)*Uedg(a,t1))
    
    A = (Uged(a,t1)*Uedg(a,tau1)*Ugde(b,tau1)*Uedg(b,t2)*Ugde(b,t2+t3)*Uedg(b,tau1)*Ugde(a,tau1))

    print(evaluate_cumulant(A, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"]))    
    
    
def print_trans_R2g_alt2():
    """

    """
    #A = (Uedg(a,t1+tau)*Ugde(b,t1+tau)*Uedg(b,t1+t2)*Ugde(b,t1+t2+t3)
    #    *Uedg(b,t1+tau)*Ugde(a,t1+tau)*Uedg(a,t1))
    
    #A = (Uged(a,t1)*Uedg(a,tau1)*Ugde(b,tau1)*Uedg(b,t2)*Ugde(b,t2+t3)*Uedg(b,tau1)*Ugde(a,tau1))
    
    A = (Uged(a,t1+tau1)*Uedg(b,t2-tau1)*Ugde(b,t2+t3-tau1)*Uegd(a,tau1))

    print(evaluate_cumulant(A, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"]))    
    
    
def generate_nth_order_R2g(states_tuple, times_tuple):
    
    order = len(states_tuple)
    if order != len(times_tuple):
        raise Exception("Wrong tuple/list length")
        
    # starting state
    a = states_tuple[0]
    # final state (can be the same as starting)
    b = states_tuple[len(states_tuple)-1]

    # final time (must be t2)
    tt = times_tuple[len(times_tuple)-1]
    
    AL = Uged(a,t1)
    Amid = Uedg(b,tt)*Ugde(b,t3+tt)
    
    filL = 1
    filR = 1
    
    for k in range(len(times_tuple)-1):
        tau = times_tuple[k]
        s1 = states_tuple[k]
        s2 = states_tuple[k+1]
        filL = filL*Uedg(s1,tau)*Ugde(s2,tau)
        filR = Uedg(s2,tau)*Ugde(s1,tau)*filR
        
    
    A = AL*filL*Amid*filR
         
    print(A)

    print(evaluate_cumulant(A, positive_times=(t1, tt, t3),
                             leading_index=a, arrays=["gg"]))   
    
    
def test():
    
    A = Uged(a,t1+t2)*Ugde(d,t3)*Uegd(a,t2)

    print(evaluate_cumulant(A, positive_times=(t1, t2, t3),
                             leading_index=a, arrays=["gg"]))    
    
def oneex_twoex():
    
    A = Uedg(f,t1)*Ugde(a,t1)
    
    print(evaluate_cumulant(A, positive_times=(t1,), leading_index=a,
                            arrays="gg"))
    
    
    
# =============================================================================
# print("R1g:")
# st_R1g = "numpy.exp("+R1g()+")"
# print(st_R1g)
# 
# print("")
# print("R2g:")
# print(R2g())
# 
# print("")
# print("R3g:")
# print(R3g())
# 
# print("")
# print("R4g:")
# print(R4g())
# 
# print("")
# print("R1fs:")
# print(R1fs())
# 
# print("")
# print("R2fs:")
# print(R2fs())
# 
# print("")
# print("R1gt")
# print_R1gt()
# 
# print("")
# print("R2gt")
# print_R2gt()
# 
# print("")
# print("R1fst")
# print_R1fst()
# 
# print("")
# print("R2fst")
# print_R2fst()
# 
# =============================================================================
#print("")
#print("Trans_R2g")
#print_trans_R2g()
#
#print("")
#print("Trans_R2g_alt")
#print_trans_R2g_alt()
#
#print("")
#print("Trans_R2g_alt2")
#print_trans_R2g_alt2()

#print("***")
#states = (a, c, b) #(a,c,b)
#times = (tau1, tau2, t2) # (tau1,tau2,t2)
#generate_nth_order_R2g(states, times)
#
#print("===")
#A = Uged(a,t1)*Uedg(a,tau1)*Ugde(c,tau1)*Uedg(c,tau2)*Ugde(b,tau2)*Uedg(b,t2)*Ugde(b,t2 + t3)*Uedg(b,tau2)*Ugde(c,tau2)*Uedg(c,tau1)*Ugde(a,tau1)
#
#print(evaluate_cumulant(A, positive_times=(t1, t2, t3),
#                             leading_index=a, arrays=["gg"]))   

#print("***")
#states = (a,b,c, d) #(a,c,b)
#times = (tau1, tau2, tau3, t2) # (tau1,tau2,t2)
#states = (a,c,b)
#times = (tau1,tau2,t2)
#generate_nth_order_R2g(states, times)
#test()

oneex_twoex()