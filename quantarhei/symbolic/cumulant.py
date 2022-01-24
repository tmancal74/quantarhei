# -*- coding: utf-8 -*-
from sympy.physics.quantum import Operator, Dagger
from sympy.physics.quantum.qexpr import QExpr 
from sympy import I, conjugate
from sympy import S
from sympy import Function, Wild, Mul, Pow
from sympy import sympify


class CumulantException(Exception):
    pass


"""
    
    Cumulant expression class


"""
class CumulantExpr(QExpr):

        
    def _eval_simplify(self, ratio, measure, rational, inverse):
        return self._evaluate_second_order_rules()
        #return A._getExpr()
        

    def getOrder(self,n):
        if n == self._calculate_order(self):
            return self
        else:
            add = self.args[0]
            A = sympify(0)
            for aa in add.args:
                if self._calculate_order(aa) == 2 :
                    A = A + aa

            return CumulantExpr(A)
            
            
    def evaluate(self,large=False):
        expr = self.expand()
        expr = expr.getOrder(2)
        C = expr.simplify()
        if large:
            C = C._make_positive(large)
            C = C._expand_in_large(large)
        expr = C._getExpr()
        return expr.simplify()        
        
    """ 
    ***************************************************************************
    Helper routines of Cumulant Expr
    ***************************************************************************
    """

        

        

    def _evaluate_second_order_rules(self):   
        """
        Evaluates second order terms in terms of the line shape function 
        
        """
        a = Wild('a')
        b = Wild('b')
        w = Wild('w')
        t1 = Wild('t1')
        t2 = Wild('t2')
        
        """
        Combinations of hh
        
        hh_plus * hh_plus --->         
        hh_minus * hh_minus --->        
        hh_minus**2 --->         
        hh_plus * hh_minus --->         
        hh_minus * hh_plus --->        
        """
        A = self.replace(w*hh_plus(a,t1)*hh_plus(b,t2), \
                         w*(gg(a,b,t1)-gg(a,b,t1-t2)+conjugate(gg(a,b,t2))))        
        A = A.replace(w*hh_plus(a,t1)**2,w*(gg(a,a,t1)+conjugate(gg(a,a,t1)))) 
        A = A.replace(w*hh_minus(a,t1)*hh_minus(b,t2), \
                         w*(conjugate(gg(b,a,t1)) \
                           -conjugate(gg(b,a,t1-t2))+gg(b,a,t2)))        
        A = A.replace(w*hh_minus(a,t1)**2,w*(gg(a,a,t1)+conjugate(gg(a,a,t1))))       
        A = A.replace(w*hh_plus(a,t1)*hh_minus(b,t2), \
                         w*(-gg(a,b,t1) \
                           +gg(a,b,t1+t2)-gg(a,b,t2)))  
        A = A.replace(w*hh_minus(a,t1)*hh_plus(b,t2), \
                         w*(conjugate(-gg(b,a,t1) \
                           +gg(b,a,t1+t2)-gg(b,a,t2))))                         
        
        """        
        Replacement rules for ggs

        First daggered ggs
        (
        and that the normal ones
        
        gg_plus ---> gg 
        \hat{g}^{(-)}_{ab}(t) ---> g^{*}_{ab}(t)
        
        """
        A = A.replace(w*Dagger(gg_plus(a,b,t1)),w*conjugate(gg(a,b,t1)))
        A = A.replace(w*Dagger(gg_minus(a,b,t1)),w*gg(a,b,t1))
        A = A.replace(w*gg_plus(a,b,t1),w*gg(a,b,t1))
        A = A.replace(w*gg_minus(a,b,t1),w*conjugate(gg(a,b,t1)))
               
        """
        Replacement rules for dVs and their combinations with hh
        """
        A = A.replace(w*dV(a,t1)*dV(b,t2),w*g2(a,b,t1-t2))
        A = A.replace(w*dV(a,t1)**2,w*g2(a,a,0))
        
        A = A.replace(w*dV(a,t1)*hh_plus(b,t2),w*(-g1(a,b,t1-t2)+g1(a,b,t1)))
        A = A.replace(w*dV(a,t1)*hh_minus(b,t2),w*(g1(a,b,t1+t2)-g1(a,b,t1)))
        A = A.replace(w*hh_plus(a,t1)*dV(b,t2), \
                      w*(g1(a,b,t1-t2)+conjugate(g1(b,a,t2))))
        A = A.replace(w*hh_minus(a,t1)*dV(b,t2),
                      w*conjugate(g1(b,a,t1+t2)-g1(b,a,t2)))        
        #A = A.replace(w*hh_plus(a,t1)*dV(b,t2), \
        #              w*(g1(a,b,t1-t2)-conjugate(g1(a,b,t2))))
        #A = A.replace(w*hh_minus(a,t1)*dV(b,t2),
        #              w*conjugate(g1(a,b,t1+t2)-g1(a,b,t2)))        
        
        return A
       
    def _make_positive(self,arg):
        a = Wild('a')
        b = Wild('b')
        w = Wild('w')
        t1 = Wild('t')
        A = self.replace(w*gg(a,b,-arg+t1),w*conjugate(gg(b,a,arg-t1)))
        A = A.replace(w*gg(a,b,-arg),w*conjugate(gg(b,a,arg)))    
        A = A.replace(w*g1(a,b,-arg+t1),-w*conjugate(g1(b,a,arg-t1)))
        A = A.replace(w*g1(a,b,-arg),-w*conjugate(g1(b,a,arg))) 
        A = A.replace(w*g2(a,b,-arg+t1),w*conjugate(g2(b,a,arg-t1)))    
        A = A.replace(w*g2(a,b,-arg),w*conjugate(g2(b,a,arg)))
        return A        
        
    def _expand_in_large(self,arg):
        a = Wild('a')
        b = Wild('b')
        w = Wild('w')
        t1 = Wild('t1')
        A = self.replace(w*gg(a,b,arg+t1),w*(gg(b,a,arg)+(dd(a,b)-I*lam(a,b))*t1))
        A = A.replace(w*g1(a,b,arg+t1),w*(dd(a,b)-I*lam(a,b)))
        A = A.replace(w*g2(a,b,arg+t1),0)
        return A
        
    def _leading_index(self,arg):
        a = Wild('a')
        w = Wild('w')
        t1 = Wild('t1')
        A = self.replace(w*gg(a,arg,t1),w*gg(arg,a,t1))  
        A = A.replace(w*g1(a,arg,t1),w*g1(arg,a,t1))  
        A = A.replace(w*g2(a,arg,t1),w*g2(arg,a,t1))  
        A = A.replace(w*dd(a,arg),w*dd(arg,a))
        A = A.replace(w*lam(a,arg),w*lam(arg,a))
        return A
        
    def _eliminate_off_diagonal(self):
        return self

    def _getExpr(self):
        return self.args[0]
        
    def _calculate_order(self,expr):
        """ 
        Calculates a perturbation order of the expression expr 
        
        """
        content = expr
        order = 0
            
        if  content.func is Mul:
            for aa in content.args:
                order += self._calculate_order(aa)    
            return order

        elif content.func is Pow:   
            sorder = self._calculate_order(content.args[0])
            return sorder*content.args[1]

        elif content.func in (hh_plus,hh_minus,gg_plus,gg_minus,dV):            
            return content.order()
            
        elif content.func is Dagger:
            return self._calculate_order(content.args[0])
        else:
            return 0
            
       
""" 
    Special Operators

"""
class Ugde(Operator):
    nargs = 2
    def _eval_rewrite_as_gg(self,a,t):
        return (1-I*hh_plus(a,t)-gg_plus(a,a,t))

class Uedg(Operator):
    nargs = 2
    def _eval_rewrite_as_gg(self,a,t):
        return (1+I*hh_plus(a,t)-Dagger(gg_plus(a,a,t)))

class Uged(Operator):
    nargs = 2
    def _eval_rewrite_as_gg(self,a,t):
        return (1+I*hh_minus(a,t)-gg_minus(a,a,t))
        
class Uegd(Operator):
    nargs = 2
    def _eval_rewrite_as_gg(self,a,t):
        return (1-I*hh_minus(a,t)-Dagger(gg_minus(a,a,t)))

class ExpdV(Operator):
    nargs = 3
    def _eval_rewrite_as_gg(self,a,t,x):
        return (1+x*dV(a,t)+x**2*dV(a,t)*dV(a,t))

class hh_plus(Operator):
    
    def order(self):
        return 1
        
class hh_minus(Operator):

    def order(self):
        return 1

class gg_plus(Operator):

    def order(self):
        return 2

class gg_minus(Operator):

    def order(self):
        return 2

class dV(Operator):
    
    def order(self):
        return 1

"""

    Lineshape function and related stuff

"""

class gg(Function):
    nargs = (1,2,3)
       
    @classmethod
    def eval(cls, a, b, x):
        if x.is_Number:
            if x is S.Zero:
                return S.Zero
        #if not (a is b):
        #    return S.Zero
""" First derivative of gg """                
class g21(Function):
    nargs = (1,2,3)
                
    @classmethod
    def eval(cls,a,b,x):
        if x.is_Number:
            if x is S.Zero:
                return S.Zero
""" First derivative of gg """                
class g1(Function):
    nargs = (1,2,3)
                
    @classmethod
    def eval(cls,a,b,x):
        if x.is_Number:
            if x is S.Zero:
                return S.Zero


 

                
""" Second derivative of gg """                
class g2(Function):
    nargs = (1,2,3)
                

""" dephasing rate """    
class dd(Function):
    nargs = (1,2)

    def _eval_is_real(self):
        return True            

""" Reorganization energy """
class lam(Function):
    nargs = (1,2)

    def _eval_is_real(self):
        return True            



