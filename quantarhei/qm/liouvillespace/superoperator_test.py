# -*- coding: utf-8 -*-
"""
    SuperOperator class for test purposes
    
    You can take this superoperator and feed it to tests. There are several
    super operators which can be created distinguished by descriptive
    names. They are listed below:
    
    Test Names
    ----------
    
    dim-2-AOA :
        This superoperator corresponds to a multiplication by a 2x2 matrix
        from both left and write (A0A), where O is an operator
        
    dim-3-AOA :
        The same as dim-2-AOA only with 3x3 matrix
        
    
    Class Details
    -------------
    
"""
import numpy

from .superoperator import SuperOperator
import quantarhei as qr

class TestSuperOperator(SuperOperator):
    """SuperOperator class for test purposes
    
    
    Parameters
    ----------
    
    name : string
        Name of the test super operator to create
    
    
    Examples
    --------
    
    >>> import quantarhei as qr
    >>> S = qr.qm.TestSuperOperator("dim-2-AOA")
    >>> A = S._def_A(2)
    >>> C = S.apply(A)
    >>> D = numpy.dot(A,numpy.dot(A,A))
    >>> print(numpy.allclose(C,D))
    True
    
    """
    
    
    def __init__(self, name=None):
        
        if name is None:
            raise Exception("Name of the test super operator "+
                            "must be specified")
            
        if name == "dim-2-AOA":
            
            dim = 2
            
            super().__init__(dim=dim)
            self._AOA(dim)
            
        elif name ==  "dim-3-AOA":
            
            dim = 3
            
            super().__init__(dim=dim)
            self._AOA(dim)  
            
        else:
            raise Exception("Unknown test name")


    def _def_A(self, dim):
        """Defines the matrix A for constuction of a super operator
        
        """
        A = numpy.zeros((dim,dim), dtype=qr.COMPLEX)
        
        for k_i in range(dim):
            A[k_i,k_i] = k_i
            
        for k_i in range(dim):
            if k_i > 0:
                A[k_i-1,k_i] = 0.1
                A[k_i,k_i-1] = 0.1
        return A
    
    
    def _AOA(self, dim):
        """ Multiplication by matrix A from left and right
        
        """
        A = self._def_A(dim)
        
        for i_i in range(dim):
            for i_j in range(dim):
                for i_k in range(dim):
                    for i_l in range(dim):
                        self.data[i_i, i_j, i_k, i_l] = \
                        A[i_i, i_k]*A[i_l, i_j]
                        
                            
        
        