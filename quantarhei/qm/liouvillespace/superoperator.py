# -*- coding: utf-8 -*-
""" Class representing superoperators
    
    
    This class represents operators on the space of Hilbert space operators.
    Usually, we refer to such operators as superoperators. 




    Class Details
    -------------

"""

# dependencies imports
import numpy

# quantarhei imports
import quantarhei as qr

# FIXME: This class should be a base class for Relaxation tensors
class SuperOperator:
    """Class representing superoperators
    
    
    Parameters
    ----------
    
    dim : int
        Dimension of the superoperator
        
    data : array
        Data of the superoperator
        
    real : bool
        Is this data real? False if they are complex
    
    """
    
    def __init__(self, dim=None, data=None, real=False):
        
        if dim is not None:
            self.dim = dim
            if real:
                self.data = numpy.zeros((dim, dim, dim, dim),
                                        dtype=qr.REAL)
            else:
                self.data = numpy.zeros((dim, dim, dim, dim), 
                                        dtype=qr.COMPLEX)
        elif data is not None:
            self.data = data
            self.dim = data.shape[0]
      

    def apply(self, oper):
        """Applies superoperator to an operator
        
        
        Parameters
        ----------
        
        oper : Operator
            Operator on which the present superoperator is applied
        
        
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> op = qr.ReducedDensityMatrix(data=[[0.0, 1.0], [1.0, 0.0]])
        >>> print(op)
        <BLANKLINE>
        quantarhei.ReducedDensityMatrix object
        ======================================
        data = 
        [[ 0.  1.]
         [ 1.  0.]]
        
        """
        
        oper.data = numpy.tensordot(self.data, oper.data)
        
        return oper
    