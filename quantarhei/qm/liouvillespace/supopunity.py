# -*- coding: utf-8 -*-
""" Class representing unity superoperator


    Class Details
    -------------

"""

import numpy

from .superoperator import SuperOperator


class SOpUnity(SuperOperator):
    """Class representing a unity superoperator
    
    
    Parameters
    ----------
    
    dim : int
        Dimension of the unity superoperator
        
    data : array
        If data is specified, only their first dimension is used to construct
        a unity superoperator
        
        
    Examples
    --------
    
    Creation of an empty operator is not allowed
    
    >>> empty = SOpUnity()
    Traceback (most recent call last):
        ...
    Exception: Dimension of the superoperator has to be defined
    
    Creating unity superoperator of defined dimensions
    
    >>> SU = SOpUnity(dim=3)
    >>> print(SU.data.shape)
    (3, 3, 3, 3)
    
    >>> import numpy
    >>> A = numpy.array([[0,    1.0, 2.0],
    ...                  [1.0, -3.0, 0.0],
    ...                  [-2.0, 0.0, 1.0]])
    
    Application of SU on A should not change the matrix, because it
    corresponds to multiplication by unity
    
    >>> B = SU.apply(A)
    >>> numpy.allclose(A,B)
    True
    
    Knowing the dimension of A, we can create corresponding unity superoperator
    using A as the `data` argument:
        
    >>> SU2 = SOpUnity(data=A)
    >>> C = SU2.apply(A)
    >>> numpy.allclose(A,C)
    True
    
    """
    
    def __init__(self, dim=None, data=None):

        if data is not None:
            dim = data.shape[0]
            
        super().__init__(dim=dim, data=None)
            
        # initialize the data
        if dim is not None:
            
            import quantarhei as qr
            self.data = numpy.zeros((dim, dim, dim, dim), qr.REAL)
            for i in range(self.dim):
                for j in range(self.dim):
                    self.data[i,j,i,j] = 1.0
                
        else:
            
            raise Exception("Dimension of the superoperator has to be defined")
                    
                
            
        
