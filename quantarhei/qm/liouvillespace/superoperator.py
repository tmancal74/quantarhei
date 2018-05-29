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
        
        
    Examples
    --------
    
    Creation of an empty SuperOperator is allowed
    
    >>> So = SuperOperator()
    >>> print(So.dim is None)
    True
    >>> print(So.data is None)
    True
    
    Creating with `dim` arguments creates a zero superoperator
    
    >>> So = SuperOperator(dim=3)
    >>> print(So.dim)
    3
    >>> print(So.data.shape)
    (3, 3, 3, 3)
    
    Creating superoperator from data checks, if the data object has the
    correct shape
    
    >>> data = numpy.zeros((3,3,3))
    >>> So = SuperOperator(data=data)
    Traceback (most recent call last):
        ...
    Exception: The data do not represent a superoperator
    
    >>> data = numpy.zeros((3,1,3,2))
    >>> So = SuperOperator(data=data)
    Traceback (most recent call last):
        ...
    Exception: `data` has to be `square` four-dimensional matrix
    
    """
    
    def __init__(self, dim=None, data=None, real=False):
        
        self.dim = dim
        self.data = data
        if dim is not None:
            if real:
                self.data = numpy.zeros((dim, dim, dim, dim),
                                        dtype=qr.REAL)
            else:
                self.data = numpy.zeros((dim, dim, dim, dim), 
                                        dtype=qr.COMPLEX)
        elif data is not None:
            if len(data.shape) != 4:
                raise Exception("The data do not represent a superoperator")
            Nd = data.shape[0]
            if numpy.any(numpy.array(data.shape) - Nd):
                raise Exception("`data` has to be `square` "+
                                "four-dimensional matrix")
            self.dim = data.shape[0]
      

    def apply(self, oper, copy=True):
        """Applies superoperator to an operator
        
        
        Parameters
        ----------
        
        oper : Operator
            Operator on which the present superoperator is applied
        
        
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> import numpy
        >>> op = qr.ReducedDensityMatrix(data=[[0.0, 1.0], [1.0, 0.0]])

        Let's take a matrix A
        
        >>> A = numpy.zeros((2,2))
        >>> A[0,1] = 1.0
        >>> A[1,0] = 1.0

        
        and create a superoperator equivalent to matrix multiplication
        
        .. math::
            
            A\\rho = \\sum_{i}A_{ik}\\rho_{kj} = \\
            \\sum_{kl}\\delta_{jl}A_{ik}\\rho_{kl} = \\
            \\sum_{kl}D_{ijkl}\\rho_{kl}

            D_{ijkl} = \\delta_{jl}A_{ik}
            
        
        >>> data = numpy.zeros((2,2,2,2))
        >>> for i in range(2):
        ...     for j in range(2):
        ...         for k in range(2):
        ...             for l in range(2):
        ...                 if j == l:
        ...                     data[i,j,k,l] = A[i,k]
        ...                 else:
        ...                     data[i,j,k,l] = 0.0
        >>> Dd = SuperOperator(data=data)
        
        Now we can check that maultiplication and application of the
        superoperator lead to the same result.
        
        >>> B1 = Dd.apply(op)
        >>> B2 = numpy.dot(A,op.data)
        >>> numpy.allclose(B1.data, B2)
        True
        
        Every time you do this, new density matrix object is created. To avoid
        this, we can apply the operation without copying the object.
        
        >>> B3 = Dd.apply(op, copy=False)
        >>> print(B3 == op)
        True
        
        `B3` is now a different reference to the same object as `op`. We can
        check that the result of the operation is still correct:
        
        >>> numpy.allclose(B3.data, B2)
        True
        """
        
        if copy:
            import copy
            oper_ven = copy.copy(oper)
            oper_ven.data = numpy.tensordot(self.data, oper.data)
            return oper_ven
        else:
            oper.data = numpy.tensordot(self.data, oper.data)
            return oper
    