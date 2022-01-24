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

from ...core.managers import BasisManaged
from ...utils.types import BasisManagedComplexArray

class SuperOperator(BasisManaged):
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
    
    But calling uninitialize data attribute raises an exception
    >>> print(So.data is None)
    Traceback (most recent call last):
        ...
    AttributeError: 'SuperOperator' object has no attribute '_data'
    
    
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
    
    data = BasisManagedComplexArray("data")
    
    def __init__(self, dim=None, data=None, real=False):
        
        # Set the currently used basis
        cb = self.manager.get_current_basis()
        self.set_current_basis(cb)
        # unless it is the basis outside any context
        if cb != 0:
            self.manager.register_with_basis(cb,self)
            
        self._data_initialized = False
        
        self.dim = dim
        if dim is not None:
            if real:
                self.data = numpy.zeros((dim, dim, dim, dim),
                                        dtype=qr.REAL)
            else:
                self.data = numpy.zeros((dim, dim, dim, dim), 
                                        dtype=qr.COMPLEX)
        elif data is not None:
            self.data = data
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


    def transform(self, SS, inv=None):
        """Transforms the superoperator to a new basis
        

        Transformation of the superoperator follows similar rules
        as the one of operators. An operar :math:`A` is transformed as 
        follows:
        
        .. math::
            
            A_{\\alpha\\beta} = \\sum_{ij}(s^{-1})_{\\alpha i}A_{ij}s_{j\\beta}
            
        and 
        
        .. math::
            
            A_{ij} = \\sum_{\\alpha\\beta} s_{i\\alpha} \\
            A_{\\alpha\\beta}(s^{-1})_{\\beta j}
            
        The transformation matrix :math:`s` is obtained by standard
        diagonalization routines and its elements can be expressed in Dirac
        notation as
        
        .. math:: 
            
            s_{i\\alpha} = \\langle i | \\alpha \\rangle
            
        The inverse matrix :math:`s^{-1}` is obtained by transposition, because
        
        .. math::
            
            \\langle i | j \\rangle = \\delta_{ij} = \\
            \sum_{\\alpha} \\langle i | \\alpha \\rangle \\
            \\langle \\alpha | j \\rangle 
            
        and we see that 
        
        .. math::
            
            (s^{-1})_{\\alpha i} = s_{i \\alpha}.
            
        
        Given an operator :math:`A` which is
        a result of an action of the superoperator :math:`R` on operar
        :math:`A` we can see that the transformation occurs as follows:
        
        .. math::
            
            A_{ij} = \\sum_{kl}R_{ijkl}B_{kl}
            
        .. math::
            
            A_{\\alpha\\beta} = \\
            \\sum_{ij}(s^{-1})_{\\alpha i}A_{ij}s_{j\\beta} = \\
            \\sum_{ijkl} (s^{-1})_{\\alpha i} s_{j\\beta} R_{ijkl} B_{kl}
            
        Using the back transformation of the operator :math:`B` we obtaine
        
        .. math::
            
            A_{\\alpha\\beta} = \\
            \\sum_{ij}(s^{-1})_{\\alpha i}A_{ij}s_{j\\beta} = \\
            \\sum_{\\gamma\\delta} \\left [  \\right ] B_{\\gamma\\delta}

        which translates into
        
        .. math::
            
            R_{\\alpha\\beta\\gamma\\delta} = \\
            \\sum_{ijkl} s_{i \\alpha} \\
            s_{j\\beta} R_{ijkl} \\
            s_{k\\gamma}s_{l\\delta}

        
        Parameters
        ----------
        
        SS : float matrix
            Transformation matrix
            
        inv : float matrix, optional
            Inverse of the transformation matrix
    
        Examples
        --------
        
        """
        if (self.manager.warn_about_basis_change):
            print("\nQr >>> SuperOperator "+
                  "'%s' changes basis" %self.name)
        
        #
        # if inverse matrix not present, we create it
        #
        if inv is None:
            S1 = numpy.linalg.inv(SS)
        else:
            S1 = inv
            
        # dimension of the transformation matrix
        dim = SS.shape[0]
        
        #
        # Dimension 4 means a single, time independent superoperator 
        #
        if self._data.ndim == 4:
            for c in range(dim):
                for d in range(dim):
                    self._data[:,:,c,d] = \
                    numpy.dot(S1,numpy.dot(self._data[:,:,c,d],SS))
                    
            for a in range(dim):
                for b in range(dim):
                    self._data[a,b,:,:] = \
                    numpy.dot(S1,numpy.dot(self._data[a,b,:,:],SS))
                    
        #
        # Larger dimension means more superoperators or time dependence
        #
        else:

            for tt in range(self._data.shape[0]):
                for c in range(dim):
                    for d in range(dim):
                        self._data[tt,:,:,c,d] = \
                            numpy.dot(S1,numpy.dot(self._data[tt,:,:,c,d],SS))
                    
                for a in range(dim):
                    for b in range(dim):
                        self._data[tt,a,b,:,:] = \
                            numpy.dot(S1,numpy.dot(self._data[tt,a,b,:,:],SS))            
    