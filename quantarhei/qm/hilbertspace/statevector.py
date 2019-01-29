# -*- coding: utf-8 -*-
"""
    Quantarhei package (http://www.github.com/quantarhei)

    statevector module


"""

import numpy

from ...utils.types import BasisManagedComplexArray
from ...core.managers import BasisManaged

class StateVector(BasisManaged):
    """Represents a quantum mechanical state vector


    Parameters
    ----------


    Examples
    --------

    >>> psi = StateVector(2)
    >>> print(psi.dim)
    2


    >>> vec = numpy.zeros((1,3), dtype=numpy.float)
    >>> psi = StateVector(data=vec)
    Traceback (most recent call last):
    ...
    Exception: Data has to be a vector



    """

    data = BasisManagedComplexArray("data")   

    def __init__(self, dim=None, data=None):

        self._initialized = False

        # check and save data
        if data is not None:

#            # list is accepted
#            if isinstance(data, list):
#                data = numpy.array(data)
            self.data = data

            if len(self.data.shape) != 1:
                raise Exception("Data has to be a vector")
            else:
                if dim is not None:
                    if self.data.shape[0] != dim:
                        raise Exception("Incompatible dim and data paramters")
                else:
                    self.dim = self.data.shape[0]
                self._initialized = True

        else:

            # check and save dim
            if dim is not None:
                self.dim = dim
                self.data = numpy.zeros(dim, dtype=numpy.float)
                self._initialized = True



    def dot(self, vec):
        """Scalar product of two StateVectors

        """

        return numpy.dot(self.data, vec.data)

    def norm(self):
        """Returns the norm of the StateVector

        """
        return numpy.sqrt(numpy.dot(self.data, self.data))


    def transform(self, SS, inv=None):
        """Transformation of the operator by a given matrix
        
        
        This function transforms the Operator into a different basis, using
        a given transformation matrix.
        
        Parameters
        ----------
        
        SS : matrix, numpy.ndarray
            transformation matrix
            
        inv : matrix, numpy.ndarray
            inverse of the transformation matrix
            
        """        
        if inv is None:
            S1 = numpy.linalg.inv(SS)
        else:
            S1 = inv

        self._data = numpy.dot(S1,self._data)
        

    def get_DensityMatrix(self):
        """Constructs DensityMatrix from the present StateVector
        
        """
        from .operators import DensityMatrix

        rho = DensityMatrix(dim=self.dim)
        
        for ii in range(self.dim):
            for jj in range(self.dim):
                rho.data[ii,jj] = self.data[ii]*numpy.conj(self.data[jj])
        
        return rho
        
        
    def __str__(self):
        out  = "\nquantarhei.StateVector object"
        out += "\n============================="
        out += "\ndata = \n"
        out += str(self.data)
        return out