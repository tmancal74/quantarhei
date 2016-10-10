# -*- coding: utf-8 -*-
"""
    Quantarhei package (http://www.github.com/quantarhei)

    statevector module


"""

import numpy

class StateVector:
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
    def __init__(self, dim=None, data=None):

        self._initialized = False

        # check and save data
        if data is not None:

            # list is accepted
            if isinstance(data, list):
                data = numpy.array(data)

            shape = data.shape
            if len(shape) != 1:
                raise Exception("Data has to be a vector")
            else:
                self.data = data
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


    def __str__(self):
        out  = "\nquantarhei.StateVector object"
        out += "\n============================="
        out += "\ndata = \n"
        out += str(self.data)
        return out