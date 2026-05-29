"""Quantarhei package (http://www.github.com/quantarhei)

statevector module


"""

from __future__ import annotations

from typing import Any

import numpy

from ... import COMPLEX
from ...core.managers import BasisManaged
from ...exceptions import QuantarheiError
from ...utils.types import BasisManagedComplexArray


class StateVector(BasisManaged):
    """Represents a quantum mechanical state vector


    Parameters
    ----------


    Examples
    --------
    >>> psi = StateVector(2)
    >>> print(psi.dim)
    2


    >>> vec = numpy.zeros((1,3), dtype=float)
    >>> psi = StateVector(data=vec)
    Traceback (most recent call last):
    ...
    quantarhei.exceptions.QuantarheiError: Data has to be a vector



    """

    data = BasisManagedComplexArray("data")
    _data: numpy.ndarray  # type: ignore[explicit-any]

    def __init__(  # type: ignore[explicit-any]
        self, dim: int | None = None, data: numpy.ndarray | None = None
    ) -> None:

        self._initialized = False

        # check and save data
        if data is not None:
            #            # list is accepted
            #            if isinstance(data, list):
            #                data = numpy.array(data)
            self.data = data

            if len(self.data.shape) != 1:
                raise QuantarheiError("Data has to be a vector")
            else:
                if dim is not None:
                    if self.data.shape[0] != dim:
                        raise QuantarheiError("Incompatible dim and data paramters")
                else:
                    self.dim = self.data.shape[0]
                self._initialized = True

        else:
            # check and save dim
            if dim is not None:
                self.dim = dim
                self.data = numpy.zeros(dim, dtype=COMPLEX)
                self._initialized = True

    def dot(self, vec: StateVector) -> Any:  # type: ignore[explicit-any]
        """Scalar product of two StateVectors"""
        return numpy.dot(self.data, vec.data)

    def norm(self) -> Any:  # type: ignore[explicit-any]
        """Returns the norm of the StateVector"""
        # vdot conjugates its first argument: <data|data> = ||data||^2
        return numpy.sqrt(numpy.vdot(self.data, self.data).real)

    def transform(self, SS: numpy.ndarray, inv: numpy.ndarray | None = None) -> None:  # type: ignore[explicit-any]
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

        self._data = numpy.dot(S1, self._data)

    def get_DensityMatrix(self) -> Any:  # type: ignore[explicit-any]
        """Constructs DensityMatrix from the present StateVector"""
        from .operators import DensityMatrix

        rho = DensityMatrix(dim=self.dim)

        for ii in range(self.dim):
            for jj in range(self.dim):
                rho.data[ii, jj] = self.data[ii] * numpy.conj(self.data[jj])

        return rho

    def __str__(self) -> str:
        out = "\nquantarhei.StateVector object"
        out += "\n============================="
        out += "\ndata = \n"
        out += str(self.data)
        return out
