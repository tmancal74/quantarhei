from __future__ import annotations

from typing import Any

import numpy

from ... import REAL
from ...exceptions import ImplementationError
from ..propagators.dmevolution import ReducedDensityMatrixEvolution


class OQSStateVectorEvolution:
    """Time evolution of a state vector under an open quantum system (OQS).

    Stores the state vector at each time step produced by an
    :class:`OQSStateVectorPropagator` and provides utilities to convert
    the result to a reduced density matrix evolution.

    Parameters
    ----------
    timeaxis : TimeAxis or None, optional
        Time grid for the propagation. Default is ``None``.
    psii : StateVector or None, optional
        Initial state vector. If provided, its dimension is used to
        allocate the storage array and it is set as the initial condition.
        Default is ``None``.
    """

    def __init__(self, timeaxis: Any = None, psii: Any = None) -> None:

        if timeaxis is not None:
            self.TimeAxis = timeaxis

            if psii is not None:
                self.dim = psii.data.shape[0]

                self.set_initial_condition(psii)
            else:
                self.dim = 0

    def set_initial_condition(self, psii: Any) -> None:
        """ """
        self.dim = psii.data.shape[0]
        self.data = numpy.zeros((self.TimeAxis.length, self.dim), dtype=numpy.float64)
        self.data[0, :] = psii.data

    def get_norm(self) -> numpy.ndarray:
        """Time dependent norm of the state vector"""
        return numpy.sum(self.data * self.data, axis=1)

    def get_ReducedDensityMatrixEvolution(
        self, decoherence: bool = False
    ) -> ReducedDensityMatrixEvolution:
        """Returns the corresponding reduced density matrix


        The present implementation cannot account for dephasing

        """
        if decoherence:
            raise ImplementationError("Decoherence not yet implemented.")

        rhoi = numpy.zeros((self.dim, self.dim), dtype=REAL)
        rhot = ReducedDensityMatrixEvolution(timeaxis=self.TimeAxis, rhoi=rhoi)
        for ii in range(self.dim):
            for jj in range(self.dim):
                rhot.data[:, ii, jj] = self.data[:, ii] * self.data[:, jj]

        return rhot
