"""Liouvillian super operator

Louvillean super operator represent the action of a commutator with
Hamiltonian

Class Details
-------------

"""

from __future__ import annotations

from typing import Any

import numpy

from .superoperator import SuperOperator


class Liouvillian(SuperOperator):
    """Class representing and creating Liouvillian super operator


    Parameters
    ----------
    ham : Hamiltonian
        Hamiltonian from which we create the Liouvillian


    Examples
    --------
    >>> import numpy
    >>> import quantarhei as qr
    >>> HH = qr.qm.Hamiltonian(data=[[0.0, 1.0], [1.0, 0.2]])
    >>> rho = qr.qm.DensityMatrix(data=[[0.9, 0.05], [0.05, 0.1]])
    >>> Li = Liouvillian(HH)
    >>> com1 = numpy.dot(HH.data, rho.data) - numpy.dot(rho.data, HH.data)
    >>> com2 = Li.apply(rho)
    >>> numpy.allclose(com1, com2.data)
    True

    """

    def __init__(self, ham: Any) -> None:
        super().__init__(dim=ham.dim)
        dim = ham.dim
        H = ham.data
        I = numpy.eye(dim)
        self.data = (numpy.kron(H, I) - numpy.kron(I, H.T)).reshape(dim, dim, dim, dim)
