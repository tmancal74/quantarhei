from __future__ import annotations

from typing import Any

import numpy

from ..exceptions import QuantarheiError


class Disorder:
    """Diagonal energetic disorder model for excitonic Hamiltonians.

    Draws random site-energy offsets from a specified distribution and
    applies them to the diagonal elements of a Hamiltonian on each
    disorder realisation.

    Parameters
    ----------
    data : numpy.ndarray
        Reference Hamiltonian data array whose diagonal entries are used as
        the unperturbed site energies.
    distribution : str, optional
        Statistical distribution for the disorder. Currently only
        ``"Gaussian"`` is supported. Default is ``"Gaussian"``.
    dtype : str, optional
        Type of disorder. Currently only ``"diagonal"`` (site-energy
        disorder) is supported. Default is ``"diagonal"``.
    seed : int or None, optional
        Seed for the NumPy random-number generator. Default is ``None``
        (random seed).

    Raises
    ------
    Exception
        If ``data`` is ``None``, or if an unknown distribution or disorder
        type is requested.
    """

    def __init__(  # type: ignore[explicit-any]
        self,
        data: numpy.ndarray | None = None,
        distribution: str = "Gaussian",
        dtype: str = "diagonal",
        seed: int | None = None,
    ) -> None:

        if data is None:
            raise QuantarheiError("Data not specified")

        self.dtype = dtype
        self.distribution = distribution
        self.seed_pool: list[Any] = []  # type: ignore[explicit-any]
        self.data = data.copy()
        self.shape = self.data.shape

        if seed is not None:
            numpy.random.seed(seed)

    def disorder_update(self, i_dis: int, H: Any, ignore_first: bool = False) -> None:  # type: ignore[explicit-any]
        """Adds disorder to an excitonic Hamiltonian"""
        if (i_dis == 0) and ignore_first:
            return

        N = H.data.shape[0] - 1

        if self.dtype == "diagonal":
            if self.distribution == "Gaussian":
                sigma = self.width / numpy.sqrt(2.0 * numpy.log(2))

                de = numpy.random.normal(0.0, sigma, N)

            else:
                raise QuantarheiError("Unknown distribution")

            # Update the Hamiltonian energies
            for i in range(N):
                H._data[1 + i, 1 + i] = self.data[1 + i, 1 + i] + de[i]

        else:
            raise QuantarheiError("Unknown disorder type")

    def set_distribution(self, distribution: str, params: dict) -> None:

        if distribution == "Gaussian":
            self.width = params["width"]

        else:
            raise QuantarheiError("Unknown distribution")
