from __future__ import annotations

from collections.abc import Callable

import numpy
import scipy.interpolate as interp

from ....core.time import TimeDependent
from ...corfunctions.correlationfunctions import c2g
from ...hilbertspace.hamiltonian import Hamiltonian
from ...liouvillespace.systembathinteraction import SystemBathInteraction


class TDFoersterRateMatrix(TimeDependent):
    """Time-dependent Förster population relaxation rate matrix.

    The returned data array has shape ``(Nt, N, N)`` and follows the
    convention ``data[:, a, b]`` is the population transfer rate from state
    ``b`` to state ``a``. Diagonal elements are depopulation rates.
    """

    def __init__(
        self,
        ham: Hamiltonian,
        sbi: SystemBathInteraction,
        initialize: bool = True,
        cutoff_time: float | None = None,
    ) -> None:

        if not isinstance(ham, Hamiltonian):
            raise Exception("First argument must be a Hamiltonian")

        if not isinstance(sbi, SystemBathInteraction):
            raise Exception("Second argument must be a SystemBathInteraction")

        self._is_initialized = False
        self._has_cutoff_time = False
        self.is_time_dependent = True

        if cutoff_time is not None:
            self.cutoff_time = cutoff_time
            self._has_cutoff_time = True

        self.ham = ham
        self.sbi = sbi

        if initialize:
            self.initialize()
            self._is_initialized = True

    def initialize(self) -> None:
        """Calculates time-dependent Förster population rates."""
        time = self.sbi.TimeAxis
        assert time is not None, "SystemBathInteraction must have a TimeAxis set"

        tt = time.data
        Nt = time.length
        HH = self.ham.data
        Na = self.ham.dim

        cc = self.sbi.CC
        assert cc is not None, "SystemBathInteraction must have CC set"

        gt = numpy.zeros((Na, Nt), dtype=numpy.complex64)
        for ii in range(1, Na):
            gt[ii, :] = c2g(time, cc.get_coft(ii - 1, ii - 1))

        ll = numpy.zeros(Na)
        for ii in range(1, Na):
            ll[ii] = cc.get_reorganization_energy(ii - 1, ii - 1)

        self.data = _td_reference_implementation(Na, Nt, HH, tt, gt, ll, _td_fintegral)

        for bb in range(Na):
            self.data[:, bb, bb] = -numpy.sum(self.data[:, :, bb], axis=1)

        self._is_initialized = True


def _td_reference_implementation(  # type: ignore[explicit-any]
    Na: int,
    Nt: int,
    HH: numpy.ndarray,
    tt: numpy.ndarray,
    gt: numpy.ndarray,
    ll: numpy.ndarray,
    fce: Callable[..., numpy.ndarray],
) -> numpy.ndarray:

    #
    # Rates between states a and b
    #
    KK = numpy.zeros((Nt, Na, Na), dtype=numpy.float64)
    for a in range(Na):
        for b in range(Na):
            if a != b:
                ed = HH[b, b]  # donor
                ea = HH[a, a]  # acceptor
                KK[:, a, b] = (HH[a, b] ** 2) * fce(
                    tt, gt[a, :], gt[b, :], ed, ea, ll[b]
                )

    return KK


def _td_fintegral(  # type: ignore[explicit-any]
    tt: numpy.ndarray,
    gtd: numpy.ndarray,
    gta: numpy.ndarray,
    ed: float,
    ea: float,
    ld: float,
) -> numpy.ndarray:
    """Time dependent Foerster integral


    Parameters
    ----------
    tt : numpy array
        Time

    gtd : numpy array
        lineshape function of the donor transition

    gta : numpy array
        lineshape function of the acceptor transition

    ed : float
        Energy of the donor transition

    ea : float
        Energy of the acceptor transition

    ld : float
        Reorganization energy of the donor

    Returns
    -------
    ret : float
        The value of the Foerster integral

    """
    prod = numpy.exp(-gtd - gta + 1j * ((ed - ea) - 2.0 * ld) * tt)

    preal = numpy.real(prod)
    pimag = numpy.imag(prod)
    splr = interp.UnivariateSpline(tt, preal, s=0).antiderivative()(tt)
    spli = interp.UnivariateSpline(tt, pimag, s=0).antiderivative()(tt)
    hoft = splr + 1j * spli

    ret = 2.0 * numpy.real(hoft)

    return ret
