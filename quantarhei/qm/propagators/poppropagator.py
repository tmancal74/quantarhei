from __future__ import annotations

from typing import Any

import numpy

from ...core.time import TimeAxis
from ...exceptions import QuantarheiError
from ..liouvillespace.rates.ratematrix import RateMatrix


class PopulationPropagator:
    """Propagator for a population vector

    The most important method of this class is the method `propagate`.
    It takes initial population vector and propagates it using the time values
    specified by TimaAxis.

    Parameters
    ----------
    timeaxis: TimeAxis
        The time interval on which we propagate.

    RateMatrix: float array of rank 2
        Rate matrix used to propagate the populations


    Examples
    --------

    """

    def __init__(self, timeaxis: TimeAxis, rate_matrix: Any = None) -> None:  # type: ignore[explicit-any]

        self.timeAxis = timeaxis
        self.Nref = 1
        self.Nt = self.timeAxis.length
        self.dt = self.timeAxis.step

        if rate_matrix is not None:
            if isinstance(rate_matrix, RateMatrix) or (
                hasattr(rate_matrix, "data")
                and not isinstance(rate_matrix, numpy.ndarray)
            ):
                self.KK = rate_matrix.data
            else:
                self.KK = rate_matrix
            self.KK = numpy.asarray(self.KK)

    def propagate(self, pini: Any) -> numpy.ndarray:  # type: ignore[explicit-any]
        """Propagates a given initional population vector"""
        if not isinstance(pini, numpy.ndarray):
            pini = numpy.array(pini)

        return self._propagate_short_exp(pini)

    def _propagate_short_exp(self, pini: numpy.ndarray, L: int = 4) -> numpy.ndarray:  # type: ignore[explicit-any]
        """Propagation of the initial pop vector by short time expansion

        Propagates initial population vector by a give rate matrix with
        constant coefficients


        """
        Nt = self.timeAxis.length
        pops = numpy.zeros((Nt, pini.shape[0]))

        pops[0, :] = pini

        rho1 = pini
        rho2 = pini

        indx = 1
        for ii in self.timeAxis.data[1 : self.Nt]:
            for jj in range(0, self.Nref):
                for ll in range(1, L + 1):
                    pref = self.dt / ll
                    rho1 = pref * numpy.dot(self.KK, rho1)

                    rho2 = rho2 + rho1
                rho1 = rho2

            pops[indx, :] = rho2
            indx += 1

        return pops

    def _split_relaxation_matrix(self) -> tuple:
        """Splits the relaxation matrix into diagonal and transfer parts"""
        N = self.KK.shape[0]
        KKD = numpy.zeros(N, dtype=numpy.float64)
        KKT = numpy.zeros((N, N), dtype=numpy.float64)
        for i in range(N):
            KKD[i] = -self.KK[i, i]
        KKT = self.KK + numpy.diag(KKD)
        return KKD, KKT

    def _integrateK0(  # type: ignore[explicit-any]
        self, timeaxis: TimeAxis, Kab: float, Ka: float, Kb: float
    ) -> numpy.ndarray:
        """Returns an integral of transfer matrix element"""
        expK = numpy.exp((Ka - Kb) * timeaxis.data)
        integ = numpy.zeros(timeaxis.length, dtype=numpy.float64)
        for i in range(timeaxis.length):
            integ[i] = numpy.sum(expK[0:i])
        return Kab * integ * timeaxis.step

    def _integrateKn(  # type: ignore[explicit-any]
        self, timeaxis: TimeAxis, Kab: float, Ka: float, Kb: float, Kbc: numpy.ndarray
    ) -> numpy.ndarray:
        """Returns an integral of transfer matrix element multiplied by
        a result of previous order of expansion


        """
        expK = numpy.exp((Ka - Kb) * timeaxis.data)
        integ = numpy.zeros(timeaxis.length, dtype=numpy.float64)
        for i in range(timeaxis.length):
            integ[i] = numpy.sum(expK[0:i] * Kbc[0:i])
        return Kab * integ * timeaxis.step

    def _time_indices(self, timeaxis: TimeAxis) -> numpy.ndarray:  # type: ignore[explicit-any]
        """Returns indices of ``timeaxis`` points on the internal time axis"""
        if not timeaxis.is_subset_of(self.timeAxis):
            raise Exception(
                "TimeAxis is not a subset of the internal TimeAxis of this propagator."
            )

        indices = numpy.rint(
            (timeaxis.data - self.timeAxis.start) / self.timeAxis.step
        ).astype(numpy.int64)
        times = self.timeAxis.start + indices * self.timeAxis.step
        if not numpy.allclose(times, timeaxis.data):
            raise Exception(
                "TimeAxis is not a subset of the internal TimeAxis of this propagator."
            )

        return indices

    def _rate_matrix_time_series(self) -> numpy.ndarray:  # type: ignore[explicit-any]
        """Returns rate matrices as an array with shape ``(Nt, N, N)``"""
        KK = numpy.asarray(self.KK)

        if len(KK.shape) == 2:
            if KK.shape[0] != KK.shape[1]:
                raise Exception("Rate matrix has to be square")
            rates = numpy.zeros(
                (self.timeAxis.length, KK.shape[0], KK.shape[1]), dtype=numpy.float64
            )
            rates[:, :, :] = KK
            return rates

        if len(KK.shape) == 3:
            if KK.shape[0] != self.timeAxis.length:
                raise Exception(
                    "Time-dependent rate matrix has to have the same length "
                    "as the internal TimeAxis"
                )
            if KK.shape[1] != KK.shape[2]:
                raise Exception("Rate matrix has to be square")
            return KK.astype(numpy.float64, copy=False)

        raise Exception("Rate matrix has to be an array of rank 2 or 3")

    def _to_matrix_time_order(self, matrix_time: numpy.ndarray) -> numpy.ndarray:  # type: ignore[explicit-any]
        """Converts ``(Nt, N, N)`` arrays to the historical ``(N, N, Nt)`` order"""
        return numpy.transpose(matrix_time, (1, 2, 0))

    def _integrate_matrix_time_series(  # type: ignore[explicit-any]
        self, matrix_time: numpy.ndarray
    ) -> numpy.ndarray:
        """Integrates a matrix-valued time series by the trapezoidal rule"""
        Nt = self.timeAxis.length
        integral = numpy.zeros_like(matrix_time)
        cumulative = numpy.zeros(matrix_time.shape[1:], dtype=matrix_time.dtype)

        for ii in range(Nt - 1):
            dt = self.timeAxis.data[ii + 1] - self.timeAxis.data[ii]
            cumulative += 0.5 * (matrix_time[ii] + matrix_time[ii + 1]) * dt
            integral[ii + 1] = cumulative.copy()

        return integral

    def _no_jump_propagator(self, rates: numpy.ndarray) -> numpy.ndarray:  # type: ignore[explicit-any]
        """Calculates the diagonal no-jump propagator U0(t)"""
        Nt = self.timeAxis.length
        N = rates.shape[1]
        U0 = numpy.zeros((Nt, N, N), dtype=numpy.float64)
        U0[0] = numpy.eye(N)

        cumulative = numpy.zeros(N, dtype=numpy.float64)
        for ii in range(Nt - 1):
            dt = self.timeAxis.data[ii + 1] - self.timeAxis.data[ii]
            diag_now = numpy.diag(rates[ii])
            diag_next = numpy.diag(rates[ii + 1])
            cumulative += 0.5 * (diag_now + diag_next) * dt
            for jj in range(N):
                U0[ii + 1, jj, jj] = numpy.exp(cumulative[jj])

        return U0

    def _interaction_picture_rates(  # type: ignore[explicit-any]
        self, rates: numpy.ndarray
    ) -> tuple[numpy.ndarray, numpy.ndarray]:
        """Returns U0(t) and interaction-picture transfer rates"""
        Nt = self.timeAxis.length
        N = rates.shape[1]
        U0 = self._no_jump_propagator(rates)

        K1 = rates.copy()
        for ii in range(Nt):
            numpy.fill_diagonal(K1[ii], 0.0)

        U0_inv = numpy.zeros_like(U0)
        for ii in range(Nt):
            for jj in range(N):
                val = U0[ii, jj, jj]
                U0_inv[ii, jj, jj] = 1.0 / val if val != 0.0 else 0.0

        KI = numpy.zeros_like(rates)
        for ii in range(Nt):
            KI[ii] = numpy.dot(U0_inv[ii], numpy.dot(K1[ii], U0[ii]))

        return U0, KI

    def _propagation_matrix_time_dependent(self, timeaxis: TimeAxis) -> numpy.ndarray:  # type: ignore[explicit-any]
        """Returns propagation matrix for time-dependent rates by RK4."""
        indices = self._time_indices(timeaxis)
        rates = self._rate_matrix_time_series()
        Nt = self.timeAxis.length
        N = rates.shape[1]

        U = numpy.zeros((Nt, N, N), dtype=numpy.float64)
        U[0] = numpy.eye(N)

        for ii in range(Nt - 1):
            dt = self.timeAxis.data[ii + 1] - self.timeAxis.data[ii]
            K1 = rates[ii]
            K4 = rates[ii + 1]
            K23 = 0.5 * (K1 + K4)

            k1 = numpy.dot(K1, U[ii])
            k2 = numpy.dot(K23, U[ii] + 0.5 * dt * k1)
            k3 = numpy.dot(K23, U[ii] + 0.5 * dt * k2)
            k4 = numpy.dot(K4, U[ii] + dt * k3)

            U[ii + 1] = U[ii] + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0

        return self._to_matrix_time_order(U[indices])

    def get_JumpExpansion(self, timeaxis: TimeAxis = None, max_order: int = 3) -> tuple:  # type: ignore[explicit-any]
        """Returns jump-expansion terms of the population propagator.

        The returned tuple contains matrices ``(U0, U1, ..., Un)`` in the
        historical Quantarhei order ``(N, N, Nt)``. ``U0`` is the no-jump
        diagonal evolution, and ``Un`` contains paths with ``n`` population
        transfers. The sum of all orders approximates the full propagator.
        """
        if max_order < 0:
            raise ValueError("max_order must be non-negative")

        if timeaxis is None:
            timeaxis = self.timeAxis

        indices = self._time_indices(timeaxis)
        rates = self._rate_matrix_time_series()
        U0, KI = self._interaction_picture_rates(rates)

        propagators = [self._to_matrix_time_order(U0[indices])]
        if max_order == 0:
            return tuple(propagators)

        Nt = self.timeAxis.length
        N = rates.shape[1]
        previous_I = numpy.zeros_like(rates)
        previous_I[:] = numpy.eye(N)

        for _order in range(1, max_order + 1):
            integrand = numpy.zeros_like(rates)
            for ii in range(Nt):
                integrand[ii] = numpy.dot(KI[ii], previous_I[ii])

            current_I = self._integrate_matrix_time_series(integrand)
            current = numpy.zeros_like(rates)
            for ii in range(Nt):
                current[ii] = numpy.dot(U0[ii], current_I[ii])

            propagators.append(self._to_matrix_time_order(current[indices]))
            previous_I = current_I

        return tuple(propagators)

    def get_KineticDecomposition(  # type: ignore[explicit-any]
        self, timeaxis: TimeAxis = None, max_order: int = 3
    ) -> tuple:
        """Alias for :meth:`get_JumpExpansion`."""
        return self.get_JumpExpansion(timeaxis=timeaxis, max_order=max_order)

    def get_PropagationMatrix(  # type: ignore[explicit-any]
        self, timeaxis: TimeAxis, corrections: int = -1, exact: bool = False
    ) -> Any:
        """Returns propagation matrix corresponding to the present propagator

        This function also returns perturbative expansion of the propagation
        in the orders of transfer rates.


        Examples
        --------
        >>> from quantarhei import TimeAxis
        >>> ta = TimeAxis(0.0, 1000, 1.0)
        >>> ts = TimeAxis(0.0, 10, 10.0)
        >>> KK = [[-1.0/100.0,  1.0/100.0],
        ...       [ 1.0/100.0, -1.0/100.0]]
        >>> prop = PopulationPropagator(ta, rate_matrix=numpy.array(KK))
        >>> U = prop.get_PropagationMatrix(ts)
        >>> print(U[:,:,0])
        [[ 1.  0.]
         [ 0.  1.]]


        """
        if timeaxis.is_subset_of(self.timeAxis):
            if len(self.KK.shape) == 3:
                U = self._propagation_matrix_time_dependent(timeaxis)
                if corrections > -1:
                    return U, self.get_JumpExpansion(timeaxis, max_order=corrections)
                return U

            N = self.KK.shape[0]
            U = numpy.zeros((N, N, timeaxis.length), dtype=numpy.float64)

            # initial condition
            U0 = numpy.eye(N)

            # diagonalization of the rate matrix
            Kd, SS = numpy.linalg.eig(self.KK)
            S1 = numpy.linalg.inv(SS)

            # calculating exp(KK*step)
            expKd_step = numpy.dot(
                SS, numpy.dot(numpy.diag(numpy.exp(Kd * timeaxis.step)), S1)
            )

            #
            # If the starts of the time axes do not coincide, and the
            # shift is not a multiple of the submitted time axis step,
            # we have to calculate the evolution matrix also with
            # the difference of the starting times
            #
            if self.timeAxis.start != timeaxis.start:
                # get the distance of the starts in timeaxis.steps
                Ns = round((timeaxis.start - self.timeAxis.start) / timeaxis.step)
                # if one can use timeaxis.step, make Ns steps
                if timeaxis.start == self.timeAxis.start + Ns * timeaxis.step:
                    for i in range(Ns):
                        U0 = numpy.dot(expKd_step, U0)
                # otherwise new exp(KK*dt) has to be calculated and applied
                else:
                    dt = timeaxis.start - self.timeAxis.start
                    expKd_dt = numpy.dot(
                        SS, numpy.dot(numpy.diag(numpy.exp(Kd * dt)), S1)
                    )
                    U0 = numpy.dot(expKd_dt, U0)

            # initial condition at the start of the submitted timeaxis
            U[:, :, 0] = U0

            # application of the exp(KK*step) on initial condition
            for i in range(1, timeaxis.length):
                U[:, :, i] = numpy.dot(expKd_step, U[:, :, i - 1])

            #
            # Calculate exact orders in transfer matrix
            #
            if corrections > -1:
                #
                # Split relaxation matrix into diagonal and transfer parts
                #
                KKD, KKT = self._split_relaxation_matrix()

                #
                # zero's order correction to the evolution matrix
                # (exact version)
                #
                Uc0 = numpy.zeros((N, N, timeaxis.length), dtype=numpy.float64)
                for i in range(N):
                    Uc0[i, i, :] = numpy.exp(-KKD[i] * timeaxis.data)

                if corrections == 0:
                    return U, (Uc0,)

            #
            # The rest of the corrections (Uc1, Uc2 ...) is not used at the moment
            #

            if (corrections > 0) and exact:
                #
                # first order correction
                # (exact version)
                #
                Uc1 = numpy.zeros((N, N, timeaxis.length), dtype=numpy.float64)
                for i in range(N):
                    for j in range(N):
                        if i != j:
                            if KKD[i] == KKD[j]:
                                Uc1[i, j, :] = KKT[i, j] * timeaxis.data
                            else:
                                Uc1[i, j, :] = (KKT[i, j] / (KKD[i] - KKD[j])) * (
                                    numpy.exp(-KKD[j] * timeaxis.data)
                                    - numpy.exp(-KKD[i] * timeaxis.data)
                                )

                if corrections == 1:
                    return U, (Uc0, Uc1)

                #
                # second order correction
                # (exact version)
                #
                Uc2 = numpy.zeros((N, N, timeaxis.length), dtype=numpy.float64)
                for i in range(N):
                    for j in range(N):
                        for k in range(N):
                            if k != j:
                                if i != j:
                                    # two parts of an expressions for i != i
                                    Uc2[i, j, :] += (
                                        KKT[i, k] * KKT[k, j] / (KKD[k] - KKD[j])
                                    ) * (
                                        (
                                            numpy.exp(-KKD[j] * timeaxis.data)
                                            - numpy.exp(-KKD[i] * timeaxis.data)
                                        )
                                        / (KKD[i] - KKD[j])
                                    )
                                    if k != i:
                                        Uc2[i, j, :] -= (
                                            KKT[i, k] * KKT[k, j] / (KKD[k] - KKD[j])
                                        ) * (
                                            (
                                                numpy.exp(-KKD[k] * timeaxis.data)
                                                - numpy.exp(-KKD[i] * timeaxis.data)
                                            )
                                            / (KKD[i] - KKD[k])
                                        )
                                else:
                                    # whole expression for i = j
                                    Uc2[i, j, :] += (
                                        KKT[i, k] * KKT[k, j] / (KKD[k] - KKD[j])
                                    ) * (
                                        numpy.exp(-KKD[i] * timeaxis.data)
                                        * timeaxis.data
                                        - (
                                            numpy.exp(-KKD[k] * timeaxis.data)
                                            - numpy.exp(-KKD[i] * timeaxis.data)
                                        )
                                        / (KKD[i] - KKD[k])
                                    )
                                    pass

                if corrections >= 2:
                    return U, (Uc0, Uc1, Uc2)

            elif (corrections > 0) and not exact:
                #
                # first order correction
                # (numerical)
                #
                Uc1 = numpy.zeros((N, N, timeaxis.length), dtype=numpy.float64)
                KKI = numpy.zeros((N, N, timeaxis.length), dtype=numpy.float64)
                for i in range(N):
                    for j in range(N):
                        KKI[i, j, :] = self._integrateK0(
                            timeaxis, KKT[i, j], KKD[i], KKD[j]
                        )
                        Uc1[i, j, :] = numpy.exp(-KKD[i] * timeaxis.data) * KKI[i, j, :]
                if corrections == 1:
                    return U, (Uc0, Uc1)

                #
                # second and higher order corrections
                # (numerical)
                #
                higher_c = [Uc0, Uc1]
                for c in range(1, corrections):
                    Uc2 = numpy.zeros((N, N, timeaxis.length), dtype=numpy.float64)
                    for i in range(N):
                        for j in range(N):
                            for k in range(N):
                                Uc2[i, j, :] += self._integrateKn(
                                    timeaxis, KKT[i, k], KKD[i], KKD[k], KKI[k, j, :]
                                )

                    KKI[:, :, :] = Uc2[
                        :, :, :
                    ]  # ! we do not want to assign a pointer !!!
                    for i in range(N):
                        for j in range(N):
                            Uc2[i, j, :] = (
                                numpy.exp(-KKD[i] * timeaxis.data) * KKI[i, j, :]
                            )
                    higher_c.append(Uc2)

                if corrections >= 2:
                    return U, tuple(higher_c)

            else:
                return U

        else:
            raise QuantarheiError(
                "TimeAxis is not a subset of the internal TimeAxis of this propagator."
            )
