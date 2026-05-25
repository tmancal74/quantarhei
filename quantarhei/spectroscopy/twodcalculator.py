from __future__ import annotations

import os
from typing import Any, Literal

import numpy

from .. import signal_NONR, signal_REPH
from ..builders.aggregates import Aggregate
from ..builders.molecules import Molecule
from ..builders.opensystem import OpenSystem
from ..core.managers import Manager, energy_units
from ..core.time import TimeAxis
from ..implementations.aceto.lab_settings import lab_settings

# deprecated class
# This is how we calculate it now
from ..spectroscopy.responses import (
    LiouvillePathway,
    NonLinearResponse,
    get_common_time_axis,
    get_single_exciton_rate_matrix,
    validate_2d_time_axes,
)
from ..utils import derived_type
from .twodresponse import TwoDResponse


def _apply_response_window(data: numpy.ndarray) -> numpy.ndarray:
    """Apply the endpoint half-weight used before Fourier transformation."""
    ret = data.copy()
    ret[:, 0] *= 0.5
    ret[0, :] *= 0.5
    return ret


def _fourier_transform_response(data: numpy.ndarray, signal: str) -> numpy.ndarray:
    """Transform a time-domain response contribution to a 2D spectrum."""
    data = _apply_response_window(data)

    if signal == signal_REPH:
        ftresp = numpy.fft.fft(data, axis=1)
    elif signal == signal_NONR:
        ftresp = numpy.fft.ifft(data, axis=1) * data.shape[1]
    else:
        raise Exception("Unknown 2D signal type: " + signal)

    ftresp = numpy.fft.ifft(ftresp, axis=0)
    return numpy.fft.fftshift(ftresp)


def _pad_response_data(
    data: numpy.ndarray, pad: int, window: numpy.ndarray | None = None
) -> numpy.ndarray:
    """Pad a response contribution in the same way as the total response."""
    if window is not None:
        size = int(len(window) / 2)
        data = data.copy()
        data[len(data) - size :, :] *= window[size:, None]
        data[:, len(data) - size :] *= window[None, size:]

    if pad > 0:
        data = numpy.hstack((data, numpy.zeros((data.shape[0], pad))))
        data = numpy.vstack((data, numpy.zeros((pad, data.shape[1]))))

    return data


def _normalize_twodtype(twodtype: str) -> str:
    if twodtype in ["2DES", "F-2DES"]:
        return twodtype
    raise Exception("Unknown type of 2D spectrum: " + twodtype)


class TwoDResponseCalculator:
    """Calculator of the 2D spectrum


    Enables setting up parameters of 2D spectrum calculation for later
    evaluation. The method `calculate` returns TwoDResponseContainer
    with a 2D spectrum.

    Parameters
    ----------
    twodtype : {"2DES", "F-2DES"}
        Type of 2D spectrum to calculate. ``"F-2DES"`` currently supports
        the ``gamma_factor`` fluorescence-detected approximation.

    gamma_factor : float
        Fluorescence-detected 2D weighting factor. ESA contributions are
        weighted by ``gamma_factor - 1``.

    population_factors : tuple
        Reserved for state-resolved fluorescence-detected 2D spectra. This
        mode is not implemented yet because it requires state-resolved ESA
        responses.


    """

    t1axis = derived_type("t1axis", TimeAxis)
    t2axis = derived_type("t2axis", TimeAxis)
    t3axis = derived_type("t3axis", TimeAxis)

    system = derived_type("system", [Molecule, Aggregate, OpenSystem])

    _has_responses = False
    _has_system = False

    def __init__(
        self,
        t1axis: Any,
        t2axis: Any,
        t3axis: Any,
        system: Any = None,
        responses: Any = None,
        dynamics: str = "secular",
        relaxation_tensor: Any = None,
        rate_matrix: Any = None,
        relaxation_theory: str | None = None,
        rate_matrix_time_dependent: bool = False,
        relaxation_cutoff_time: float | None = None,
        rate_matrix_options: dict[str, Any] | None = None,
        effective_hamiltonian: Any = None,
        twodtype: str = "2DES",
        gamma_factor: float | None = None,
        population_factors: Any = None,
        jump_order: int = 0,
        jump_time_graining: int = 1,
        jump_kernel_cutoff: float = 0.0,
        jump_kernel_zero_cutoff: float = 0.0,
    ) -> None:

        twodtype = _normalize_twodtype(twodtype)
        if twodtype == "F-2DES":
            if gamma_factor is None and population_factors is None:
                raise Exception("Not enough parameters for F-2DES")
            if population_factors is not None:
                raise NotImplementedError(
                    "F-2DES with population_factors requires state-resolved "
                    "ESA responses and is not implemented yet."
                )

        if jump_order not in (0, 1):
            raise ValueError("2D response jump_order has to be 0 or 1")
        if jump_time_graining < 1:
            raise ValueError("jump_time_graining has to be a positive integer")
        if jump_kernel_cutoff < 0.0:
            raise ValueError("jump_kernel_cutoff has to be non-negative")
        if jump_kernel_zero_cutoff < 0.0:
            raise ValueError("jump_kernel_zero_cutoff has to be non-negative")

        self.t1axis = t1axis
        self.t2axis = t2axis
        self.t3axis = t3axis
        self.twodtype = twodtype
        self.gamma_factor = gamma_factor
        self.population_factors = population_factors
        self.jump_order = jump_order
        self.jump_time_graining = jump_time_graining
        self.jump_kernel_cutoff = jump_kernel_cutoff
        self.jump_kernel_zero_cutoff = jump_kernel_zero_cutoff
        self.response_diagnostics: list[dict[str, Any]] = []

        # FIXME: check the compatibility of the axes

        if system is not None:
            self.system = system
            self._has_system = True
        else:
            self._has_system = False

        if responses is not None:
            self.resp_fcions = responses
            self._has_responses = True
        else:
            self._has_responses = False

        # FIXME: properties to be protected
        self.dynamics = dynamics

        # unprotected properties
        self.data: numpy.ndarray | None = None

        self.responses: list[Any] = []

        self._relaxation_tensor = None
        self._rate_matrix = None
        self._response_rate_matrix: Any = None
        self._relaxation_hamiltonian = None
        self._population_time_axis = None
        self._has_relaxation_tensor = False
        self._has_rate_matrix = False
        self.relaxation_theory = relaxation_theory
        self.rate_matrix_time_dependent = rate_matrix_time_dependent
        self.relaxation_cutoff_time = relaxation_cutoff_time
        self.rate_matrix_options = (
            {} if rate_matrix_options is None else rate_matrix_options
        )
        if relaxation_tensor is not None:
            self._relaxation_tensor = relaxation_tensor
            self._has_relaxation_tensor = True
        if effective_hamiltonian is not None:
            self._relaxation_hamiltonian = effective_hamiltonian
        if rate_matrix is not None:
            self._rate_matrix = rate_matrix
            self._has_rate_matrix = True

        #
        # after bootstrap information
        #
        self.sys: Any = None
        self.lab: Any = None
        self.t1s: Any = None
        self.t3s: Any = None
        self.rmin: Any = None
        self.rwa: Any = None
        self.oa1: Any = None
        self.oa3: Any = None
        self.Uee: Any = None
        self.Uc0: Any = None

        self.tc = 0

    def _detection_weight(self, resp: Any) -> float:
        """Returns the detection weight for a response contribution."""
        if self.twodtype == "2DES":
            return 1.0

        if self.twodtype == "F-2DES":
            if not isinstance(resp, NonLinearResponse):
                raise NotImplementedError(
                    "F-2DES requires NonLinearResponse metadata; predefined "
                    "LiouvillePathway responses are not supported."
                )

            if self.gamma_factor is not None:
                if resp.process == "ESA":
                    return self.gamma_factor - 1.0
                return 1.0

        raise Exception("Unknown type of 2D spectrum: " + self.twodtype)

    def _vprint(self, *args: Any, **kwargs: Any) -> None:
        """Prints a string if the self.verbose attribute is True"""
        if self.verbose:
            print(*args, **kwargs)

    def bootstrap(
        self,
        rwa: float = 0.0,
        pad: int = 0,
        lab: Any = None,
        verbose: bool = False,
        write_resp: bool | str = False,
        keep_resp: bool = False,
        **kwargs: Any,
    ) -> None:
        """Sets up the environment for 2D calculation
        write_resp takes a string, creates a directory with the name of
        the string and saves the respoonses and time axis as a npz file

        keep_resp saves the responses as a list of dictionaries. The
        list goes through the time points in t2.

        """
        self.verbose = verbose
        self.pad = pad
        self.write_resp = write_resp
        self.keep_resp = keep_resp
        self.rwa = Manager().convert_energy_2_internal_u(rwa)

        with energy_units("int"):
            if self.write_resp:
                try:
                    if isinstance(write_resp, str):
                        os.mkdir(write_resp)
                except OSError:
                    print(
                        "Creation of the directory failed, "
                        "it either already exists "
                        "or you didn't give a string"
                    )

            if self._has_system:
                if isinstance(self.system, (Aggregate, OpenSystem)):
                    pass

                else:
                    raise Exception("Molecule 2D not implememted")

                sys = self.system
                sys.diagonalize()

                #
                # Construct band_system object
                #
                Nb = 3
                Ns = numpy.zeros(Nb, dtype=numpy.int32)
                Ns[0] = sys.Nb[0]  # 1
                Ns[1] = sys.Nb[1]  # agg.nmono
                Ns[2] = sys.Nb[2]  # Ns[1]*(Ns[1]-1)/2

                #
                # Relaxation rates
                #
                relaxation_requested = (
                    self._has_rate_matrix or self.relaxation_theory is not None
                )
                if relaxation_requested:
                    validate_2d_time_axes(self.t1axis, self.t2axis, self.t3axis)
                    self._population_time_axis = get_common_time_axis(
                        self.t1axis, self.t2axis, self.t3axis
                    )

                if self._has_rate_matrix:
                    KK = self._rate_matrix
                elif self.relaxation_theory is not None:
                    KK = sys.get_RateMatrix(
                        relaxation_theory=self.relaxation_theory,
                        time_dependent=self.rate_matrix_time_dependent,
                        relaxation_cutoff_time=self.relaxation_cutoff_time,
                        **self.rate_matrix_options,
                    )
                else:
                    KK = None

                # relaxation rate in single exciton band
                if KK is None:
                    Kr = None
                    self._response_rate_matrix = None
                else:
                    Kr = get_single_exciton_rate_matrix(sys, KK)  # *10.0
                    self._response_rate_matrix = Kr
                # print(1.0/KK.data)

                # FIXME: we need also 2 exciton rates
                #

                #
                # Lineshape functions
                #
                sbi = sys.get_SystemBathInteraction()
                cfm = sbi.CC
                cfm.create_double_integral()
                sys.get_lineshape_functions(self.jump_order)

                #
                #  This section will also be removed - It goes to the new Response class
                #

            ###############################################################################

            #
            # bootstrap responses
            #
            if self._has_responses:
                for rsp in self.resp_fcions:
                    rsp.set_rwa(self.rwa)

            elif self._has_system:
                # FIXME: create responses from system

                pass
                # FIXME: set _has_responses to True after they are calculated

            #
            # define lab settings
            #
            if lab is None:
                self.lab = lab_settings(lab_settings.FOUR_WAVE_MIXING)
                X = numpy.array([1.0, 0.0, 0.0], dtype=numpy.float64)
                self.lab.set_laser_polarizations(X, X, X, X)
            else:
                self.lab = lab

            #
            # Other parameters
            #
            # dt = self.t1axis.step
            self.rmin = 0.0001
            self.t1s = self.t1axis.data
            self.t3s = self.t3axis.data

            atype = self.t1axis.atype
            self.t1axis.atype = "complete"
            self.oa1 = self.t1axis.get_FrequencyAxis()
            self.oa1.data += self.rwa
            self.oa1.start += self.rwa
            # print(self.oa1.start, self.oa1.data[0])
            self.t1axis.atype = atype

            atype = self.t3axis.atype
            self.t3axis.atype = "complete"
            self.oa3 = self.t3axis.get_FrequencyAxis()
            self.oa3.data += self.rwa
            self.oa3.start += self.rwa
            # print(self.oa3.start, self.oa3.data[0])
            self.t3axis.atype = atype

            self.tc = 0

    def reset_t2_time(self) -> None:
        """Resets the population time of the calculations"""
        self.tc = 0

    def calculate_next(self) -> Any:
        """Calculate next population time of a 2D spectrum"""
        sone = self.calculate_one(self.tc)
        self.tc += 1
        return sone

    def calculate_one(self, tc: int) -> Any:
        """Calculate one population time"""
        try:
            tt2 = self.t2axis.data[tc]
        except (IndexError, AttributeError):
            print(
                "Time axis error:\n  perhaps tc =",
                tc,
                " (representing t2 population time) is outside range?",
            )
            print(
                "You can reset automatic calculation along the population time axes by calling:"
            )
            print("> twodcalc.reset_t2_time() ")
            print("where 'twodcalc' is the TwoDResponseCalculator object.")

        Nt2 = self.t2axis.length
        Nr1 = self.t1axis.length
        Nr3 = self.t3axis.length

        # FIXME: on which axis we should be looking for it2 ???
        (it2, err) = self.t2axis.locate(tt2)
        self._vprint(
            "t2 = " + str(tt2) + "fs (it2 = " + str(it2) + " of " + str(Nt2) + ")",
            end="\r",
        )

        #
        # Initialize response storage
        #
        # if _have_aceto and self._has_system:
        #    order = 'F'
        # else:
        order: Literal["C", "F"] = "C"

        ntype = numpy.complex128

        # FIXME:  Fix the axis of time
        # the order of axis is wrong. 2D code works only if it is Nr3, Nr1

        resp_r = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_n = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Rgsb = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Ngsb = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)

        resp_Rse = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Nse = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Resa = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Nesa = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Rsewt = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Nsewt = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Resawt = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        resp_Nesawt = numpy.zeros((Nr1, Nr3), dtype=ntype, order=order)
        response_pieces: list[tuple[Any, numpy.ndarray]] = []

        if self._has_system and not self._has_responses:
            #
            # Calculating all responses from the system
            #
            self.resp_fcions = []
            response_kwargs: dict[str, Any] = dict(
                rate_matrix=self._response_rate_matrix,
                population_time_axis=self._population_time_axis,
                jump_time_graining=self.jump_time_graining,
                jump_kernel_cutoff=self.jump_kernel_cutoff,
                jump_kernel_zero_cutoff=self.jump_kernel_zero_cutoff,
            )

            # basic pathways
            Nr1g = NonLinearResponse(
                self.lab,
                self.system,
                "R1g",
                self.t1axis,
                self.t2axis,
                self.t3axis,
                **response_kwargs,
            )
            Nr2g = NonLinearResponse(
                self.lab,
                self.system,
                "R2g",
                self.t1axis,
                self.t2axis,
                self.t3axis,
                **response_kwargs,
            )
            Nr3g = NonLinearResponse(
                self.lab,
                self.system,
                "R3g",
                self.t1axis,
                self.t2axis,
                self.t3axis,
                **response_kwargs,
            )
            Nr4g = NonLinearResponse(
                self.lab,
                self.system,
                "R4g",
                self.t1axis,
                self.t2axis,
                self.t3axis,
                **response_kwargs,
            )

            self.resp_fcions.append(Nr1g)
            self.resp_fcions.append(Nr2g)
            self.resp_fcions.append(Nr3g)
            self.resp_fcions.append(Nr4g)

            if self.system.mult > 1:
                # ESA (if mult > 1)
                Nr1f = NonLinearResponse(
                    self.lab,
                    self.system,
                    "R1f",
                    self.t1axis,
                    self.t2axis,
                    self.t3axis,
                    **response_kwargs,
                )
                Nr2f = NonLinearResponse(
                    self.lab,
                    self.system,
                    "R2f",
                    self.t1axis,
                    self.t2axis,
                    self.t3axis,
                    **response_kwargs,
                )

                self.resp_fcions.append(Nr1f)
                self.resp_fcions.append(Nr2f)

            # relaxation (if relax neq 0)
            Nr1g_scM0g = NonLinearResponse(
                self.lab,
                self.system,
                "R1g_scM0g",
                self.t1axis,
                self.t2axis,
                self.t3axis,
                **response_kwargs,
            )
            Nr2g_scM0g = NonLinearResponse(
                self.lab,
                self.system,
                "R2g_scM0g",
                self.t1axis,
                self.t2axis,
                self.t3axis,
                **response_kwargs,
            )

            KK = Nr1g_scM0g.KK
            if numpy.all(numpy.isclose(KK, 0.0, atol=1e-9)):
                pass  # we avoid calculating relaxation if the matrix is zero

            else:
                # print("Including relaxation")
                self.resp_fcions.append(Nr1g_scM0g)
                self.resp_fcions.append(Nr2g_scM0g)
                if Nr1g_scM0g._uses_single_jump_storage():
                    self.resp_fcions.append(
                        NonLinearResponse(
                            self.lab,
                            self.system,
                            "R1g_scM1g",
                            self.t1axis,
                            self.t2axis,
                            self.t3axis,
                            **response_kwargs,
                        )
                    )
                    self.resp_fcions.append(
                        NonLinearResponse(
                            self.lab,
                            self.system,
                            "R2g_scM1g",
                            self.t1axis,
                            self.t2axis,
                            self.t3axis,
                            **response_kwargs,
                        )
                    )

            if self.system.mult > 1:
                Nr1f_scM0g = NonLinearResponse(
                    self.lab,
                    self.system,
                    "R1f_scM0g",
                    self.t1axis,
                    self.t2axis,
                    self.t3axis,
                    **response_kwargs,
                )
                Nr2f_scM0g = NonLinearResponse(
                    self.lab,
                    self.system,
                    "R2f_scM0g",
                    self.t1axis,
                    self.t2axis,
                    self.t3axis,
                    **response_kwargs,
                )
                Nr1f_scM0e = NonLinearResponse(
                    self.lab,
                    self.system,
                    "R1f_scM0e",
                    self.t1axis,
                    self.t2axis,
                    self.t3axis,
                    **response_kwargs,
                )
                Nr2f_scM0e = NonLinearResponse(
                    self.lab,
                    self.system,
                    "R2f_scM0e",
                    self.t1axis,
                    self.t2axis,
                    self.t3axis,
                    **response_kwargs,
                )

                self.resp_fcions.append(Nr1f_scM0g)
                self.resp_fcions.append(Nr2f_scM0g)
                self.resp_fcions.append(Nr1f_scM0e)
                self.resp_fcions.append(Nr2f_scM0e)
                if Nr1f_scM0g._uses_single_jump_storage():
                    self.resp_fcions.append(
                        NonLinearResponse(
                            self.lab,
                            self.system,
                            "R1f_scM1g",
                            self.t1axis,
                            self.t2axis,
                            self.t3axis,
                            **response_kwargs,
                        )
                    )
                    self.resp_fcions.append(
                        NonLinearResponse(
                            self.lab,
                            self.system,
                            "R2f_scM1g",
                            self.t1axis,
                            self.t2axis,
                            self.t3axis,
                            **response_kwargs,
                        )
                    )
                    self.resp_fcions.append(
                        NonLinearResponse(
                            self.lab,
                            self.system,
                            "R1f_scM1e",
                            self.t1axis,
                            self.t2axis,
                            self.t3axis,
                            **response_kwargs,
                        )
                    )
                    self.resp_fcions.append(
                        NonLinearResponse(
                            self.lab,
                            self.system,
                            "R2f_scM1e",
                            self.t1axis,
                            self.t2axis,
                            self.t3axis,
                            **response_kwargs,
                        )
                    )

            self._has_responses = True

        if self._has_responses:
            #
            # Calculation from predefined non-linear responses
            #

            for resp in self.resp_fcions:
                if isinstance(resp, NonLinearResponse):
                    data = resp.calculate_matrix(tt2)
                    self.response_diagnostics.append(dict(resp.diagnostics))
                    response_pieces.append((resp, data))
                    if resp.rtype == "R":
                        if resp.process == "GSB":
                            resp_Rgsb += data
                        elif resp.process == "SE":
                            if resp.is_transfer:
                                resp_Rsewt += data
                            else:
                                resp_Rse += data
                        elif resp.process == "ESA":
                            if resp.is_transfer:
                                resp_Resawt += data
                            else:
                                resp_Resa += data
                        else:
                            raise Exception("Unknown response process")

                    elif resp.rtype == "NR":
                        if resp.process == "GSB":
                            resp_Ngsb += data
                        elif resp.process == "SE":
                            if resp.is_transfer:
                                resp_Nsewt += data
                            else:
                                resp_Nse += data
                        elif resp.process == "ESA":
                            if resp.is_transfer:
                                resp_Nesawt += data
                            else:
                                resp_Nesa += data
                        else:
                            raise Exception("Unknown response process")

                    else:
                        raise Exception("Unknown response type")

                elif isinstance(resp, LiouvillePathway):
                    resp_any: Any = resp
                    if resp.rtype == "R":
                        data = resp_any.calculate_matrix(
                            self.lab, None, tt2, self.t1s, self.t3s, self.rwa
                        )
                        response_pieces.append((resp, data))
                        resp_Rgsb += data

                    elif resp.rtype == "NR":
                        data = resp_any.calculate_matrix(
                            self.lab, None, tt2, self.t1s, self.t3s, self.rwa
                        )
                        response_pieces.append((resp, data))
                        resp_Ngsb += data
                    else:
                        raise Exception("Unknown response type")

        else:
            raise Exception("Calculation method not implemented")

        # only for Aceto we need the sum
        #
        # FIXME: discontinue Aceto and remove the sum (and the code above)
        #
        resp_r = resp_Rgsb + resp_Rse + resp_Resa + resp_Rsewt + resp_Resawt
        resp_n = resp_Ngsb + resp_Nse + resp_Nesa + resp_Nsewt + resp_Nesawt

        #
        # Calculate corresponding 2D spectrum
        #
        onetwod = TwoDResponse()

        # pad is set to 0 by default. If changed in the bootstrap,
        # responses are padded with 0s and the time axis is lengthened
        t13Pad = TimeAxis(
            self.t1axis.start, self.t1axis.length + self.pad, self.t1axis.step
        )
        response_window = None
        if self.pad > 0:
            self._vprint("padding by - " + str(self.pad))

            t13Pad.atype = "complete"
            t13PadFreq = t13Pad.get_FrequencyAxis()
            t13PadFreq.data += self.rwa
            t13PadFreq.start += self.rwa

            onetwod.set_axis_1(t13PadFreq)
            onetwod.set_axis_3(t13PadFreq)

            # Sloping the end of the data down to 0 so there isn't a hard
            # cutoff at the end of the data
            from scipy.signal import windows as sig

            window = 20
            response_window = sig.tukey(window * 2, 1, sym=False)

            resp_r = _pad_response_data(resp_r, self.pad, response_window)
            resp_n = _pad_response_data(resp_n, self.pad, response_window)
            resp_Rgsb = _pad_response_data(resp_Rgsb, self.pad, response_window)
            resp_Ngsb = _pad_response_data(resp_Ngsb, self.pad, response_window)
            resp_Rse = _pad_response_data(resp_Rse, self.pad, response_window)
            resp_Nse = _pad_response_data(resp_Nse, self.pad, response_window)
            resp_Resa = _pad_response_data(resp_Resa, self.pad, response_window)
            resp_Nesa = _pad_response_data(resp_Nesa, self.pad, response_window)
            resp_Rsewt = _pad_response_data(resp_Rsewt, self.pad, response_window)
            resp_Nsewt = _pad_response_data(resp_Nsewt, self.pad, response_window)
            resp_Resawt = _pad_response_data(resp_Resawt, self.pad, response_window)
            resp_Nesawt = _pad_response_data(resp_Nesawt, self.pad, response_window)

        else:
            onetwod.set_axis_1(self.oa1)
            onetwod.set_axis_3(self.oa3)

        # FIXME: Make a decision, if this is to be kept
        # Right now the code does not distinguish different response types, except rephasing and non-rephasing
        # If we decide to remove the detained storage, we can remove a lot of functionality from TwoDResponse.
        # This would discart a let of information useful for inspection.
        #
        # Likely the best solution, is the allow storage of details, only if the user asks
        #
        if self.keep_resp:
            resp = {
                "time": self.t1axis.data,
                "time_pad": t13Pad.data,
                "rTot": resp_r,
                "nTot": resp_n,
                "rGSB": resp_Rgsb,
                "nGSB": resp_Ngsb,
                "rSE": resp_Rse,
                "nSE": resp_Nse,
                "rESA": resp_Resa,
                "nESA": resp_Nesa,
                "rSEWT": resp_Rsewt,
                "nSEWT": resp_Nsewt,
                "rESAWT": resp_Resawt,
                "nESAWT": resp_Nesawt,
            }
            self.responses.append(resp)

        if self.write_resp and isinstance(self.write_resp, str):
            numpy.savez(
                "./" + self.write_resp + "/respT" + str(int(tt2)) + ".npz",
                time=self.t1axis.data,
                time_pad=t13Pad.data,
                rTot=resp_r,
                nTot=resp_n,
                rGSB=resp_Rgsb,
                nGSB=resp_Ngsb,
                rSE=resp_Rse,
                nSE=resp_Nse,
                rESA=resp_Resa,
                nESA=resp_Nesa,
                rSEWT=resp_Rsewt,
                nSEWT=resp_Nsewt,
                rESAWT=resp_Resawt,
                nESAWT=resp_Nesawt,
            )

        for resp, data in response_pieces:
            data = _pad_response_data(data, self.pad, response_window)

            if isinstance(resp, NonLinearResponse):
                signal = resp.signal
                dtype = resp.storage_type
                resolution = "types"
            elif isinstance(resp, LiouvillePathway):
                signal = signal_REPH if resp.rtype == "R" else signal_NONR
                dtype = signal
                resolution = "signals"
            else:
                raise Exception("Unknown response object")

            spect_data = _fourier_transform_response(data, signal)
            spect_data *= self._detection_weight(resp)
            onetwod._add_data(spect_data, resolution=resolution, dtype=dtype)

        onetwod.set_t2(self.t2axis.data[tc])

        return onetwod

    def calculate(self) -> Any:
        """Returns 2D spectrum

        Calculates and returns TwoDResponseContainer containing 2D spectrum
        based on the parameters specified in this object.


        """
        # FIXME: we will later use only one branch below
        from .twodcontainer import TwoDResponseContainer

        self.response_diagnostics = []

        # if _have_aceto and self._has_system:

        #     twods = TwoDResponseContainer(self.t2axis)

        #     teetoos = self.t2axis.data

        #     kk = 0
        #     Nk = teetoos.shape[0]

        #     for tt2 in teetoos:

        #         self._vprint_r(" Calculating t2 =", tt2, "fs (",kk,"of",Nk,")")
        #         onetwod = self.calculate_next()
        #         twods.set_spectrum(onetwod)
        #         kk += 1

        #     return twods

        if self._has_responses or self._has_system:
            # calculate user defined responses

            twods = TwoDResponseContainer(self.t2axis)

            teetoos = self.t2axis.data

            kk = 0
            Nk = teetoos.shape[0]
            for tt2 in teetoos:
                onetwod = self.calculate_next()
                twods.set_spectrum(onetwod)
                kk += 1

            return twods

        raise Exception("2D calculation in this mode not implemented.")

    def reset_evaluation_functions(self, fcions: list[Any]) -> None:
        """Resets the evaluation functions used by the reseponse functions"""
        if self._has_responses:
            for rsp, fce in zip(self.resp_fcions, fcions):
                rsp.set_evaluation_function(fce)
                print(rsp)

        else:
            raise Exception("Calculatore has no responses defined.")
