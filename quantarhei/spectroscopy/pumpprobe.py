from __future__ import annotations

import copy
from typing import Any

import matplotlib.pyplot as plt
import numpy
import scipy

from .. import COMPLEX, signal_TOTL
from ..builders.aggregates import Aggregate
from ..builders.molecules import Molecule
from ..core.dfunction import DFunction
from ..core.frequency import FrequencyAxis
from ..core.managers import Manager, energy_units
from ..core.time import TimeAxis
from ..core.units import convert
from ..utils import derived_type
from .mocktwodcalculator import MockTwoDResponseCalculator as MockTwoDSpectrumCalculator
from .responses import NonLinearResponse
from .twodcontainer import TwoDSpectrumContainer


class PumpProbeSpectrum(DFunction):
    """Class representing a pump-probe spectrum.

    A one-dimensional spectrum obtained at a fixed waiting time ``t2``
    from a pump-probe experiment. Data and axis are set via ``set_data``
    and ``set_axis`` after construction.
    """

    # data = None

    def __init__(self) -> None:
        self.t2 = -1.0
        self.data = None
        self._has_imag = False

    def set_axis(self, axis: Any) -> None:
        self.xaxis = axis
        self.axis = self.xaxis

    def set_data(self, data: Any) -> None:
        self.data = data

    def get_PumpProbeSpectrum(self) -> PumpProbeSpectrum:
        """Returns self

        This method is here to override the one inherented from TwoDSpectrum

        """
        return self

    def set_t2(self, t2: float) -> None:
        """Sets the t2 (waiting time) of the spectrum"""
        self.t2 = t2

    def get_t2(self) -> float:
        """Returns the t2 (waiting time) of the spectrum"""
        return self.t2

    def _add_data(self, data: Any) -> None:
        if self.data is None:
            self.set_data(data)
        else:
            if self.data.size != len(data):
                raise OSError("Added data length does not match the currentone")
            self.data += data

    # FIXME: Add function _add_data (if data None = set_data, else add)


class _RWAOverrideSystem:
    """Delegates to a system while overriding the response-backend RWA."""

    def __init__(self, system: Any, rwa: Any, lineshape_timeaxis: Any = None) -> None:
        self._system = system
        self._rwa = rwa
        self._lineshape_timeaxis = lineshape_timeaxis

    def get_RWA_suggestion(self) -> Any:
        return self._rwa

    def get_lineshape_functions(self, config: dict | int | None = None) -> Any:
        return self._system.get_lineshape_functions(
            config=config, timeaxis=self._lineshape_timeaxis
        )

    def __getattr__(self, name: str) -> Any:
        return getattr(self._system, name)


class PumpProbeSpectrumContainer(TwoDSpectrumContainer):
    """Container for a set of pump-probe spectra indexed by waiting time.

    Parameters
    ----------
    t2axis : TimeAxis or None, optional
        Waiting-time axis shared by all stored spectra. Default is ``None``.
    """

    def __init__(self, t2axis: Any = None) -> None:

        self.t2axis = t2axis
        self.spectra: dict[Any, Any] = {}

    def plot(self) -> None:

        plt.clf()
        spctr = self.get_spectra()
        for sp in spctr:
            plt.plot(sp.xaxis.data, sp.data)

    def set_spectrum(self, spec: Any, tag: Any = None) -> None:
        self.spectra[tag] = spec

    def amax(self, spart: Any = None) -> Any:
        mxs = []
        for s in self.get_spectra():
            spect = numpy.real(s.data)
            mx = numpy.amax(spect)
            mxs.append(mx)
        return numpy.amax(numpy.array(mxs))

    def amin(self) -> Any:
        mxs = []
        for s in self.get_spectra():
            spect = numpy.real(s.data)
            mx = numpy.amin(spect)
            mxs.append(mx)
        return numpy.amin(numpy.array(mxs))

    def plot2D(
        self,
        axis: Any = None,
        units: str = "nm",
        zero_centered: bool = True,
        lines: Any = None,
    ) -> None:

        t2ax = self.t2axis
        freqax = self.spectra[t2ax.data[0]].axis
        # use only positive frequency axis
        min_ind = numpy.min(numpy.where(freqax.data > 0.0))
        with energy_units(units):
            # Prepare array of pump-probe spectra
            ppspec2D = numpy.zeros((t2ax.length, freqax.length))
            count = 0
            for T2 in t2ax.data:
                ppspec2D[count] += self.spectra[T2].data
                count += 1

            # plot pump-probe spectra
            X, Y = numpy.meshgrid(t2ax.data, freqax.data[min_ind:])
            fig, ax = plt.subplots(figsize=(18, 9))

            if zero_centered:
                p = ax.contourf(
                    X,
                    Y,
                    ppspec2D.T[min_ind:],
                    60,
                    cmap=plt.cm.jet,  # type: ignore[attr-defined]
                    vmin=-abs(ppspec2D).max(),
                    vmax=abs(ppspec2D).max(),
                )  # , vmin=ppspec2D.min(), vmax=ppspec2D.max())
            else:
                p = ax.contourf(
                    X,
                    Y,
                    ppspec2D.T[min_ind:],
                    60,
                    cmap=plt.cm.jet,  # type: ignore[attr-defined]
                    vmin=ppspec2D.min(),
                    vmax=ppspec2D.max(),
                )

            if isinstance(lines, list):
                for line_freq in lines:
                    plt.plot(
                        t2ax.data, numpy.ones(t2ax.length) * line_freq, "w", linewidth=2
                    )

            if axis is not None:
                ax.axes.set_ylim(axis[1][0], axis[1][1])
                ax.axes.set_xlim(axis[0][0], axis[0][1])
            else:
                ax.axes.set_ylim(400, 750)
                ax.axes.set_xlim(0, 500)
            plt.xlabel("Time [fs]")
            plt.ylabel("Wavelength [" + units + "]")

            ax.legend()
            plt.colorbar(p, ax=ax)
            fig.savefig("PP_2D_spectra.png", format="png", dpi=1200)

    def plot_slices(
        self, freqs: Any, expRes: Any = None, units: str = "nm"
    ) -> numpy.ndarray:

        # Initialize frequency cuts of PP spectra
        spectra_freq = numpy.zeros((len(freqs), self.t2axis.length), dtype="f8")
        count = 0
        for T2 in self.t2axis.data:
            sp = self.spectra[T2]
            sp._has_imag = False
            sp._set_splines()
            freq_int = convert(freqs, units, "int")
            spectra_freq[:, count] = self.spectra[T2].at(freq_int)
            count += 1

        with energy_units(units):
            fig = plt.figure(figsize=(18, 9))
            Nsp = len(freqs)
            plt_num = numpy.arange(Nsp) + 1
            plt_num = plt_num.reshape((Nsp // 3, 3)).T.reshape(Nsp)
            plt.xlabel("Delay time [ps]")
            plt.ylabel("Intensity [arb. u.]")
            for ii in range(Nsp):
                plt.subplot(Nsp // 3, 3, plt_num[ii])
                plt.title(str(numpy.round(freqs[ii])) + " " + units)
                plt.plot(self.t2axis.data / 1000, spectra_freq[ii])
                if expRes is not None:
                    PP_exp_freq = expRes["frequency"]
                    PP_exp_spec = expRes["intensity"]
                    plt.plot(PP_exp_freq[ii], PP_exp_spec[ii])

                plt.xscale("symlog", linthreshx=(0.1))
            plt.subplots_adjust(hspace=0.6, wspace=0.6)
            plt.show()
            fig.savefig("PP_slices_spectra.png", format="png", dpi=1200)
        return spectra_freq

    def make_movie(
        self,
        filename: str,
        window: Any = None,
        stype: Any = None,
        spart: Any = None,
        cmap: Any = None,
        Npos_contours: int = 10,
        vmax: Any = None,
        vmin_ratio: float = 0.5,
        xlabel: Any = None,
        ylabel: Any = None,
        axis_label_font: Any = None,
        start: Any = None,
        end: Any = None,
        frate: int = 20,
        dpi: int = 100,
        show_states: Any = None,
        show_states_func: Any = None,
        label: Any = None,
        label_func: Any = None,
        text_loc: Any = None,
        progressbar: bool = False,
        use_t2: bool = False,
        title: str = "",
        comment: str = "",
        axis: Any = None,
        vmin: Any = None,
        **kwargs: Any,
    ) -> None:

        import matplotlib.animation as manimation
        import matplotlib.pyplot as plt

        FFMpegWriter = manimation.writers["ffmpeg"]

        metadata = dict(
            title="Test Movie", artist="Matplotlib", comment="Movie support!"
        )
        writer = FFMpegWriter(fps=frate, metadata=metadata)

        fig = plt.figure()

        spctr = self.get_spectra()
        l = len(spctr)
        last_t2 = spctr[l - 1].get_t2()
        first_t2 = spctr[0].get_t2()

        if vmax is None:
            mx = self.amax()
        else:
            mx = vmax

        if vmin is None:
            mn = self.amin()
        else:
            mn = vmin

        mxabs = max(numpy.abs(mx), numpy.abs(mn))
        mx = mx + 0.05 * mxabs
        mn = mn - 0.05 * mxabs

        if start is None:
            start = first_t2
        if end is None:
            end = last_t2

        with writer.saving(fig, filename, dpi):
            k = 0
            # Initial call to print 0% progress
            sp2write = self.get_spectra(start=start, end=end)
            l = len(sp2write)
            if progressbar:
                self._printProgressBar(
                    0, l, prefix="Progress:", suffix="Complete", length=50
                )
            for sp in self.get_spectra(start=start, end=end):
                # FIXME: this does not work as it should yet
                sp.plot(
                    show=False, fig=fig, axis=axis, label="T=" + str(sp.get_t2()) + "fs"
                )  # , vmax=mx, vmin=mn,
                # )
                writer.grab_frame()
                if progressbar:
                    self._printProgressBar(
                        k + 1, l, prefix="Progress:", suffix="Complete", length=50
                    )

                k += 1


class PumpProbeSpectrumCalculator:
    """Calculator for pump-probe spectra derived from 2D response pathways.

    Parameters
    ----------
    t2axis : TimeAxis
        Waiting-time axis over which spectra are calculated.
    t3axis : TimeAxis
        Detection-time axis used to build the frequency axis.
    system : Molecule or Aggregate, optional
        Quantum system for which spectra are calculated.
    dynamics : str, optional
        Dynamics type, e.g. ``"secular"``. Default is ``"secular"``.
    relaxation_tensor : optional
        Relaxation tensor to include population relaxation.
    rate_matrix : optional
        Rate matrix alternative to a full relaxation tensor.
    effective_hamiltonian : optional
        Effective Hamiltonian used in the relaxation basis.
    separate_relax_pwy : bool, optional
        If ``True``, pathways with relaxation steps are treated separately.
        Default is ``True``.
    """

    t2axis = derived_type("t2axis", TimeAxis)
    t3axis = derived_type("t3axis", TimeAxis)

    system = derived_type("system", [Molecule, Aggregate])

    def __init__(
        self,
        t2axis: Any,
        t3axis: Any,
        system: Any = None,
        dynamics: str = "secular",
        relaxation_tensor: Any = None,
        rate_matrix: Any = None,
        effective_hamiltonian: Any = None,
        separate_relax_pwy: bool = True,
        population_propagator: Any = None,
        density_matrix_propagator: Any = None,
        density_matrix_trajectory: Any = None,
        population_time_axis: Any = None,
        include_nonsecular_remainder: bool = True,
        include_remainder: bool = True,
        dipole_normalization_tol: float = 1.0e-12,
    ) -> None:

        self.t2axis = t2axis
        self.t3axis = t3axis

        external_count = sum(
            item is not None
            for item in (
                population_propagator,
                density_matrix_propagator,
                density_matrix_trajectory,
            )
        )
        if external_count > 1:
            raise ValueError(
                "population_propagator, density_matrix_propagator, and "
                "density_matrix_trajectory are mutually exclusive"
            )

        if system is not None:
            self.system = copy.deepcopy(system)

        # FIXME: properties to be protected
        self.dynamics = dynamics

        # whether to treat differently pwy with jumps or not
        self._separate_relax = separate_relax_pwy

        # unprotected properties
        self.data = None

        self._relaxation_tensor = None
        self._rate_matrix = None
        self._relaxation_hamiltonian = None
        self._has_relaxation_tensor = False
        self._is_adiabatic = False
        self._adiabatic_noBath = False
        if relaxation_tensor is not None:
            self._relaxation_tensor = relaxation_tensor
            self._has_relaxation_tensor = True
        if effective_hamiltonian is not None:
            self._relaxation_hamiltonian = effective_hamiltonian
        if rate_matrix is not None:
            self._rate_matrix = rate_matrix
            self._has_rate_matrix = True
        self.goft_matrix: Any = None
        self.reorg_matrix: Any = None

        # after bootstrap information
        self.lab: Any = None
        self.t3s: Any = None
        self.pathways: Any = None
        self.rwa: Any = None
        self.density_matrix_trajectory: Any = density_matrix_trajectory
        self.population_propagator: Any = population_propagator
        self.density_matrix_propagator: Any = density_matrix_propagator
        self.population_time_axis: Any = population_time_axis
        self.include_nonsecular_remainder = include_nonsecular_remainder
        self.include_remainder = include_remainder
        self.dipole_normalization_tol = dipole_normalization_tol

        self.verbose = False

        self.tc = 0

    def bootstrap(
        self,
        rwa: float = 0.0,
        pathways: Any = None,
        lab: Any = None,
        verbose: bool = False,
        adiabatic: Any = None,
    ) -> None:
        """Sets up the environment for pump-probe calculation"""
        self.verbose = verbose

        if isinstance(self.system, Aggregate):
            pass

        else:
            raise Exception("Molecule pump-probe not implememted")

        self.verbose = verbose
        self.rwa = Manager().convert_energy_2_internal_u(rwa)
        self.pathways = pathways

        if adiabatic is not None:
            if adiabatic != False:
                self._is_adiabatic = True
            else:
                self._is_adiabatic = False

            if (adiabatic == "SubtractBath") or (adiabatic == "NoBath"):
                self._adiabatic_noBath = True

        with energy_units("int"):
            # atype = self.t3axis.atype
            # self.t3axis.atype = 'complete'
            # self.oa3 = self.t3axis.get_FrequencyAxis()
            # self.oa3.data += self.rwa
            # self.oa3.start += self.rwa
            # self.t3axis.atype = atype

            # we only want to retain the upper half of the spectrum
            freq = self.t3axis.get_FrequencyAxis()
            freq.data += self.rwa
            Nt = len(freq.data) // 2
            do = freq.data[1] - freq.data[0]
            st = freq.data[Nt // 2]
            # we represent the Frequency axis anew
            self.oa3 = FrequencyAxis(st, Nt, do)

        self.tc = 0
        self.lab = lab

    def set_pathways(self, pathways: Any) -> None:
        self.pathways = pathways

    def set_density_matrix_trajectory(
        self, density_matrix_trajectory: Any, timeaxis: Any = None
    ) -> None:
        """Set an externally calculated density-matrix trajectory."""
        self.density_matrix_trajectory = density_matrix_trajectory
        self.population_propagator = None
        self.density_matrix_propagator = None
        self.population_time_axis = timeaxis

    def set_population_propagator(
        self, population_propagator: Any, timeaxis: Any = None
    ) -> None:
        """Set an externally calculated one-exciton population propagator."""
        self.population_propagator = population_propagator
        self.density_matrix_trajectory = None
        self.density_matrix_propagator = None
        self.population_time_axis = timeaxis

    def set_density_matrix_propagator(
        self, density_matrix_propagator: Any, timeaxis: Any = None
    ) -> None:
        """Set an externally calculated one-exciton density-matrix propagator."""
        self.density_matrix_propagator = density_matrix_propagator
        self.density_matrix_trajectory = None
        self.population_propagator = None
        self.population_time_axis = timeaxis

    def set_dynamics(
        self,
        density_matrix_trajectory: Any = None,
        population_propagator: Any = None,
        density_matrix_propagator: Any = None,
        timeaxis: Any = None,
    ) -> None:
        """Set exactly one external dynamics object for response-backend PP."""
        dynamics_count = sum(
            item is not None
            for item in (
                density_matrix_trajectory,
                population_propagator,
                density_matrix_propagator,
            )
        )
        if dynamics_count != 1:
            raise ValueError(
                "Exactly one of density_matrix_trajectory, population_propagator, "
                "or density_matrix_propagator has to be supplied"
            )
        if density_matrix_trajectory is not None:
            self.set_density_matrix_trajectory(density_matrix_trajectory, timeaxis)
        elif population_propagator is not None:
            self.set_population_propagator(population_propagator, timeaxis)
        else:
            self.set_density_matrix_propagator(density_matrix_propagator, timeaxis)

    def bath_reorg(self, cfm: Any, indx: Any) -> float:
        coft = cfm.cfuncs[cfm.get_index_by_where((indx, indx))]
        reorg_bath = 0.0
        for parm in coft.params:
            if parm["ftype"] == "OverdampedBrownian":
                reorg_bath += parm["reorg"]
        return reorg_bath

    def _excitonic_reorg_diag(
        self, SS: numpy.ndarray, subtract_bath: bool = True
    ) -> numpy.ndarray:
        """Returns the reorganisation energy of an exciton state"""
        # SystemBathInteraction
        sbi = self.system.get_SystemBathInteraction()
        AG = self.system
        # CorrelationFunctionMatrix
        cfm = sbi.CC

        reorg_exct = numpy.zeros(AG.Ntot)

        # electronic states corresponding to single excited states
        elst = numpy.where(AG.which_band == 1)[0]
        for n in range(1, AG.Nb[1] + 1):
            for el1 in elst:
                if subtract_bath:
                    reorgB = self.bath_reorg(cfm, el1 - 1)
                else:
                    reorgB = 0.0

                reorg = cfm.get_reorganization_energy(el1 - 1, el1 - 1) - reorgB
                for kk in AG.vibindices[el1]:
                    reorg_exct[n] += (SS[kk, n] ** 2) * (SS[kk, n] ** 2) * reorg

        elst_dbl1 = numpy.where(AG.which_band == 2)[0]
        elst_dbl2 = numpy.where(AG.which_band == 2)[0]
        for n in range(AG.Nb[0] + AG.Nb[1], AG.Ntot):
            for el1 in elst_dbl1:
                for el2 in elst_dbl2:
                    if subtract_bath:
                        reorgB = self.bath_reorg(cfm, el1 - 1)
                    else:
                        reorgB = 0.0

                    reorg = cfm.get_reorganization_energy(el1 - 1, el2 - 1) - reorgB

                    for kk in AG.vibindices[el1]:
                        for ll in AG.vibindices[el2]:
                            reorg_exct[n] += (SS[kk, n] ** 2) * (SS[ll, n] ** 2) * reorg

        return reorg_exct

    def _site_reorg_diag(self, subtract_bath: bool = True) -> numpy.ndarray:
        """Returns the reorganisation energy of an exciton state"""
        # SystemBathInteraction
        sbi = self.system.get_SystemBathInteraction()
        AG = self.system
        # CorrelationFunctionMatrix
        cfm = sbi.CC

        reorg_site = numpy.zeros(AG.Ntot)

        # electronic states corresponding to single excited states
        elst = numpy.where(AG.which_band == 1)[0]
        for el1 in elst:
            if subtract_bath:
                reorgB = self.bath_reorg(cfm, el1 - 1)
            else:
                reorgB = 0.0

            reorg = cfm.get_reorganization_energy(el1 - 1, el1 - 1) - reorgB
            for kk in AG.vibindices[el1]:
                reorg_site[kk] += reorg

        elst_dbl1 = numpy.where(AG.which_band == 2)[0]
        for el1 in elst_dbl1:
            if subtract_bath:
                reorgB = self.bath_reorg(cfm, el1 - 1)
            else:
                reorgB = 0.0

            reorg = cfm.get_reorganization_energy(el1 - 1, el1 - 1) - reorgB

            for kk in AG.vibindices[el1]:
                reorg_site[kk] += reorg

        return reorg_site

    def calculate_all_system_approx(
        self,
        sys: Any,
        rdmt: Any,
        lab: Any,
        show_progress: bool = False,
        approx: Any = None,
        spec: Any = None,
        gsb_mode: str = "hole",
        orientational_averaging: str = "legacy",
    ) -> Any:
        """Calculates all 2D spectra for a system and reduced density matrix
        evolution. The approach assumes no diiference between pathways with
        jumps and without the jumps.

        """
        if spec is None:
            spec = ["Full"]
        # Check if the magic angle polarization is used
        if not numpy.isclose(lab.F4eM4[1:], [0, 0], atol=1e-6).all():
            message = (
                "Lab is not set to the magic angle polarization which is"
                "the only supported measurement setting for this calculation. Set "
                "polarization angle between first two pulses and the last two to"
                "54.7356103 deg. and repeat the calculation."
            )
            raise OSError(message)

        self.system = copy.deepcopy(sys)

        if self.system._diagonalized and self._is_adiabatic:
            raise Warning(
                "Not possible to use adiabatic eigenstate with diagonalized afggregate"
            )

        if self._is_adiabatic:
            # SS = numpy.identity(self.system.Ntot)
            reorg_site = self._site_reorg_diag(subtract_bath=self._adiabatic_noBath)
            # reorg_site = self._excitonic_reorg_diag(SS, subtract_bath=self._adiabatic_noBath)
            # print("site reorg :",numpy.isclose(reorg_site2,reorg_site).all())
            # print("site reorg2:",reorg_site2==reorg_site)
            for kk in range(self.system.Ntot):
                self.system.HH[kk, kk] -= reorg_site[kk]

        if not self.system._diagonalized:
            self.system.diagonalize()

        if self._is_adiabatic:
            # get exciton reorganization energy
            reorg_excit = self._excitonic_reorg_diag(
                self.system.SS, subtract_bath=self._adiabatic_noBath
            )
            # reorg_excit = self._site2excit_reorg(reorg_site,self.system.SS)

            # shift the diagonal of the exciton hamiltonian
            for ii in range(self.system.Ntot):
                self.system.HH[ii, ii] += reorg_excit[ii]

        tcont = PumpProbeSpectrumContainer(t2axis=self.t2axis)

        kk = 0
        Nk = self.t2axis.length

        rdm0 = rdmt.data[0, :, :].copy()

        printProgressBar(
            kk, Nk, prefix="     - Progress:", suffix="Complete", length=50
        )

        for T2 in self.t2axis.data:
            if show_progress:
                print(" - calculating", kk, "of", Nk, "at t2 =", T2, "fs")

            rdm = rdmt.data[kk, :, :].copy()

            if approx == "Novoderezhkin":
                ppspec1 = self.calculate_pathways_rdm_novoderezhkin(
                    rdm0,
                    rdm,
                    T2,
                    lab,
                    ptol=1.0e-6,
                    spec=spec,
                    gsb_mode=gsb_mode,
                    orientational_averaging=orientational_averaging,
                )
            else:
                ppspec1 = self.calculate_pathways_rdm(
                    rdm0,
                    rdm,
                    T2,
                    lab,
                    spec=spec,
                    gsb_mode=gsb_mode,
                    orientational_averaging=orientational_averaging,
                )

            tcont.set_spectrum(ppspec1, tag=T2)

            kk += 1
            printProgressBar(
                kk, Nk, prefix="     - Progress:", suffix="Complete", length=50
            )

        return tcont

    def _response_backend_trace(self, diagrams: list[str], tau: float, lab: Any) -> Any:
        """Calculate selected response-backend diagrams as a t3 trace at t1 = 0."""
        t1axis = TimeAxis(0.0, 1, self.t3axis.step)
        backend_system = _RWAOverrideSystem(
            self.system, self.rwa, lineshape_timeaxis=[t1axis, self.t3axis]
        )

        response = numpy.zeros(self.t3axis.length, dtype=numpy.complex128)
        for diagram in diagrams:
            rsp = NonLinearResponse(
                lab,
                backend_system,
                diagram,
                t1axis,
                self.t2axis,
                self.t3axis,
            )
            response += numpy.ravel(numpy.asarray(rsp.calculate_matrix(tau)))

        return response

    def _full_dipole_rdm_weight(
        self, rdm: Any, ii: int, jj: int, tol: float = 1.0e-12
    ) -> Any:
        """Return RDM element normalized by the two preparation dipole lengths."""
        norm_i = numpy.linalg.norm(self.system.DD[ii, 0, :])
        norm_j = numpy.linalg.norm(self.system.DD[jj, 0, :])
        norm = norm_i * norm_j
        if numpy.abs(norm) <= tol:
            raise ValueError(
                "Cannot normalize density matrix element by zero transition "
                "dipole length"
            )
        return rdm[ii, jj] / norm

    def _four_dipole_prefactors(self, lab: Any) -> dict[str, Any]:
        """Return full orientational prefactors used by response functions."""
        prefactors = {}
        for key in ("abba", "baba"):
            prefactors[key] = numpy.einsum(
                "i,abi->ab", lab.F4eM4, self.system.get_F4d(key)
            )
        if self.system.mult > 1:
            for key in ("fbfaba", "fafbba"):
                prefactors[key] = numpy.einsum(
                    "i,fabi->fab", lab.F4eM4, self.system.get_F4d(key)
                )
        return prefactors

    def _response_backend_to_pump_probe(
        self, response: numpy.ndarray, tau: float
    ) -> Any:
        """Fourier transform a t1=0 response trace into a pump-probe spectrum."""
        onepp = PumpProbeSpectrum()
        onepp.set_axis(self.oa3)

        ppspec = -numpy.asarray(response, dtype=numpy.complex128)
        ft = numpy.fft.hfft(ppspec) * self.t3axis.step
        ft = numpy.fft.fftshift(ft)
        ft = numpy.flipud(ft)
        Nt = self.t3axis.length

        data = numpy.real(ft[Nt // 2 : Nt + Nt // 2])
        data = self.oa3.data * data

        onepp._add_data(data)
        onepp.set_t2(tau)
        return onepp

    def calculate_all_system_approx_response_backend(
        self,
        sys: Any,
        rdmt: Any = None,
        lab: Any = None,
        show_progress: bool = False,
        spec: Any = None,
        population_propagator: Any = None,
        density_matrix_propagator: Any = None,
        population_time_axis: Any = None,
        include_nonsecular_remainder: bool = True,
        include_remainder: bool = True,
        dipole_normalization_tol: float = 1.0e-12,
    ) -> Any:
        """Calculate pump-probe spectra through the modern response backend.

        This method is a pump-probe-like ``t1 = 0`` calculation.  The supplied
        dynamics can be a reduced-density-matrix trajectory, a one-exciton
        population propagator, or a one-exciton density-matrix propagator.
        """
        if spec is None:
            spec = ["Full"]
        if lab is None:
            raise ValueError("A lab setup has to be supplied")
        dynamics_count = sum(
            item is not None
            for item in (rdmt, population_propagator, density_matrix_propagator)
        )
        if dynamics_count != 1:
            raise ValueError(
                "Exactly one of rdmt, population_propagator, or "
                "density_matrix_propagator has to be supplied"
            )
        if not numpy.isclose(lab.F4eM4[1:], [0, 0], atol=1e-6).all():
            message = (
                "Lab is not set to the magic angle polarization which is"
                "the only supported measurement setting for this calculation. Set "
                "polarization angle between first two pulses and the last two to"
                "54.7356103 deg. and repeat the calculation."
            )
            raise OSError(message)

        self.system = copy.deepcopy(sys)
        if not self.system._diagonalized:
            self.system.diagonalize()

        t1axis = TimeAxis(0.0, 1, self.t3axis.step)
        backend_system = _RWAOverrideSystem(
            self.system, self.rwa, lineshape_timeaxis=[t1axis, self.t3axis]
        )
        response_types = []
        if "Full" in spec or "SE" in spec:
            response_types.extend(["R1g", "R2g"])
        if "Full" in spec or "GSB" in spec:
            response_types.extend(["R3g", "R4g"])
        if self.system.mult > 1 and ("Full" in spec or "ESA" in spec):
            response_types.extend(["R1f", "R2f"])

        if population_time_axis is None:
            population_time_axis = self.t2axis
        response_kwargs: dict[str, Any] = dict(
            population_time_axis=population_time_axis
        )
        if rdmt is not None:
            response_kwargs["density_matrix_trajectory"] = rdmt
            response_kwargs["dipole_normalization_tol"] = dipole_normalization_tol
        elif population_propagator is not None:
            response_kwargs["population_propagator"] = population_propagator
        else:
            response_kwargs["density_matrix_propagator"] = density_matrix_propagator
            response_kwargs["include_nonsecular_remainder"] = (
                include_nonsecular_remainder
            )

        if include_remainder and rdmt is None:
            if "Full" in spec or "SE" in spec:
                response_types.extend(["R1g_scM0g", "R2g_scM0g"])
            if self.system.mult > 1 and ("Full" in spec or "ESA" in spec):
                response_types.extend(
                    ["R1f_scM0g", "R2f_scM0g", "R1f_scM0e", "R2f_scM0e"]
                )

        responses = [
            NonLinearResponse(
                lab,
                backend_system,
                diagram,
                t1axis,
                self.t2axis,
                self.t3axis,
                **response_kwargs,
            )
            for diagram in response_types
        ]

        tcont = PumpProbeSpectrumContainer(t2axis=self.t2axis)

        kk = 0
        Nk = self.t2axis.length
        printProgressBar(
            kk, Nk, prefix="     - Progress:", suffix="Complete", length=50
        )

        for T2 in self.t2axis.data:
            if show_progress:
                print(" - calculating", kk, "of", Nk, "at t2 =", T2, "fs")

            response = numpy.zeros(self.t3axis.length, dtype=numpy.complex128)
            for rsp in responses:
                data = numpy.asarray(rsp.calculate_matrix(T2))
                response += numpy.ravel(data)

            ppspec1 = self._response_backend_to_pump_probe(response, T2)
            tcont.set_spectrum(ppspec1, tag=T2)

            kk += 1
            printProgressBar(
                kk, Nk, prefix="     - Progress:", suffix="Complete", length=50
            )

        return tcont

    def calculate(
        self,
        system: Any = None,
        lab: Any = None,
        density_matrix_trajectory: Any = None,
        population_propagator: Any = None,
        density_matrix_propagator: Any = None,
        population_time_axis: Any = None,
        method: str = "response",
        show_progress: bool = False,
        spec: Any = None,
        include_nonsecular_remainder: bool | None = None,
        include_remainder: bool | None = None,
        dipole_normalization_tol: float | None = None,
        approx: Any = None,
        gsb_mode: str = "hole",
        orientational_averaging: str = "legacy",
    ) -> Any:
        """Calculate pump-probe spectra with the selected backend.

        Parameters
        ----------
        system : Aggregate, optional
            System to calculate. If omitted, the calculator's stored system is
            used.
        lab : LabSetup, optional
            Laboratory setup. If omitted, the setup supplied to ``bootstrap`` is
            used.
        density_matrix_trajectory, population_propagator, density_matrix_propagator
            Optional external dynamics. If omitted, dynamics previously set by
            ``set_dynamics`` or one of the specialized setters are used.
        method : {"response", "legacy"}
            ``"response"`` uses the modern nonlinear-response backend.
            ``"legacy"`` uses the old RDM pump-probe implementation.
        """
        if system is None:
            system = self.system
        if lab is None:
            lab = self.lab
        if lab is None:
            raise ValueError("A lab setup has to be supplied or bootstrapped")

        if population_time_axis is None:
            population_time_axis = self.population_time_axis
        if density_matrix_trajectory is None:
            density_matrix_trajectory = self.density_matrix_trajectory
        if population_propagator is None:
            population_propagator = self.population_propagator
        if density_matrix_propagator is None:
            density_matrix_propagator = self.density_matrix_propagator
        if include_nonsecular_remainder is None:
            include_nonsecular_remainder = self.include_nonsecular_remainder
        if include_remainder is None:
            include_remainder = self.include_remainder
        if dipole_normalization_tol is None:
            dipole_normalization_tol = self.dipole_normalization_tol

        if method in ("response", "modern"):
            return self.calculate_all_system_approx_response_backend(
                system,
                rdmt=density_matrix_trajectory,
                lab=lab,
                show_progress=show_progress,
                spec=spec,
                population_propagator=population_propagator,
                density_matrix_propagator=density_matrix_propagator,
                population_time_axis=population_time_axis,
                include_nonsecular_remainder=include_nonsecular_remainder,
                include_remainder=include_remainder,
                dipole_normalization_tol=dipole_normalization_tol,
            )

        if method == "legacy":
            if density_matrix_trajectory is None:
                raise ValueError(
                    "Legacy pump-probe calculation requires density_matrix_trajectory"
                )
            if (
                population_propagator is not None
                or density_matrix_propagator is not None
            ):
                raise ValueError(
                    "Legacy pump-probe calculation does not accept propagators"
                )
            return self.calculate_all_system_approx(
                system,
                density_matrix_trajectory,
                lab,
                show_progress=show_progress,
                approx=approx,
                spec=spec,
                gsb_mode=gsb_mode,
                orientational_averaging=orientational_averaging,
            )

        raise ValueError("method has to be 'response' or 'legacy'")

    def calculate_all_system(
        self, sys: Any, eUt: Any, lab: Any, show_progress: bool = False
    ) -> Any:
        """Calculates all 2D spectra for a system and evolution superoperator"""
        tcont = PumpProbeSpectrumContainer(t2axis=self.t2axis)

        kk = 1
        Nk = self.t2axis.length

        printProgressBar(0, Nk, prefix="     - Progress:", suffix="Complete", length=50)
        for T2 in self.t2axis.data:
            if show_progress:
                print(" - calculating", kk, "of", Nk, "at t2 =", T2, "fs")

            ppspec1 = self.calculate_one_system(T2, sys, eUt, lab)

            tcont.set_spectrum(ppspec1, tag=T2)

            printProgressBar(
                kk, Nk, prefix="     - Progress:", suffix="Complete", length=50
            )
            kk += 1

        return tcont

    def calculate_one_system(
        self, t2: float, sys: Any, eUt: Any, lab: Any, pways: Any = None
    ) -> Any:
        """Returns pump-probe spectrum at t2 for a system and evolution
        superoperator

        """
        try:
            # print(t2)
            Uin = eUt.at(t2)
            # print(Uin.data.shape)
        except (AttributeError, IndexError):
            Uin = eUt
            # print("False")

        temp = sys.get_temperature()
        # FIXME: Delete zero temperature
        # temp = 0.0
        rho0 = sys.get_DensityMatrix(condition_type="thermal", temperature=temp)

        # if the Hamiltonian is larger than eUt, we will calculate ESA
        has_ESA = True
        H = self.system.get_Hamiltonian()

        # get Liouville pathways
        if has_ESA:
            pws = sys.liouville_pathways_3T(
                ptype=(
                    "R1g",
                    "R2g",  # SE
                    "R3g",
                    "R4g",  # GSB
                    "R1f*",
                    "R2f*",  # ESA
                    "R1f*E",
                    "R2f*E",
                ),  # Pathways with relaxation to the ground state
                # "R1gE", "R2gE"),
                eUt=eUt,
                ham=H,
                t2=t2,
                lab=lab,
            )
        else:
            pws = sys.liouville_pathways_3T(
                ptype=("R1g", "R2g", "R3g", "R4g", "R1f*E", "R2f*E"),
                eUt=eUt,
                ham=H,
                t2=t2,
                lab=lab,
            )

        self.set_pathways(pws)

        SS = self.system.SS.copy()
        self.goft_matrix = self._SE_excitonic_gofts(SS, self.system, tau=0.0)

        if pways is not None:
            pways[str(t2)] = pws

        pprobe1 = self.calculate_next(t2)

        return pprobe1

    def calculate_next(self, t2: float) -> Any:

        sone = self.calculate_one(self.tc, t2)
        # print(self.tc, sone)
        self.tc += 1
        return sone

    def calculate_one(self, tc: int, t2: float) -> Any:
        """Calculate the 2D spectrum for all pathways"""
        # import time

        onepp = PumpProbeSpectrum()
        onepp.set_axis(self.oa3)

        # start = time.time()
        data = self.calculate_pathways(self.pathways, t2)
        # end = time.time()
        # print("Single spectra calculation:", end - start)
        onepp._add_data(data)

        onepp.set_t2(self.t2axis.data[tc])

        return onepp

    def _c2g(self, timeaxis: Any, coft: numpy.ndarray) -> numpy.ndarray:
        """Converts correlation function to lineshape function

        Explicit numerical double integration of the correlation
        function to form a lineshape function.

        Parameters
        ----------
        timeaxis : cu.oqs.time.TimeAxis
            TimeAxis of the correlation function

        coft : complex numpy array
            Values of correlation function given at points specified
            in the TimeAxis object


        """
        ta = timeaxis
        rr = numpy.real(coft)
        ri = numpy.imag(coft)
        sr = scipy.interpolate.UnivariateSpline(ta.data, rr, s=0).antiderivative()(
            ta.data
        )
        sr = scipy.interpolate.UnivariateSpline(ta.data, sr, s=0).antiderivative()(
            ta.data
        )
        si = scipy.interpolate.UnivariateSpline(ta.data, ri, s=0).antiderivative()(
            ta.data
        )
        si = scipy.interpolate.UnivariateSpline(ta.data, si, s=0).antiderivative()(
            ta.data
        )
        gt = sr + 1j * si
        return gt

    def _excitonic_coft(self, SS: numpy.ndarray, AG: Any, n: int) -> numpy.ndarray:
        """Returns energy gap correlation function data of an exciton state n"""
        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC

        # get number of monomeric basis states
        Na = 0
        for monomer in AG.monomers:
            Na += monomer.nel - 1

        ct = numpy.zeros((self.t3axis.length), dtype=numpy.complex128)

        # electronic states corresponding to single excited states
        elst = numpy.where(AG.which_band == 1)[0]
        for el1 in elst:
            for el2 in elst:
                coft = DFunction(cfm.timeAxis, cfm.get_coft(el1, el2))
                ct3 = coft.at(self.t3axis.data)
                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        ct += (SS[kk, n] ** 2) * (SS[ll, n] ** 2) * ct3
        #        coft = DFunction(cfm.timeAxis,cfm.get_coft(2,2))
        #        ct3 = coft.at(self.t3axis.data)
        #        print(numpy.isclose(ct,ct3).all())
        return ct

    def _SE_excitonic_cofts(
        self, SS: numpy.ndarray, AG: Any, tau: float = 0
    ) -> tuple[numpy.ndarray, numpy.ndarray]:
        """Returns energy gap correlation function data of an exciton state n"""
        c0 = AG.monomers[0].get_egcf((0, 1))
        Nt = len(c0)

        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC
        ctimeAxis = cfm.timeAxis

        if self.t3axis.max + tau > ctimeAxis.max:
            raise OSError(
                "Correlation function should be defined on interval (0, t2_max+t3_max)."
            )

        # get number of monomeric basis states
        Na = 0
        for monomer in AG.monomers:
            Na += monomer.nel - 1

        ct3 = numpy.zeros(
            (AG.Ntot, AG.Ntot, self.t3axis.length), dtype=numpy.complex128
        )
        ct3tau = numpy.zeros(
            (AG.Ntot, AG.Ntot, self.t3axis.length), dtype=numpy.complex128
        )

        # electronic states corresponding to single excited states
        elst = numpy.where(AG.which_band == 1)[0]
        # print(elst)
        for el1 in elst:
            for el2 in elst:
                # get_coft starts from the excited state (ground not included in indexes)
                coft = DFunction(cfm.timeAxis, cfm.get_coft(el1 - 1, el2 - 1))
                coft_t3 = coft.at(self.t3axis.data)
                coft_t3tau = coft.at(self.t3axis.data + tau)
                # FIXME: both at time t3 and t3+tau

                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        for n in range(AG.Nb[0], AG.Nb[0] + AG.Nb[1]):
                            for m in range(n, AG.Nb[0] + AG.Nb[1]):
                                ct3[n, m] += (
                                    (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * coft_t3
                                )
                                ct3tau[n, m] += (
                                    (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * coft_t3tau
                                )
        for n in range(AG.Nb[0], AG.Nb[0] + AG.Nb[1]):
            for m in range(n + 1, AG.Nb[0] + AG.Nb[1]):
                ct3[m, n] += ct3[n, m]
                ct3tau[m, n] += ct3tau[n, m]

        # electronic states corresponding to double excited states
        elstd = numpy.where(AG.which_band == 2)[0]
        for el1 in elstd:
            for el2 in elstd:
                # print(el1,el2)
                coft = DFunction(cfm.timeAxis, cfm.get_coft(el1 - 1, el2 - 1))
                coft_t3 = coft.at(self.t3axis.data)
                coft_t3tau = coft.at(self.t3axis.data + tau)
                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        for n in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                            for m in range(n, numpy.sum(AG.Nb[0:3])):
                                ct3[n, m] += (
                                    (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * coft_t3
                                )
                                ct3tau[n, m] += (
                                    (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * coft_t3tau
                                )
        for n in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
            for m in range(n + 1, numpy.sum(AG.Nb[0:3])):
                ct3[m, n] += ct3[n, m]
                ct3tau[m, n] += ct3tau[n, m]

        # Mixed sigle-double excited state correlation functions
        for el1 in elst:
            for el2 in elstd:
                equal = numpy.where(AG.elsigs[el1] == AG.elsigs[el2])[0]
                if equal.size == 1:
                    coft = DFunction(cfm.timeAxis, cfm.get_coft(el1 - 1, el1 - 1))
                    coft_t3 = coft.at(self.t3axis.data)
                    coft_t3tau = coft.at(self.t3axis.data + tau)
                    for kk in AG.vibindices[el1]:
                        for ll in AG.vibindices[el2]:
                            for n in range(AG.Nb[0], AG.Nb[0] + AG.Nb[1]):
                                for m in range(
                                    numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])
                                ):
                                    ct3[n, m] += (
                                        (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * coft_t3
                                    )
                                    ct3tau[n, m] += (
                                        (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * coft_t3tau
                                    )
        for n in range(AG.Nb[0], AG.Nb[0] + AG.Nb[1]):
            for m in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                ct3[m, n] += ct3[n, m]
                ct3tau[m, n] += ct3tau[n, m]

        return ct3, ct3tau

    def _SE_excitonic_cofts_test(
        self, SS: numpy.ndarray, AG: Any, tau: float = 0
    ) -> tuple[numpy.ndarray, numpy.ndarray]:
        """Returns energy gap correlation function data of an exciton state n"""
        c0 = AG.monomers[0].get_egcf((0, 1))
        Nt = len(c0)

        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC
        ctimeAxis = cfm.timeAxis

        if self.t3axis.max + tau > ctimeAxis.max:
            raise OSError(
                "Correlation function should be defined on interval (0, t2_max+t3_max)."
            )

        # get number of monomeric basis states
        Na = 0
        for monomer in AG.monomers:
            Na += monomer.nel - 1

        ct = numpy.zeros(
            (AG.Ntot, AG.Ntot, cfm.timeAxis.length), dtype=numpy.complex128
        )
        gt3 = numpy.zeros(
            (AG.Ntot, AG.Ntot, self.t3axis.length), dtype=numpy.complex128
        )
        gt3tau = numpy.zeros(
            (AG.Ntot, AG.Ntot, self.t3axis.length), dtype=numpy.complex128
        )

        # electronic states corresponding to single excited states
        elst = numpy.where(AG.which_band == 1)[0]
        # print(elst)
        for el1 in elst:
            for el2 in elst:
                # get_coft starts from the excited state (ground not included in indexes)
                coft = cfm.get_coft(el1 - 1, el2 - 1)
                goft = DFunction(cfm.timeAxis, self._c2g(cfm.timeAxis, coft))
                goft_t3 = goft.at(self.t3axis.data)
                goft_t3tau = goft.at(self.t3axis.data + tau)
                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        for n in range(AG.Nb[0], AG.Nb[0] + AG.Nb[1]):
                            for m in range(n, AG.Nb[0] + AG.Nb[1]):
                                # ct[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*coft)
                                gt3[n, m] += (
                                    (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * goft_t3
                                )
                                gt3tau[n, m] += (
                                    (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * goft_t3tau
                                )
        for n in range(AG.Nb[0], AG.Nb[0] + AG.Nb[1]):
            for m in range(n + 1, AG.Nb[0] + AG.Nb[1]):
                # ct[m,n] += ct[n,m]
                gt3[m, n] += gt3[n, m]
                gt3tau[m, n] += gt3tau[n, m]

        # electronic states corresponding to double excited states
        elstd = numpy.where(AG.which_band == 2)[0]
        for el1 in elstd:
            for el2 in elstd:
                coft = cfm.get_coft(el1 - 1, el2 - 1)  # DFunction(cfm.timeAxis,)
                goft = DFunction(cfm.timeAxis, self._c2g(cfm.timeAxis, coft))
                goft_t3 = goft.at(self.t3axis.data)
                goft_t3tau = goft.at(self.t3axis.data + tau)
                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        for n in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                            for m in range(n, numpy.sum(AG.Nb[0:3])):
                                # ct[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*coft)
                                gt3[n, m] += (
                                    (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * goft_t3
                                )
                                gt3tau[n, m] += (
                                    (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * goft_t3tau
                                )
        for n in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
            for m in range(n + 1, numpy.sum(AG.Nb[0:3])):
                # ct[m,n] += ct[n,m]
                gt3[m, n] += gt3[n, m]
                gt3tau[m, n] += gt3tau[n, m]

        # Mixed sigle-double excited state correlation functions
        for el1 in elst:
            for el2 in elstd:
                equal = numpy.where(AG.elsigs[el1] == AG.elsigs[el2])[0]
                if equal.size == 1:
                    coft = cfm.get_coft(el1 - 1, el1 - 1)  # DFunction(cfm.timeAxis,)
                    goft = DFunction(cfm.timeAxis, self._c2g(cfm.timeAxis, coft))
                    goft_t3 = goft.at(self.t3axis.data)
                    goft_t3tau = goft.at(self.t3axis.data + tau)
                    for kk in AG.vibindices[el1]:
                        for ll in AG.vibindices[el2]:
                            for n in range(AG.Nb[0], AG.Nb[0] + AG.Nb[1]):
                                for m in range(
                                    numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])
                                ):
                                    # ct[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*coft)
                                    gt3[n, m] += (
                                        (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * goft_t3
                                    )
                                    gt3tau[n, m] += (
                                        (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * goft_t3tau
                                    )
        for n in range(AG.Nb[0], AG.Nb[0] + AG.Nb[1]):
            for m in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                # ct[m,n] += ct[n,m]
                gt3[m, n] += gt3[n, m]
                gt3tau[m, n] += gt3tau[n, m]

        if 0:
            # check if gts correctly defined
            for ii in range(AG.Nb[0], AG.Ntot):
                for jj in range(AG.Nb[0], AG.Ntot):
                    goft = DFunction(cfm.timeAxis, self._c2g(cfm.timeAxis, ct[ii, jj]))
                    if not numpy.isclose(gt3[ii, jj], goft.at(self.t3axis.data)).all():
                        print(ii, jj, False)
                        print(gt3[ii, jj] - goft.at(self.t3axis.data))

        return gt3, gt3tau

    def _SE_excitonic_gofts(
        self,
        SS: numpy.ndarray,
        AG: Any,
        tau: float = 0,
        _diag_double_only: bool = False,
    ) -> numpy.ndarray:
        """Returns energy gap correlation function data of an exciton state n"""
        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC
        ctimeAxis = cfm.timeAxis

        if self.t3axis.max + tau > ctimeAxis.max:
            raise OSError(
                "Correlation function should be defined on interval (0, t2_max+t3_max)."
            )

        # get number of monomeric basis states
        Na = 0
        for monomer in AG.monomers:
            Na += monomer.nel - 1

        gt3tau = numpy.zeros(
            (AG.Ntot, AG.Ntot, self.t3axis.length), dtype=numpy.complex128
        )

        # electronic states corresponding to single excited states
        elst = numpy.where(AG.which_band == 1)[0]
        # print(elst)
        for el1 in elst:
            for el2 in elst:
                # get_coft starts from the excited state (ground not included in indexes)
                coft = cfm.get_coft(el1 - 1, el2 - 1)
                goft = DFunction(cfm.timeAxis, self._c2g(cfm.timeAxis, coft))
                goft_t3tau = goft.at(self.t3axis.data + tau)
                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        for n in range(AG.Nb[0], AG.Nb[0] + AG.Nb[1]):
                            for m in range(n, AG.Nb[0] + AG.Nb[1]):
                                gt3tau[n, m] += (
                                    (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * goft_t3tau
                                )
        for n in range(AG.Nb[0], AG.Nb[0] + AG.Nb[1]):
            for m in range(n + 1, AG.Nb[0] + AG.Nb[1]):
                gt3tau[m, n] += gt3tau[n, m]

        # electronic states corresponding to double excited states
        elstd = numpy.where(AG.which_band == 2)[0]
        for el1 in elstd:
            for el2 in elstd:
                coft = cfm.get_coft(el1 - 1, el2 - 1)  # DFunction(cfm.timeAxis,)
                goft = DFunction(cfm.timeAxis, self._c2g(cfm.timeAxis, coft))
                goft_t3tau = goft.at(self.t3axis.data + tau)
                for kk in AG.vibindices[el1]:
                    for ll in AG.vibindices[el2]:
                        for n in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                            if _diag_double_only:
                                gt3tau[n, n] += (
                                    (SS[kk, n] ** 2) * (SS[ll, n] ** 2) * goft_t3tau
                                )
                            else:
                                for m in range(n, numpy.sum(AG.Nb[0:3])):
                                    gt3tau[n, m] += (
                                        (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * goft_t3tau
                                    )
        if not _diag_double_only:
            for n in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                for m in range(n + 1, numpy.sum(AG.Nb[0:3])):
                    gt3tau[m, n] += gt3tau[n, m]

        # Mixed sigle-double excited state correlation functions
        for el1 in elst:
            for el2 in elstd:
                # equal = numpy.where(AG.elsigs[el1] == AG.elsigs[el2])[0]
                # if equal.size == 1:
                # coft = cfm.get_coft(el1-1,el1-1) # DFunction(cfm.timeAxis,)
                # goft = DFunction(cfm.timeAxis,self._c2g(cfm.timeAxis,coft))
                # goft_t3tau = goft.at(self.t3axis.data + tau)
                # for kk in AG.vibindices[el1]:
                # for ll in AG.vibindices[el2]:
                # for n in range(AG.Nb[0], AG.Nb[0]+AG.Nb[1]):
                # for m in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                # gt3tau[n,m] += ((SS[kk,n]**2)*(SS[ll,m]**2)*goft_t3tau)
                if el1 in AG.twoex_indx[el2]:
                    if AG.twoex_indx[el2, 0] == el1:
                        el_exct = AG.twoex_indx[el2, 0]
                        coft = cfm.get_coft(
                            el_exct - 1, el_exct - 1
                        )  # DFunction(cfm.timeAxis,)
                    elif AG.twoex_indx[el2, 1] == el1:
                        el_exct = AG.twoex_indx[el2, 1]
                        coft = cfm.get_coft(
                            el_exct - 1, el_exct - 1
                        )  # DFunction(cfm.timeAxis,)
                    # else:
                    # coft = cfm.get_coft(el1-1,el1-1) # DFunction(cfm.timeAxis,)
                    goft = DFunction(cfm.timeAxis, self._c2g(cfm.timeAxis, coft))
                    goft_t3tau = goft.at(self.t3axis.data + tau)
                    for kk in AG.vibindices[el1]:
                        for ll in AG.vibindices[el2]:
                            for n in range(AG.Nb[0], AG.Nb[0] + AG.Nb[1]):
                                for m in range(
                                    numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])
                                ):
                                    gt3tau[n, m] += (
                                        (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * goft_t3tau
                                    )
        for n in range(AG.Nb[0], AG.Nb[0] + AG.Nb[1]):
            for m in range(numpy.sum(AG.Nb[0:2]), numpy.sum(AG.Nb[0:3])):
                gt3tau[m, n] += gt3tau[n, m]

        return gt3tau

    def _excitonic_reorg_energy(
        self, SS: numpy.ndarray, AG: Any
    ) -> tuple[numpy.ndarray, numpy.ndarray]:
        """Returns the reorganisation energy of an exciton state"""
        # SystemBathInteraction
        sbi = AG.get_SystemBathInteraction()
        # CorrelationFunctionMatrix
        cfm = sbi.CC

        reorg_exct = numpy.zeros(AG.Nb[0] + AG.Nb[1])
        reorg_exct_sd = numpy.zeros((AG.Nb[0] + AG.Nb[1], AG.Ntot))

        # electronic states corresponding to single excited states
        elst = numpy.where(AG.which_band == 1)[0]
        for n in range(reorg_exct.size):
            for el1 in elst:
                reorg = cfm.get_reorganization_energy(el1 - 1, el1 - 1)
                for kk in AG.vibindices[el1]:
                    reorg_exct[n] += (SS[kk, n] ** 2) * (SS[kk, n] ** 2) * reorg

        elst_sgl = numpy.where(AG.which_band == 1)[0]
        elst_dbl = numpy.where(AG.which_band == 2)[0]
        for n in range(AG.Nb[0] + AG.Nb[1]):
            for m in range(AG.Nb[0] + AG.Nb[1], AG.Ntot):
                for el1 in elst_sgl:
                    for el2 in elst_dbl:
                        reorg = cfm.get_reorganization_energy(el1 - 1, el2 - 1)
                        for kk in AG.vibindices[el1]:
                            for ll in AG.vibindices[el2]:
                                reorg_exct_sd[n, m] += (
                                    (SS[kk, n] ** 2) * (SS[ll, m] ** 2) * reorg
                                )

        return reorg_exct, reorg_exct_sd

    def calculate_pathways(self, pathways: Any, tau: float) -> numpy.ndarray:
        """Calculate the shape of a Liouville pathway"""
        # we can calculate empty pathway
        if pathways is None:
            N3 = self.oa3.length
            ppspec: numpy.ndarray = numpy.zeros(N3, dtype=COMPLEX)
            return ppspec

        N3 = self.oa3.length

        SS = self.system.SS.copy()

        # precalculate single excited state correlation functions
        #        start = time.time()
        #        ct3,ct3tau = self._SE_excitonic_cofts(SS,self.system, tau = tau)
        #        gt3s,gt3tau = self._SE_excitonic_cofts_test(SS,self.system,tau = tau)
        if self.goft_matrix is not None:
            gt3s = self.goft_matrix
        else:
            gt3s = self._SE_excitonic_gofts(
                SS, self.system, tau=0.0, _diag_double_only=True
            )
            self.goft_matrix = gt3s
        gt3tau = self._SE_excitonic_gofts(
            SS, self.system, tau=tau, _diag_double_only=True
        )
        #        end = time.time()
        #        print("Calculation of coft:", end - start)

        #        start = time.time()
        #        # convert correlation functions to lineshape functions
        #        ngs = self.system.Nb[0]
        #        nes = self.system.Nb[1]
        #        nfs = self.system.Nb[2]
        #        gt3s = numpy.zeros((ngs+nes+nfs,ngs+nes+nfs,self.t3axis.length),dtype=numpy.complex128)
        #        for ii in range(ngs, ngs + nes + nfs): #self.system.Ntot): # self.system.Nb[0]+self.system.Nb[1]):
        #            for jj in range(ii, ngs + nes + nfs):# self.system.Ntot): # self.system.Nb[0]+self.system.Nb[1]):
        #                gt3s[ii,jj,:] = self._c2g(self.t3axis,ct3[ii,jj])
        #                if ii!=jj:
        #                    gt3s[jj,ii,:] = gt3s[ii,jj,:]
        #
        #        end = time.time()
        #        print("Conversion coft 2 goft:", end - start)

        #        start = time.time()
        ppspec = numpy.zeros(self.t3axis.length, dtype=numpy.complex128)
        for pwy in self.pathways:
            _is_relax = False
            _is_jump = True
            if len(pwy.relaxations) == 1:
                _is_relax = True
                if pwy.relaxations[0][0] == pwy.relaxations[0][1]:
                    _is_jump = False

            om = pwy.frequency[-2] - self.rwa
            pref = pwy.pref
            if _is_relax:
                omtau = pwy.frequency[-3]
            else:
                omtau = pwy.frequency[1]

            if pwy.pathway_name not in ["R1f*", "R2f*"]:
                #                pass
                if (
                    _is_relax and _is_jump and self._separate_relax
                ):  # Pathways with relaxation
                    n = pwy.states[-2][0]
                    # print(om,n,pwy.states[2],pwy.relaxations[0],pwy.pathway_name)
                    ppspec += pref * numpy.exp(-gt3s[n, n] - 1j * om * self.t3axis.data)

                else:  # Pathways without relaxation
                    if pwy.pathway_name in ["R1g", "R2g"]:  # Stimulated emission
                        if _is_relax:
                            state = pwy.states[-3]
                        else:
                            state = pwy.states[1]
                        ft = -1j * om * self.t3axis.data - 1j * omtau * tau
                        #                        coft = (ct3tau[state[0],state[0]]
                        #                                 - numpy.conj(ct3tau[state[0],state[1]])
                        #                                 + numpy.conj(ct3[state[1],state[0]]))
                        #                        goft = self._c2g(self.t3axis,coft)
                        gt1 = gt3tau[state[0], state[0]] - numpy.conj(
                            gt3tau[state[0], state[1]]
                        )
                        gt2 = (
                            numpy.conj(gt3tau[state[1], state[1], 0])
                            - gt3tau[state[0], state[1], 0]
                        )
                        gt3 = numpy.conj(gt3s[state[1], state[0]])
                        ft -= gt1 + gt2 + gt3

                        #                        gt1 = self._c2g(self.t3axis,ct3tau[state[0],state[0]]
                        #                                    - numpy.conj(ct3tau[state[0],state[1]]))
                        #                        gt2 = self._c2g(self.t3axis,numpy.conj(ct3tau[state[1],state[1]])
                        #                                    - ct3tau[state[0],state[1]])
                        #                        gt3 = numpy.conj(gt3s[state[1],state[0]])
                        ##                        print(numpy.isclose(goft,gt1+gt3).all())
                        #                        ft -= gt1 + gt2[0] + gt3
                        ##
                        ##                        ft -= goft + gt2[0]
                        ppspec += pref * numpy.exp(ft)

                    elif pwy.pathway_name in ["R3g", "R4g"]:  # Ground state bleach
                        ft = -1j * om * self.t3axis.data - 1j * omtau * tau
                        state = pwy.states[-2]
                        n = max(state)
                        ft -= gt3s[n, n]
                        ppspec += pref * numpy.exp(ft)
            else:
                if (
                    _is_relax and _is_jump and self._separate_relax
                ):  # Pathways with relaxation
                    # Ecxited state absorption
                    n = pwy.states[-2][0]
                    m = pwy.states[-2][1]

                    ft = -1j * om * self.t3axis.data
                    ft -= gt3s[n, n] + numpy.conj(gt3s[m, m])
                    ft += gt3s[n, m] + numpy.conj(gt3s[m, n])
                    ppspec += pref * numpy.exp(ft)

                else:  # Pathways without relaxation
                    # Ecxited state absorption
                    Fl = pwy.states[-2][0]
                    Ek = pwy.states[-2][1]
                    if _is_relax:
                        Ej = pwy.states[-3][0]
                    else:
                        Ej = pwy.states[1][0]

                    #                    ######### TEST  ###########
                    #                    sbi = self.system.get_SystemBathInteraction()
                    #                    # CorrelationFunctionMatrix
                    #                    cfm = sbi.CC
                    #                    coft = cfm.get_coft(0,0)
                    #                    goft1 = DFunction(cfm.timeAxis,self._c2g(cfm.timeAxis,coft))
                    #                    coft = cfm.get_coft(1,1)
                    #                    goft2 = DFunction(cfm.timeAxis,self._c2g(cfm.timeAxis,coft))
                    #                    ###########################

                    ft = -1j * om * self.t3axis.data - 1j * omtau * tau

                    # Faster, but might cause numerical errors.
                    gt1 = gt3s[Fl, Fl] - gt3s[Fl, Ek] + gt3s[Ej, Ek] - gt3s[Ej, Fl]

                    gt2 = (
                        gt3tau[Ej, Ej, 0]
                        - numpy.conj(gt3tau[Ej, Ek, 0])
                        - gt3tau[Fl, Ej, 0]
                        + numpy.conj(gt3tau[Fl, Ek, 0])
                    )

                    gt3 = (
                        numpy.conj(gt3tau[Ek, Ek])
                        - gt3tau[Ek, Ej]
                        + gt3tau[Fl, Ej]
                        - numpy.conj(gt3tau[Fl, Ek])
                    )
                    ft -= gt1 + gt2 + gt3

                    #                    if tau == 0.0:
                    #                        print(Ej,Ek,Fl)
                    #
                    #                        ft2 = (goft1.data[:ft.size] - SS[1,2]**2*numpy.conj(goft1.data[:ft.size])) * SS[2,2]**2
                    #                        ft2 += (goft2.data[:ft.size] - SS[2,2]**2*numpy.conj(goft2.data[:ft.size])) * SS[1,2]**2
                    #                        print("real:",numpy.isclose(numpy.real(gt1 + gt2 + gt3),numpy.real(ft2)))
                    #                        print("imag:",numpy.isclose(numpy.imag(gt1 + gt2 + gt3),numpy.imag(ft2)))
                    #                        if Ej==2 and Ek==2:
                    #                            print("real:",numpy.real(gt1 + gt2 + gt3))
                    #                            print("real:",numpy.real(ft2))
                    #                            print("imag:",numpy.imag(gt1 + gt2 + gt3))
                    #                            print("imag:",numpy.imag(ft2))
                    #                            print(gt2)
                    #                            print(gt3tau.shape)
                    #                            print(goft2.data.size)
                    #                            print(ft.size)
                    #                            print(om)

                    ppspec += pref * numpy.exp(ft)

        #        end = time.time()
        #        print("Calculation of the pathways:", end - start)

        ppspec = -ppspec

        # Fourier transform the result
        ft = numpy.fft.hfft(ppspec) * self.t3axis.step
        ft = numpy.fft.fftshift(ft)
        # invert the order because hfft is a transform with -i
        ft = numpy.flipud(ft)
        # cut the center of the spectrum
        Nt = self.t3axis.length  # len(ta.data)

        data = numpy.real(ft[Nt // 2 : Nt + Nt // 2])

        data = self.oa3.data * data

        return data

    def calculate_pathways_rdm(
        self,
        rdm0: Any,
        rdm: Any,
        tau: float,
        lab: Any,
        ptol: float = 1.0e-6,
        spec: Any = None,
        gsb_mode: str = "hole",
        orientational_averaging: str = "legacy",
    ) -> Any:
        """Calculate the shape of a Liouville pathway
        so far implemented only for electronic
        aggregate.
        """
        if spec is None:
            spec = ["Full"]
        if gsb_mode not in ("hole", "response"):
            raise ValueError("gsb_mode has to be 'hole' or 'response'")
        if orientational_averaging not in ("legacy", "four_dipole"):
            raise ValueError(
                "orientational_averaging has to be 'legacy' or 'four_dipole'"
            )
        onepp = PumpProbeSpectrum()
        onepp.set_axis(self.oa3)

        SS = self.system.SS.copy()

        # precalculate single excited state correlation functions
        if self.goft_matrix is not None:
            gt3s = self.goft_matrix
        else:
            gt3s = self._SE_excitonic_gofts(
                SS, self.system, tau=0.0, _diag_double_only=True
            )
            self.goft_matrix = gt3s
        gt3tau = self._SE_excitonic_gofts(
            SS, self.system, tau=tau, _diag_double_only=True
        )

        # initialize the spectra
        ppspec = numpy.zeros(self.t3axis.length, dtype=numpy.complex128)

        # Compute the spectra
        # Make sure that the aggregate was diagonalized or we are working in
        # the eigenbais => self.system.DD[nf,ni,:] would be proper transition
        # dipoles between eigenstates.
        dim = self.system.Nb[1] + self.system.Nb[0]
        N0 = self.system.Nb[0]
        N1 = self.system.Nb[1]
        full_pref = None
        if orientational_averaging == "four_dipole":
            full_pref = self._four_dipole_prefactors(lab)

        # GSB (Ground state bleach)
        if "Full" in spec or "GSB" in spec:
            if gsb_mode == "response":
                ppspec += self._response_backend_trace(["R3g", "R4g"], tau, lab)
            else:
                for jj in range(1, dim):
                    pref_GSB = lab.F4eM4[0]
                    pref_GSB *= 2  # There are two pathways leading to the same results R3 and R4 (therefore twice)
                    pref_GSB *= numpy.sum(
                        numpy.diag(rdm0)
                    )  # The GSB signal is dependent only on the last state
                    # => prefactor can be computed as a sum before
                    pref_GSB *= numpy.dot(
                        self.system.DD[jj, 0, :], self.system.DD[jj, 0, :]
                    )

                    om = self.system.HH[jj, jj] - self.system.HH[0, 0] - self.rwa
                    omtau = 0  # during the t2 time bra and ket both in the ground state

                    ft = -1j * om * self.t3axis.data - 1j * omtau * tau

                    ft -= gt3s[jj, jj]
                    ppspec += pref_GSB * numpy.exp(ft)

        # SE
        if "Full" in spec or "SE" in spec:
            for ii in range(1, dim):
                om = self.system.HH[ii, ii] - self.system.HH[0, 0] - self.rwa

                for jj in range(1, dim):
                    if rdm[ii, jj] < ptol:
                        continue

                    state = [ii, jj]
                    omtau = self.system.HH[ii, ii] - self.system.HH[jj, jj]

                    if orientational_averaging == "four_dipole":
                        assert full_pref is not None
                        a = ii - N0
                        b = jj - N0
                        pref_SE = self._full_dipole_rdm_weight(rdm, ii, jj)
                        pref_SE *= full_pref["abba"][b, a]
                        pref_SE += (
                            self._full_dipole_rdm_weight(rdm, ii, jj)
                            * full_pref["baba"][b, a]
                        )
                    else:
                        pref_SE = lab.F4eM4[0]
                        pref_SE *= 2  # There are two pathways leading to the same results R1 and R2 (therefore twice)
                        pref_SE *= rdm[
                            ii, jj
                        ]  # It should include excitation weighted by evolution
                        pref_SE *= numpy.dot(
                            self.system.DD[ii, 0, :], self.system.DD[jj, 0, :]
                        )

                    ft = -1j * om * self.t3axis.data - 1j * omtau * tau

                    gt1 = gt3tau[state[0], state[0]] - numpy.conj(
                        gt3tau[state[0], state[1]]
                    )
                    gt2 = (
                        numpy.conj(gt3tau[state[1], state[1], 0])
                        - gt3tau[state[0], state[1], 0]
                    )
                    gt3 = numpy.conj(gt3s[state[1], state[0]])
                    ft -= gt1 + gt2 + gt3

                    ppspec += pref_SE * numpy.exp(ft)

        # ESA
        if "Full" in spec or "ESA" in spec:
            for ii in range(1, dim):
                for jj in range(1, dim):
                    if numpy.abs(rdm[ii, jj]) < ptol:
                        continue

                    for ll in range(dim, self.system.Ntot):
                        # Ek Ek      =   jj jj
                        # Fl Ek      =   ll jj
                        # Ei Ek      =   ii jj

                        om = self.system.HH[ll, ll] - self.system.HH[jj, jj] - self.rwa
                        omtau = self.system.HH[ii, ii] - self.system.HH[jj, jj]

                        if orientational_averaging == "four_dipole":
                            assert full_pref is not None
                            a = ii - N0
                            b = jj - N0
                            f = ll - N0 - N1
                            pref_ESA = self._full_dipole_rdm_weight(rdm, ii, jj)
                            pref_ESA *= full_pref["fbfaba"][f, b, a]
                            pref_ESA += (
                                self._full_dipole_rdm_weight(rdm, ii, jj)
                                * full_pref["fafbba"][f, b, a]
                            )
                        else:
                            pref_ESA = lab.F4eM4[0]
                            pref_ESA *= 2  # There are two pathways leading to the same results R1* and R2* (therefore twice)
                            pref_ESA *= rdm[
                                ii, jj
                            ]  # It should include excitation weighted by evolution
                            pref_ESA *= numpy.dot(
                                self.system.DD[ll, ii, :], self.system.DD[jj, ll, :]
                            )

                        Fl = ll
                        Ek = jj
                        Ej = ii

                        ft = -1j * om * self.t3axis.data - 1j * omtau * tau

                        gt1 = gt3s[Fl, Fl] - gt3s[Fl, Ek] + gt3s[Ej, Ek] - gt3s[Ej, Fl]

                        gt2 = (
                            gt3tau[Ej, Ej, 0]
                            - numpy.conj(gt3tau[Ej, Ek, 0])
                            - gt3tau[Fl, Ej, 0]
                            + numpy.conj(gt3tau[Fl, Ek, 0])
                        )

                        gt3 = (
                            numpy.conj(gt3tau[Ek, Ek])
                            - gt3tau[Ek, Ej]
                            + gt3tau[Fl, Ej]
                            - numpy.conj(gt3tau[Fl, Ek])
                        )
                        ft -= gt1 + gt2 + gt3

                        ppspec -= pref_ESA * numpy.exp(ft)

        ppspec = -ppspec

        # Fourier transform the result
        ft = numpy.fft.hfft(ppspec) * self.t3axis.step
        ft = numpy.fft.fftshift(ft)
        # invert the order because hfft is a transform with -i
        ft = numpy.flipud(ft)
        # cut the center of the spectrum
        Nt = self.t3axis.length  # len(ta.data)

        data = numpy.real(ft[Nt // 2 : Nt + Nt // 2])

        data = self.oa3.data * data

        onepp._add_data(data)
        onepp.set_t2(tau)

        return onepp

    def calculate_pathways_rdm_novoderezhkin(
        self,
        rdm0: Any,
        rdm: Any,
        tau: float,
        lab: Any,
        ptol: float = 1.0e-6,
        spec: Any = None,
        gsb_mode: str = "hole",
        orientational_averaging: str = "legacy",
    ) -> Any:
        """Calculate the shape of a Liouville pathway
        so far implemented only for electronic
        aggregate.
        """
        if spec is None:
            spec = ["Full"]
        if gsb_mode not in ("hole", "response"):
            raise ValueError("gsb_mode has to be 'hole' or 'response'")
        if orientational_averaging not in ("legacy", "four_dipole"):
            raise ValueError(
                "orientational_averaging has to be 'legacy' or 'four_dipole'"
            )
        onepp = PumpProbeSpectrum()
        onepp.set_axis(self.oa3)

        SS = self.system.SS.copy()

        # precalculate single excited state correlation functions
        if self.goft_matrix is not None:
            gt3s = self.goft_matrix
        else:
            gt3s = self._SE_excitonic_gofts(
                SS, self.system, tau=0.0, _diag_double_only=True
            )
            self.goft_matrix = gt3s

        if self.reorg_matrix is not None:
            reorg_exct = self.reorg_matrix[0]
            reorg_exct_sd = self.reorg_matrix[1]
        else:
            reorg_exct, reorg_exct_sd = self._excitonic_reorg_energy(SS, self.system)
            self.reorg_matrix = [reorg_exct, reorg_exct_sd]

        # initialize the spectra
        ppspec = numpy.zeros(self.t3axis.length, dtype=numpy.complex128)

        # Compute the spectra
        # Make sure that the aggregate was diagonalized or we are working in
        # the eigenbais => self.system.DD[nf,ni,:] would be proper transition
        # dipoles between eigenstates.
        dim = self.system.Nb[1] + self.system.Nb[0]
        N0 = self.system.Nb[0]
        N1 = self.system.Nb[1]
        full_pref = None
        if orientational_averaging == "four_dipole":
            full_pref = self._four_dipole_prefactors(lab)

        # GSB (Ground state bleach)
        if "Full" in spec or "GSB" in spec:
            if gsb_mode == "response":
                ppspec += self._response_backend_trace(["R3g", "R4g"], tau, lab)
            else:
                for jj in range(1, dim):
                    pref_GSB = lab.F4eM4[0]
                    pref_GSB *= 2  # There are two pathways leading to the same results R3 and R4 (therefore twice)
                    pref_GSB *= numpy.sum(
                        numpy.diag(rdm0)
                    )  # The GSB signal is dependent only on the last state
                    # => prefactor can be computed as a sum before
                    pref_GSB *= numpy.dot(
                        self.system.DD[jj, 0, :], self.system.DD[jj, 0, :]
                    )

                    om = self.system.HH[jj, jj] - self.system.HH[0, 0] - self.rwa

                    ft = -1j * om * self.t3axis.data

                    ft -= gt3s[jj, jj]
                    ppspec += pref_GSB * numpy.exp(ft)

        # SE
        if "Full" in spec or "SE" in spec:
            for ii in range(1, dim):
                om = self.system.HH[ii, ii] - self.system.HH[0, 0] - self.rwa
                state = [ii, ii]

                if orientational_averaging == "four_dipole":
                    assert full_pref is not None
                    a = ii - N0
                    weight = self._full_dipole_rdm_weight(rdm, ii, ii)
                    pref_SE = weight * (
                        full_pref["abba"][a, a] + full_pref["baba"][a, a]
                    )
                else:
                    pref_SE = lab.F4eM4[0]
                    pref_SE *= 2  # There are two pathways leading to the same results R1 and R2 (therefore twice)
                    pref_SE *= rdm[
                        ii, ii
                    ]  # It should include excitation weighted by evolution
                    pref_SE *= numpy.dot(
                        self.system.DD[ii, 0, :], self.system.DD[ii, 0, :]
                    )

                ft = (
                    -1j * om * self.t3axis.data
                    + 2 * 1j * reorg_exct[ii] * self.t3axis.data
                )
                ft -= numpy.conj(gt3s[state[0], state[0]])

                ppspec += pref_SE * numpy.exp(ft)

        # ESA
        if "Full" in spec or "ESA" in spec:
            for ii in range(1, dim):
                if numpy.abs(rdm[ii, ii]) < ptol:
                    continue
                for ll in range(dim, self.system.Ntot):
                    # Ek Ek      =   jj jj
                    # Fl Ek      =   ll jj
                    # Ei Ek      =   ii jj

                    om = self.system.HH[ll, ll] - self.system.HH[ii, ii] - self.rwa

                    if orientational_averaging == "four_dipole":
                        assert full_pref is not None
                        a = ii - N0
                        f = ll - N0 - N1
                        weight = self._full_dipole_rdm_weight(rdm, ii, ii)
                        pref_ESA = weight * (
                            full_pref["fbfaba"][f, a, a] + full_pref["fafbba"][f, a, a]
                        )
                    else:
                        pref_ESA = lab.F4eM4[0]
                        pref_ESA *= 2  # There are two pathways leading to the same results R1* and R2* (therefore twice)
                        pref_ESA *= rdm[
                            ii, ii
                        ]  # It should include excitation weighted by evolution
                        pref_ESA *= numpy.dot(
                            self.system.DD[ll, ii, :], self.system.DD[ii, ll, :]
                        )

                    Fl = ll
                    Ek = ii

                    ft = (
                        -1j * om * self.t3axis.data
                        + 2
                        * 1j
                        * (reorg_exct_sd[Ek, Fl] - reorg_exct[Ek])
                        * self.t3axis.data
                    )
                    ft -= gt3s[Ek, Ek] + gt3s[Fl, Fl] - 2 * gt3s[Fl, Ek]

                    ppspec -= pref_ESA * numpy.exp(ft)

        ppspec = -ppspec

        # Fourier transform the result
        ft = numpy.fft.hfft(ppspec) * self.t3axis.step
        ft = numpy.fft.fftshift(ft)
        # invert the order because hfft is a transform with -i
        ft = numpy.flipud(ft)
        # cut the center of the spectrum
        Nt = self.t3axis.length  # len(ta.data)

        data = numpy.real(ft[Nt // 2 : Nt + Nt // 2])

        data = self.oa3.data * data

        onepp._add_data(data)
        onepp.set_t2(tau)

        return onepp


def calculate_from_2D(twod: Any) -> PumpProbeSpectrum:
    """Calculates pump-probe spectrum from 2D spectrum

    Calculates pump-probe spectrum from 2D spectrum
    using projection theorem

    Parameters
    ----------
    twod: TwoDSpectrum
        2D spectrum from which pump-probe will be calculated

    """
    pp = PumpProbeSpectrum()

    # waiting time
    t2 = twod.get_t2()
    pp.set_t2(t2)

    # frequency step for integration
    xaxis = twod.xaxis
    dw = xaxis.step

    # real part of the total 2D spectrum or response
    if hasattr(twod, "d__data"):
        twod.set_data_flag(signal_TOTL)
        tddata = numpy.real(twod.d__data)
    else:
        tddata = numpy.real(twod.data)

    # integration over omega_1 axis with the detection-frequency prefactor
    ppdata = -twod.yaxis.data * numpy.sum(tddata, axis=1) * dw / (2.0 * numpy.pi)

    # setting pump-probe data
    pp.set_data(ppdata)
    pp.set_axis(twod.yaxis)

    return pp


class MockPumpProbeSpectrumCalculator(MockTwoDSpectrumCalculator):
    """Effective line-shape pump-probe spectrum calculator.

    Uses Gaussian or Lorentzian lineshapes evaluated directly in the
    frequency domain to produce pump-probe spectra without computing the
    full time-domain response.

    Parameters
    ----------
    t1axis : TimeAxis
        Coherence-time axis (used by parent class; not directly needed for
        pump-probe).
    t2axis : TimeAxis
        Waiting-time axis over which spectra are calculated.
    t3axis : TimeAxis
        Detection-time axis used to build the frequency axis.
    temp : float or None, optional
        Temperature in Kelvin. Default is ``None`` (temperature-independent
        lineshapes).
    """

    def calculate_all_system(
        self,
        sys: Any,
        eUt: Any,
        lab: Any,
        selection: Any = None,
        show_progress: bool = False,
        dtol: float = 1.0e-12,
        H: Any = None,
        **kwargs: Any,
    ) -> Any:
        """Calculate all pump-probe spectra for all t2 times in the t2axis.

        Parameters
        ----------
        sys : Aggregate
            Molecular aggregate system.
        eUt : evolution superoperator
            Evolution superoperator evaluated at each t2 time.
        lab : LabSetup
            Laboratory setup with polarization information.
        selection : optional
            Pathway selection; not currently used.
        show_progress : bool, optional
            If ``True``, print progress messages. Default is ``False``.
        dtol : float, optional
            Tolerance for dipole moment contributions. Default is ``1e-12``.
        H : optional
            Hamiltonian; not currently used.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        PumpProbeSpectrumContainer
            Container with all calculated pump-probe spectra.
        """
        temporary_fix = True

        if temporary_fix:
            # calculation via 2D spectrum
            tcont = super().calculate_all_system(sys, eUt, lab)
            return tcont.get_PumpProbeSpectrumContainer()

    def calculate_one_system(
        self,
        t2: float,
        sys: Any,
        eUt: Any,
        lab: Any,
        selection: Any = None,
        pways: Any = None,
        dtol: float = 1.0e-12,
        H: Any = None,
        **kwargs: Any,
    ) -> Any:
        """Calculate a single pump-probe spectrum at the specified t2 time.

        Parameters
        ----------
        t2 : float
            Waiting time at which to calculate the spectrum.
        sys : Aggregate
            Molecular aggregate system.
        eUt : evolution superoperator
            Evolution superoperator.
        lab : LabSetup
            Laboratory setup with polarization information.
        selection : optional
            Pathway selection; not currently used.
        pways : optional
            Pathway storage; not currently used.
        dtol : float, optional
            Tolerance for dipole moment contributions. Default is ``1e-12``.
        H : optional
            Hamiltonian; not currently used.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        PumpProbeSpectrum
            Calculated pump-probe spectrum at t2.
        """
        temporary_fix = True

        if temporary_fix:
            # calculation via 2D spectrum
            tcont = super().calculate_one_system(t2, sys, eUt, lab)
            return tcont.get_PumpProbeSpectrum()


def printProgressBar(
    iteration: int,
    total: int,
    prefix: str = "",
    suffix: str = "",
    decimals: int = 1,
    length: int = 100,
    fill: str = "█",
) -> None:  # █ = U-219
    """Call n a loop to create terminal progress bar

    Parameters
    --------
    iteration: int (Required)
        Current iteration
    total: int (Required)
        Total interactions
    prefix: string (optional)
        Prefix string
    suffix: string (optional)
        Suffix string
    decimal: int (optional)
        positive number of decimals in percent complete
    length: int (optional)
        Character length of bar
    fill: str (optional)
        Fill bar character
    """
    percent = ("{:0." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledlength = int(length * iteration // total)
    bar = fill * filledlength + "-" * (length - filledlength)
    print(f"\r{prefix} |{bar}| {percent}% {suffix}", end="\r")
    # print New line after completition
    if iteration == total:
        print()
