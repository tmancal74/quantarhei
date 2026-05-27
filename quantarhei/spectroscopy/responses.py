from __future__ import annotations

from typing import Any

import numpy

from .. import COMPLEX, REAL, signal_NONR, signal_REPH
from ..core.managers import Manager
from ..core.time import TimeAxis
from ..exceptions import QuantarheiError
from ..qm.propagators.poppropagator import PopulationPropagator
from .response_implementations import get_implementation

"""
    This packege contains two methodologies for calculating non-linear
    response. The one based on the LiouvillPathway class is now deprecated.

    The new methodology is based on NonLinearResponse class





"""


RESPONSE_DIAGRAMS = {
    # Ground-state bleach
    # Corresponds to the GSB calls in the old Aceto/twod22 workflow.
    "R3g": {
        "rtype": "R",
        "signal": signal_REPH,
        "process": "GSB",
        "storage_type": "R3g",
        "transfer": False,
    },
    "R4g": {
        "rtype": "NR",
        "signal": signal_NONR,
        "process": "GSB",
        "storage_type": "R4g",
        "transfer": False,
    },
    # Stimulated emission
    "R1g": {
        "rtype": "NR",
        "signal": signal_NONR,
        "process": "SE",
        "storage_type": "R1g",
        "transfer": False,
    },
    "R2g": {
        "rtype": "R",
        "signal": signal_REPH,
        "process": "SE",
        "storage_type": "R2g",
        "transfer": False,
    },
    # Excited-state absorption.  The modern implementation names these
    # diagrams R1f/R2f; older storage comments sometimes use R1fs/R2fs.
    "R1f": {
        "rtype": "NR",
        "signal": signal_NONR,
        "process": "ESA",
        "storage_type": "R1f",
        "transfer": False,
    },
    "R2f": {
        "rtype": "R",
        "signal": signal_REPH,
        "process": "ESA",
        "storage_type": "R2f",
        "transfer": False,
    },
    # Relaxation/transfer corrections to SE-like response terms.
    "R1g_scM0g": {
        "rtype": "NR",
        "signal": signal_NONR,
        "process": "SE",
        "storage_type": "R1g_scM0g",
        "transfer": True,
        "transfer_channel": "scM0g",
    },
    "R2g_scM0g": {
        "rtype": "R",
        "signal": signal_REPH,
        "process": "SE",
        "storage_type": "R2g_scM0g",
        "transfer": True,
        "transfer_channel": "scM0g",
    },
    # Relaxation/transfer corrections to ESA-like response terms.
    "R1f_scM0g": {
        "rtype": "NR",
        "signal": signal_NONR,
        "process": "ESA",
        "storage_type": "R1f_scM0g",
        "transfer": True,
        "transfer_channel": "scM0g",
    },
    "R2f_scM0g": {
        "rtype": "R",
        "signal": signal_REPH,
        "process": "ESA",
        "storage_type": "R2f_scM0g",
        "transfer": True,
        "transfer_channel": "scM0g",
    },
    "R1f_scM0e": {
        "rtype": "NR",
        "signal": signal_NONR,
        "process": "ESA",
        "storage_type": "R1f_scM0e",
        "transfer": True,
        "transfer_channel": "scM0e",
    },
    "R2f_scM0e": {
        "rtype": "R",
        "signal": signal_REPH,
        "process": "ESA",
        "storage_type": "R2f_scM0e",
        "transfer": True,
        "transfer_channel": "scM0e",
    },
    "R1g_scM1g": {
        "rtype": "NR",
        "signal": signal_NONR,
        "process": "SE",
        "storage_type": "R1g_scM1g",
        "transfer": True,
        "transfer_channel": "scM1g",
    },
    "R2g_scM1g": {
        "rtype": "R",
        "signal": signal_REPH,
        "process": "SE",
        "storage_type": "R2g_scM1g",
        "transfer": True,
        "transfer_channel": "scM1g",
    },
    "R1f_scM1g": {
        "rtype": "NR",
        "signal": signal_NONR,
        "process": "ESA",
        "storage_type": "R1f_scM1g",
        "transfer": True,
        "transfer_channel": "scM1g",
    },
    "R2f_scM1g": {
        "rtype": "R",
        "signal": signal_REPH,
        "process": "ESA",
        "storage_type": "R2f_scM1g",
        "transfer": True,
        "transfer_channel": "scM1g",
    },
    "R1f_scM1e": {
        "rtype": "NR",
        "signal": signal_NONR,
        "process": "ESA",
        "storage_type": "R1f_scM1e",
        "transfer": True,
        "transfer_channel": "scM1e",
    },
    "R2f_scM1e": {
        "rtype": "R",
        "signal": signal_REPH,
        "process": "ESA",
        "storage_type": "R2f_scM1e",
        "transfer": True,
        "transfer_channel": "scM1e",
    },
}


def get_response_diagram_info(diagram: str) -> dict[str, Any]:
    """Returns explicit metadata for a nonlinear response diagram."""
    try:
        return RESPONSE_DIAGRAMS[diagram].copy()
    except KeyError:
        raise Exception("Unknown response diagram: " + diagram)


def _rate_matrix_data(rate_matrix: Any) -> numpy.ndarray:
    """Returns numerical data from a rate matrix object or array."""
    if hasattr(rate_matrix, "data"):
        return numpy.asarray(rate_matrix.data)
    return numpy.asarray(rate_matrix)


def get_single_exciton_rate_matrix(system: Any, rate_matrix: Any) -> numpy.ndarray:
    """Returns the single-exciton block of a rate matrix.

    ``OpenSystem.get_RateMatrix`` returns rates indexed by the full Hamiltonian
    basis. Response implementations, however, use single-exciton indices
    shifted to start at zero. If a user already supplies a single-exciton rate
    matrix, it is returned unchanged.
    """
    data = _rate_matrix_data(rate_matrix)

    dim = system.Nb[1]
    if data.shape[-2:] == (dim, dim):
        return data

    band1 = system.get_band(1)
    if len(data.shape) == 2:
        return data[numpy.ix_(band1, band1)]
    if len(data.shape) == 3:
        return data[:, band1, :][:, :, band1]

    raise Exception("Rate matrix has to be an array of rank 2 or 3")


def get_common_time_axis(*axes: Any) -> Any:
    """Returns a common zero-start TimeAxis covering all submitted axes."""
    dt = min(axis.step for axis in axes)
    tmax = max(axis.max for axis in axes)
    length = int(numpy.rint(tmax / dt)) + 1
    return TimeAxis(0.0, length, dt)


def validate_2d_time_axes(t1s: Any, t2s: Any, t3s: Any) -> None:
    """Validates 2D time axes for relaxation-enabled calculations."""
    if not numpy.isclose(t1s.step, t3s.step):
        raise Exception("t1 and t3 axes must have the same time step")

    nstep = numpy.rint(t2s.step / t1s.step)
    if not numpy.isclose(nstep * t1s.step, t2s.step):
        raise Exception("t2 time step must be a multiple of t1/t3 time step")


def _axis_indices(axis: Any, base_axis: Any) -> numpy.ndarray:
    """Returns integer indices of ``axis`` values on ``base_axis``."""
    if not axis.is_subset_of(base_axis):
        raise Exception("TimeAxis is not a subset of the population time axis")

    indices = numpy.rint((axis.data - base_axis.start) / base_axis.step).astype(
        numpy.int64
    )
    times = base_axis.start + indices * base_axis.step
    if not numpy.allclose(times, axis.data):
        raise Exception("TimeAxis is not a subset of the population time axis")

    return indices


class NonLinearResponse:
    """Non-linear response function for a specific Feynman diagram type.

    Parameters
    ----------
    lab : LabSetup
        Laboratory setup holding pulse polarization information.
    system : OpenSystem
        The quantum system for which the response is calculated.
    diagram : str
        Feynman diagram type (e.g. ``"R1g"``, ``"R2g"``).
    t1s : TimeAxis
        Time axis for the first coherence period.
    t2s : TimeAxis
        Time axis for the waiting period.
    t3s : TimeAxis
        Time axis for the second coherence period.
    """

    def __init__(
        self,
        lab: Any,
        system: Any,
        diagram: str,
        t1s: Any,
        t2s: Any,
        t3s: Any,
        rate_matrix: Any = None,
        relaxation_theory: str | None = None,
        rate_matrix_time_dependent: bool = False,
        relaxation_cutoff_time: float | None = None,
        rate_matrix_options: dict[str, Any] | None = None,
        population_time_axis: Any = None,
        population_propagator: Any = None,
        density_matrix_propagator: Any = None,
        density_matrix_trajectory: Any = None,
        population_dynamics_mode: str | None = None,
        include_nonsecular_remainder: bool = True,
        dipole_normalization_tol: float = 1.0e-12,
        jump_time_graining: int = 1,
        jump_kernel_cutoff: float = 0.0,
        jump_kernel_zero_cutoff: float = 0.0,
    ) -> None:

        # info about pulse polarizations
        self.lab = lab

        # info about energies, dipolemoments and rwa
        self.sys = system

        # which response to calculate; the function to calculate the respose
        self.diag = diagram
        self.diagram_info = get_response_diagram_info(self.diag)
        self.rtype = self.diagram_info["rtype"]
        self.signal = self.diagram_info["signal"]
        self.process = self.diagram_info["process"]
        self.storage_type = self.diagram_info["storage_type"]
        self.is_transfer = self.diagram_info["transfer"]
        self.transfer_channel = self.diagram_info.get("transfer_channel", None)

        self.func = get_implementation(self.diag)

        # what times to calculate for
        self.t1s = t1s
        self.t2s = t2s
        self.t3s = t3s
        if population_dynamics_mode is None:
            population_dynamics_mode = (
                "full_conditional"
                if population_propagator is not None
                or density_matrix_propagator is not None
                or density_matrix_trajectory is not None
                else "jump_decomposition"
            )

        self.population_time_axis = population_time_axis
        self.population_dynamics_mode = population_dynamics_mode
        self.include_nonsecular_remainder = include_nonsecular_remainder
        self.dipole_normalization_tol = dipole_normalization_tol
        self.jump_time_graining = jump_time_graining
        self.jump_kernel_cutoff = jump_kernel_cutoff
        self.jump_kernel_zero_cutoff = jump_kernel_zero_cutoff
        self.KK: numpy.ndarray
        self.U0_t1: numpy.ndarray
        self.U0_t2: numpy.ndarray
        self.U0_t3: numpy.ndarray
        self.U0_population: numpy.ndarray
        self.U0fe_t3: numpy.ndarray | None = None
        self.Uee: numpy.ndarray
        self.U1_t2: numpy.ndarray
        self.Usingle_t2: numpy.ndarray
        self.Ujump_t2: tuple[numpy.ndarray, ...]
        self.Uremainder_t2: numpy.ndarray
        self.Utransfer_t2: numpy.ndarray
        self.diagnostics: dict[str, Any] = {}
        self._previous_response_matrix: numpy.ndarray | None = None

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
        if external_count > 0 and (
            rate_matrix is not None or relaxation_theory is not None
        ):
            raise ValueError(
                "External propagators cannot be combined with rate_matrix "
                "or relaxation_theory"
            )

        if density_matrix_propagator is not None:
            self.set_density_matrix_propagator(
                density_matrix_propagator,
                population_time_axis=population_time_axis,
                mode=population_dynamics_mode,
                include_nonsecular_remainder=include_nonsecular_remainder,
            )
            return

        if density_matrix_trajectory is not None:
            self.set_density_matrix_trajectory(
                density_matrix_trajectory,
                population_time_axis=population_time_axis,
                mode=population_dynamics_mode,
                dipole_normalization_tol=dipole_normalization_tol,
            )
            return

        if population_propagator is not None:
            self.set_population_propagator(
                population_propagator,
                population_time_axis=population_time_axis,
                mode=population_dynamics_mode,
            )
            return

        if rate_matrix is not None:
            KK = get_single_exciton_rate_matrix(system, rate_matrix)
        elif relaxation_theory is not None:
            options = {} if rate_matrix_options is None else rate_matrix_options
            KK = system.get_RateMatrix(
                relaxation_theory=relaxation_theory,
                time_dependent=rate_matrix_time_dependent,
                relaxation_cutoff_time=relaxation_cutoff_time,
                **options,
            )
            KK = get_single_exciton_rate_matrix(system, KK)
        else:
            KK = None

        if KK is None:
            self._set_identity_dynamics(system.Nb[1])
        else:
            self.set_rate_matrix(KK, population_time_axis=population_time_axis)

    def calculate_matrix(self, t2: float) -> Any:
        """Calculate the matrix of response values over t1 and t3 times.

        Parameters
        ----------
        t2 : float
            Waiting time for which the response is calculated.

        Returns
        -------
        numpy.ndarray
            2D complex array of response values with shape
            ``(len(t1s), len(t3s))``.
        """
        # identify the index of the present t2 time
        out = self.t2s.locate(t2)
        t2i = out[0]
        pt2i = t2i
        pzeroi = 0
        if self.population_time_axis is not None:
            try:
                pzeroi = self.population_time_axis.locate(0.0)[0]
            except Exception:
                pzeroi = 0
            pt2i = self.population_time_axis.locate(t2)[0]

        # population decay factors at t2
        U0t2 = self.U0_t2[:, t2i]

        metadata: dict[str, Any] = {
            "jump_time_axis": self.population_time_axis,
            "jump_time_zero_index": pzeroi,
            "jump_time_t2_index": pt2i,
            "jump_time_graining": self.jump_time_graining,
            "jump_zero_propagator": self.U0_population,
            "jump_kernel_cutoff": self.jump_kernel_cutoff,
            "jump_kernel_zero_cutoff": self.jump_kernel_zero_cutoff,
        }
        if hasattr(self, "Udm_zero_t2"):
            metadata["density_matrix_endpoint_t2"] = self.Udm_zero_t2[:, :, t2i]

        # here we specify evolution matrices
        evol = (
            self.U0_t1,
            self.U0_t3,
            U0t2,
            self._get_transfer_matrix(t2i),
            metadata,
        )

        data = self.func(
            t2,
            self.t1s.data,
            self.t3s.data,
            self.lab,
            self.sys,
            evol,
            self.KK,
        )
        self._update_diagnostics(t2, data, metadata)
        return data

    def _update_diagnostics(
        self, t2: float, data: numpy.ndarray, metadata: dict[str, Any]
    ) -> None:
        """Update response diagnostics after a t2 calculation."""
        norm = float(numpy.linalg.norm(data))
        previous_norm = 0.0
        relative_change = None
        if self._previous_response_matrix is not None:
            previous_norm = float(numpy.linalg.norm(self._previous_response_matrix))
            denom = max(norm, previous_norm, 1.0e-300)
            relative_change = float(
                numpy.linalg.norm(data - self._previous_response_matrix) / denom
            )

        self.diagnostics = {
            "diagram": self.diag,
            "t2": float(t2),
            "response_norm": norm,
            "previous_response_norm": previous_norm,
            "relative_change": relative_change,
            "is_transfer": self.is_transfer,
            "transfer_channel": self.transfer_channel,
        }
        if self.transfer_channel is not None and self.transfer_channel.startswith(
            "scM0"
        ):
            self.diagnostics["remainder_relative_change"] = relative_change
        if "jump_diagnostics" in metadata:
            self.diagnostics["jump"] = metadata["jump_diagnostics"].get(self.diag, {})

        self._previous_response_matrix = data.copy()

    def set_rwa(self, rwa: float) -> None:
        """Sets rotating wave approximation frequency"""
        pass  # rwa is set through the system class, at least for now

    def _set_identity_dynamics(self, dim: int) -> None:
        """Sets relaxation-free population and coherence dynamics."""
        self.KK = numpy.zeros((dim, dim), dtype=REAL)
        self.U0_t1 = numpy.ones((dim, self.t1s.length), dtype=REAL)
        self.U0_t2 = numpy.ones((dim, self.t2s.length), dtype=REAL)
        self.U0_t3 = numpy.ones((dim, self.t3s.length), dtype=REAL)
        self.U0_population = numpy.ones((dim, self.t2s.length), dtype=REAL)

        self.Uee = numpy.zeros((dim, dim, self.t2s.length), dtype=REAL)
        for ii in range(self.t2s.length):
            self.Uee[:, :, ii] = numpy.eye(dim, dtype=REAL)
        self.U1_t2 = self.Uee.copy()
        self.Usingle_t2 = numpy.zeros_like(self.Uee)
        self.Ujump_t2 = (self.U1_t2,)
        self.Uremainder_t2 = self._get_remainder_propagator(self.Ujump_t2)
        self.Utransfer_t2 = self.Uremainder_t2

    def _uses_single_jump_storage(self) -> bool:
        """Returns ``True`` when line-shape storage is one-jump-ready."""
        if not hasattr(self.sys, "get_lineshape_functions"):
            return False
        try:
            gg = self.sys.get_lineshape_functions()
        except Exception:
            return False
        return hasattr(gg, "time_mapping") and "s2" in gg.time_mapping

    def _get_jump_expansion_order(self) -> int:
        """Returns the highest explicit jump order to calculate."""
        if self.population_dynamics_mode != "jump_decomposition":
            return 0
        return 1 if self._uses_single_jump_storage() else 0

    def _get_transfer_matrix(self, t2i: int) -> numpy.ndarray:
        """Returns transfer propagator for the current response channel."""
        if self.transfer_channel is not None and self.transfer_channel.startswith(
            "scM1"
        ):
            return self.Usingle_t2[:, :, t2i]
        return self.Uremainder_t2[:, :, t2i]

    def _get_remainder_propagator(
        self, jumps: tuple[numpy.ndarray, ...]
    ) -> numpy.ndarray:
        """Returns propagation not covered by explicit jump contributions."""
        jump_sum = numpy.zeros_like(self.Uee)
        for jump in jumps:
            jump_sum += jump
        return self.Uee - jump_sum

    def _population_propagator_data(self, population_propagator: Any) -> numpy.ndarray:
        """Returns population propagator data in ``(N, N, Nt)`` order."""
        data = numpy.asarray(
            population_propagator.data
            if hasattr(population_propagator, "data")
            else population_propagator
        )
        if len(data.shape) != 3:
            raise Exception("Population propagator has to be an array of rank 3")

        dim = self.sys.Nb[1]
        if data.shape[0] == dim and data.shape[1] == dim:
            return data.astype(REAL, copy=False)
        if data.shape[1] == dim and data.shape[2] == dim:
            return numpy.transpose(data, (1, 2, 0)).astype(REAL, copy=False)

        raise Exception(
            "Population propagator has to have shape (N, N, Nt) "
            "or (Nt, N, N), where N is the number of single-exciton states"
        )

    def _density_matrix_propagator_data(
        self, density_matrix_propagator: Any
    ) -> numpy.ndarray:
        """Returns density-matrix propagator data in ``(N, N, N, N, Nt)`` order."""
        data = numpy.asarray(
            density_matrix_propagator.data
            if hasattr(density_matrix_propagator, "data")
            else density_matrix_propagator
        )
        if len(data.shape) != 5:
            raise Exception("Density-matrix propagator has to be an array of rank 5")

        dim = self.sys.Nb[1]
        if data.shape[:4] == (dim, dim, dim, dim):
            return data.astype(COMPLEX, copy=False)
        if data.shape[1:] == (dim, dim, dim, dim):
            return numpy.transpose(data, (1, 2, 3, 4, 0)).astype(COMPLEX, copy=False)

        raise Exception(
            "Density-matrix propagator has to have shape (N, N, N, N, Nt) "
            "or (Nt, N, N, N, N), where N is the number of single-exciton states"
        )

    def _density_matrix_trajectory_data(
        self, density_matrix_trajectory: Any
    ) -> numpy.ndarray:
        """Returns density-matrix trajectory data in ``(N, N, Nt)`` order."""
        data = numpy.asarray(
            density_matrix_trajectory.data
            if hasattr(density_matrix_trajectory, "data")
            else density_matrix_trajectory
        )
        if len(data.shape) != 3:
            raise Exception("Density-matrix trajectory has to be an array of rank 3")

        dim = self.sys.Nb[1]
        if data.shape[0] == dim and data.shape[1] == dim:
            return data.astype(COMPLEX, copy=False)
        if data.shape[1] == dim and data.shape[2] == dim:
            return numpy.transpose(data, (1, 2, 0)).astype(COMPLEX, copy=False)

        band1 = self.sys.get_band(1)
        if data.shape[0] == self.sys.Ntot and data.shape[1] == self.sys.Ntot:
            return data[numpy.ix_(band1, band1)].astype(COMPLEX, copy=False)
        if data.shape[1] == self.sys.Ntot and data.shape[2] == self.sys.Ntot:
            subset = data[:, band1, :][:, :, band1]
            return numpy.transpose(subset, (1, 2, 0)).astype(COMPLEX, copy=False)

        raise Exception(
            "Density-matrix trajectory has to have shape (N, N, Nt), "
            "(Nt, N, N), or the corresponding full-system shape"
        )

    def _ground_exciton_dipole_lengths(self) -> numpy.ndarray:
        """Returns lengths of ground-to-one-exciton transition dipoles."""
        self.sys.diagonalize()
        band0 = self.sys.get_band(0)
        if len(band0) != 1:
            raise Exception(
                "Density-matrix trajectory normalization requires one ground state"
            )
        ground = band0[0]
        band1 = self.sys.get_band(1)
        lengths = numpy.zeros(self.sys.Nb[1], dtype=REAL)
        for aa, state in enumerate(band1):
            lengths[aa] = numpy.linalg.norm(self.sys.DD[state, ground, :])
        return lengths

    def set_density_matrix_trajectory(
        self,
        density_matrix_trajectory: Any,
        population_time_axis: Any = None,
        mode: str = "full_conditional",
        dipole_normalization_tol: float = 1.0e-12,
    ) -> None:
        """Set an externally propagated RDM trajectory as endpoint weights.

        The trajectory is intended for pump-probe-like calculations with
        ``t1 = 0``.  It is divided by the lengths of the two transition dipoles
        that prepared the initial excited-state density matrix.
        """
        if mode != "full_conditional":
            raise ValueError(
                "Density-matrix trajectories currently support only "
                "population_dynamics_mode='full_conditional'"
            )
        if self.t1s.length != 1 or not numpy.isclose(self.t1s.data[0], 0.0):
            raise ValueError(
                "Density-matrix trajectory weighting is only allowed for "
                "pump-probe-like calculations with t1 = 0"
            )

        if population_time_axis is None:
            population_time_axis = self.population_time_axis
        if population_time_axis is None:
            population_time_axis = self.t2s
        self.population_time_axis = population_time_axis
        self.population_dynamics_mode = mode
        self.dipole_normalization_tol = dipole_normalization_tol

        if not self.t2s.is_subset_of(population_time_axis):
            raise Exception("t2 TimeAxis is not compatible with population dynamics")

        rho = self._density_matrix_trajectory_data(density_matrix_trajectory)
        dim = rho.shape[0]
        if rho.shape[2] != population_time_axis.length:
            raise Exception(
                "Density-matrix trajectory time dimension has to match the "
                "population TimeAxis length"
            )

        dipoles = self._ground_exciton_dipole_lengths()
        norm = dipoles[:, None] * dipoles[None, :]
        if numpy.any(numpy.abs(norm) <= dipole_normalization_tol):
            raise ValueError(
                "Cannot normalize density matrix trajectory by zero transition "
                "dipole length"
            )

        weights = rho / norm[:, :, None]
        it2 = _axis_indices(self.t2s, population_time_axis)

        self.KK = numpy.zeros((dim, dim), dtype=REAL)
        self.U0_t1 = numpy.ones((dim, self.t1s.length), dtype=REAL)
        self.U0_t3 = numpy.ones((dim, self.t3s.length), dtype=REAL)
        self.U0_t2 = numpy.ones((dim, self.t2s.length), dtype=REAL)
        self.U0_population = numpy.ones((dim, population_time_axis.length), dtype=REAL)

        self.Udm_zero_t2 = weights[:, :, it2]
        self.Uee = numpy.zeros((dim, dim, self.t2s.length), dtype=COMPLEX)
        self.U1_t2 = numpy.zeros_like(self.Uee)
        for ii in range(self.t2s.length):
            self.U1_t2[:, :, ii] = numpy.diag(numpy.diag(self.Udm_zero_t2[:, :, ii]))
        self.Usingle_t2 = numpy.zeros_like(self.Uee)
        self.Ujump_t2 = (self.U1_t2,)
        self.Uremainder_t2 = numpy.zeros_like(self.Uee)
        self.Utransfer_t2 = self.Uremainder_t2

    def set_density_matrix_propagator(
        self,
        density_matrix_propagator: Any,
        population_time_axis: Any = None,
        mode: str = "full_conditional",
        include_nonsecular_remainder: bool = True,
    ) -> None:
        """Set an externally calculated one-exciton RDM superoperator.

        The endpoint-preserving elements ``U[a,b,a,b,t]`` weight ordinary
        response expressions.  Endpoint-changing elements are collapsed into
        an endpoint remainder when ``include_nonsecular_remainder`` is true.
        """
        if mode != "full_conditional":
            raise ValueError(
                "External density-matrix propagators currently support only "
                "population_dynamics_mode='full_conditional'"
            )

        if population_time_axis is None:
            population_time_axis = self.population_time_axis
        if population_time_axis is None:
            population_time_axis = self.t2s
        self.population_time_axis = population_time_axis
        self.population_dynamics_mode = mode
        self.include_nonsecular_remainder = include_nonsecular_remainder

        if not self.t2s.is_subset_of(population_time_axis):
            raise Exception("t2 TimeAxis is not compatible with population dynamics")

        Ufull = self._density_matrix_propagator_data(density_matrix_propagator)
        dim = Ufull.shape[0]
        if Ufull.shape[4] != population_time_axis.length:
            raise Exception(
                "Density-matrix propagator time dimension has to match the "
                "population TimeAxis length"
            )

        it2 = _axis_indices(self.t2s, population_time_axis)

        self.KK = numpy.zeros((dim, dim), dtype=REAL)
        self.U0_t1 = numpy.ones((dim, self.t1s.length), dtype=REAL)
        self.U0_t3 = numpy.ones((dim, self.t3s.length), dtype=REAL)
        self.U0_t2 = numpy.ones((dim, self.t2s.length), dtype=REAL)
        self.U0_population = numpy.ones((dim, population_time_axis.length), dtype=REAL)

        Uzero = numpy.zeros((dim, dim, population_time_axis.length), dtype=COMPLEX)
        Uremainder = numpy.zeros_like(Uzero)
        for aa in range(dim):
            for bb in range(dim):
                Uzero[aa, bb, :] = Ufull[aa, bb, aa, bb, :]
                if include_nonsecular_remainder:
                    Uremainder[aa, bb, :] = numpy.sum(
                        Ufull[aa, bb, :, :, :], axis=(0, 1)
                    )
                    Uremainder[aa, bb, :] -= Uzero[aa, bb, :]

        self.Udm_zero_t2 = Uzero[:, :, it2]
        self.Uee = numpy.zeros((dim, dim, self.t2s.length), dtype=COMPLEX)
        for ii, tindex in enumerate(it2):
            self.Uee[:, :, ii] = Ufull[:, :, :, :, tindex].trace(axis1=0, axis2=2)
        self.U1_t2 = numpy.zeros_like(self.Uee)
        for ii in range(self.t2s.length):
            self.U1_t2[:, :, ii] = numpy.diag(numpy.diag(self.Udm_zero_t2[:, :, ii]))

        self.U0_t2[:, :] = 1.0
        self.Usingle_t2 = numpy.zeros_like(self.Uee)
        self.Ujump_t2 = (self.U1_t2,)
        self.Uremainder_t2 = Uremainder[:, :, it2]
        self.Utransfer_t2 = self.Uremainder_t2

    def set_population_propagator(
        self,
        population_propagator: Any,
        population_time_axis: Any = None,
        mode: str = "full_conditional",
    ) -> None:
        """Set an externally calculated one-exciton population propagator.

        The ``full_conditional`` mode classifies the supplied conditional
        propagator by endpoints: diagonal elements weight the ordinary
        no-transfer response expressions, and off-diagonal elements weight the
        transfer/remainder response expressions.  No jump-order information is
        inferred from this propagator.
        """
        if mode != "full_conditional":
            raise ValueError(
                "External population propagators currently support only "
                "population_dynamics_mode='full_conditional'"
            )

        if population_time_axis is None:
            population_time_axis = self.population_time_axis
        if population_time_axis is None:
            population_time_axis = self.t2s
        self.population_time_axis = population_time_axis
        self.population_dynamics_mode = mode

        if not self.t2s.is_subset_of(population_time_axis):
            raise Exception("t2 TimeAxis is not compatible with population dynamics")

        Ufull = self._population_propagator_data(population_propagator)
        dim = Ufull.shape[0]
        if Ufull.shape[2] != population_time_axis.length:
            raise Exception(
                "Population propagator time dimension has to match the "
                "population TimeAxis length"
            )

        it2 = _axis_indices(self.t2s, population_time_axis)

        self.KK = numpy.zeros((dim, dim), dtype=REAL)
        self.U0_t1 = numpy.ones((dim, self.t1s.length), dtype=REAL)
        self.U0_t3 = numpy.ones((dim, self.t3s.length), dtype=REAL)
        self.U0_t2 = numpy.ones((dim, self.t2s.length), dtype=REAL)
        self.U0_population = numpy.ones((dim, population_time_axis.length), dtype=REAL)

        Udiag_full = numpy.zeros_like(Ufull)
        for ii in range(population_time_axis.length):
            diag = numpy.diag(Ufull[:, :, ii])
            if numpy.any(diag < -1.0e-12):
                raise Exception("Population propagator diagonal cannot be negative")
            Udiag_full[:, :, ii] = numpy.diag(diag)

        self.Uee = Ufull[:, :, it2]
        self.U1_t2 = Udiag_full[:, :, it2]
        for ii in range(self.t2s.length):
            diag = numpy.maximum(numpy.diag(self.U1_t2[:, :, ii]), 0.0)
            self.U0_t2[:, ii] = numpy.sqrt(diag)

        self.Usingle_t2 = numpy.zeros_like(self.Uee)
        self.Ujump_t2 = (self.U1_t2,)
        self.Uremainder_t2 = self.Uee - self.U1_t2
        self.Utransfer_t2 = self.Uremainder_t2

    def set_rate_matrix(
        self, KK: numpy.ndarray, population_time_axis: Any = None
    ) -> None:
        """Set the rate matrix and pre-compute the evolution coefficients.

        Parameters
        ----------
        KK : numpy.ndarray
            Rate matrix of shape ``(N, N)`` or time-dependent rate matrix of
            shape ``(Nt, N, N)``, where ``N`` is the number of single-excited
            states.

        Raises
        ------
        Exception
            If ``KK`` is not a square matrix or if any depopulation rate is
            positive.
        """
        KK = numpy.asarray(KK)

        if population_time_axis is None:
            population_time_axis = self.population_time_axis
        if population_time_axis is None:
            population_time_axis = get_common_time_axis(self.t1s, self.t2s, self.t3s)
        self.population_time_axis = population_time_axis

        if not self.t1s.is_subset_of(population_time_axis):
            raise Exception("t1 TimeAxis is not compatible with population dynamics")
        if not self.t2s.is_subset_of(population_time_axis):
            raise Exception("t2 TimeAxis is not compatible with population dynamics")
        if not self.t3s.is_subset_of(population_time_axis):
            raise Exception("t3 TimeAxis is not compatible with population dynamics")

        if len(KK.shape) == 2:
            dim = KK.shape[0]
            if KK.shape[0] != KK.shape[1]:
                raise Exception("Square matrix must be submitted")
        elif len(KK.shape) == 3:
            dim = KK.shape[1]
            if KK.shape[1] != KK.shape[2]:
                raise Exception("Square matrix must be submitted")
            if KK.shape[0] != population_time_axis.length:
                raise Exception(
                    "Time-dependent rate matrix has to have the same length "
                    "as the population TimeAxis"
                )
        else:
            raise Exception("Rate matrix has to be an array of rank 2 or 3")

        self.KK = KK

        self.U0_t2 = numpy.zeros((dim, self.t2s.length), dtype=REAL)
        self.U0_t1 = numpy.zeros((dim, self.t1s.length), dtype=REAL)
        self.U0_t3 = numpy.zeros((dim, self.t3s.length), dtype=REAL)
        self.U0_population = numpy.zeros((dim, population_time_axis.length), dtype=REAL)

        it1 = _axis_indices(self.t1s, population_time_axis)
        it2 = _axis_indices(self.t2s, population_time_axis)
        it3 = _axis_indices(self.t3s, population_time_axis)

        if dim != self.sys.Ntot:
            if dim == self.sys.Nb[1]:
                if self.sys.mult == 2:
                    self.U0fe_t3 = numpy.zeros(
                        (self.sys.Nb[2], dim, self.t3s.length), dtype=REAL
                    )
            elif self.sys.mult == 2:
                self.U0fe_t3 = numpy.zeros(
                    (self.sys.Nb[2], dim, self.t3s.length), dtype=REAL
                )
            else:
                raise Exception("Relaxation matrix has a wrong size: " + str(dim))

        # time independent rate matrix
        if len(KK.shape) == 2:
            #
            # Relaxation caused dephasing for single excitons
            #
            for aa in range(KK.shape[0]):
                if KK[aa, aa] <= 0.0:
                    self.U0_population[aa, :] = numpy.exp(
                        KK[aa, aa] * population_time_axis.data
                    )
                    coherence_amplitude = numpy.sqrt(self.U0_population[aa, :])
                    self.U0_t1[aa, :] = coherence_amplitude[it1]
                    self.U0_t2[aa, :] = coherence_amplitude[it2]
                    self.U0_t3[aa, :] = coherence_amplitude[it3]
                else:
                    raise QuantarheiError("Depopulation rate must be negative.")

            #
            # Relaxation caused dephasing for double-excitons
            #

            #
            # Finding population evolution matrix
            #

            # FIXME: Make sure it works with all t2s
            prop = PopulationPropagator(population_time_axis, self.KK)

            self.Uee = prop.get_PropagationMatrix(self.t2s)
            jumps = prop.get_JumpExpansion(
                self.t2s, max_order=self._get_jump_expansion_order()
            )

            self.U1_t2 = jumps[0]
            self.Usingle_t2 = jumps[1] if len(jumps) > 1 else numpy.zeros_like(self.Uee)
            self.Ujump_t2 = jumps
            self.Uremainder_t2 = self._get_remainder_propagator(jumps)
            self.Utransfer_t2 = self.Uremainder_t2

        if len(KK.shape) == 3:
            # time dependent rate matrix
            for aa in range(dim):
                if numpy.all(KK[:, aa, aa] <= 0.0):
                    cumulative: numpy.ndarray = numpy.zeros(
                        population_time_axis.length, dtype=REAL
                    )
                    for ii in range(population_time_axis.length - 1):
                        dt = (
                            population_time_axis.data[ii + 1]
                            - population_time_axis.data[ii]
                        )
                        cumulative[ii + 1] = (
                            cumulative[ii]
                            + 0.5 * (KK[ii, aa, aa] + KK[ii + 1, aa, aa]) * dt
                        )
                    self.U0_population[aa, :] = numpy.exp(cumulative)
                    coherence_amplitude = numpy.sqrt(self.U0_population[aa, :])
                    self.U0_t1[aa, :] = coherence_amplitude[it1]
                    self.U0_t2[aa, :] = coherence_amplitude[it2]
                    self.U0_t3[aa, :] = coherence_amplitude[it3]
                else:
                    raise Exception("Depopulation rate must be negative.")

            prop = PopulationPropagator(population_time_axis, self.KK)
            self.Uee = prop.get_PropagationMatrix(self.t2s)
            jumps = prop.get_JumpExpansion(
                self.t2s, max_order=self._get_jump_expansion_order()
            )

            self.U1_t2 = jumps[0]
            self.Usingle_t2 = jumps[1] if len(jumps) > 1 else numpy.zeros_like(self.Uee)
            self.Ujump_t2 = jumps
            self.Uremainder_t2 = self._get_remainder_propagator(jumps)
            self.Utransfer_t2 = self.Uremainder_t2


###############################################################################
#
#   DEPRECATED CODE BELOW
#
###############################################################################


class LiouvillePathway:
    """A single Liouville pathway (deprecated).

    Parameters
    ----------
    ptype : str
        Type of the Liouville pathway (e.g. ``"R1g"``, ``"R2g"``).
    """

    def __init__(self, ptype: str) -> None:

        self.ptype = ptype

        self.sign = 1.0
        #
        # Determine the response type
        # Types are:
        #   R ... rephasing
        #   NR .. non-rephasing
        #
        if ptype == "R2g":
            self.rtype = "R"
        elif ptype == "R1g":
            self.rtype = "NR"
        elif ptype == "R3g":
            self.rtype = "R"
        elif ptype == "R4g":
            self.rtype = "NR"
        elif ptype == "R1f":
            self.rtype = "R"
            self.sign = -1.0
        elif ptype == "R2f":
            self.rtype = "NR"
            self.sign = -1.0
        else:
            # unknown response type
            self.rtype = "U"

        self._frequencies_set = False
        self._rwa_set = False

        self.F4n = numpy.zeros(3)

    def set_dipoles(
        self,
        d1: numpy.ndarray,
        d2: numpy.ndarray | None = None,
        d3: numpy.ndarray | None = None,
        d4: numpy.ndarray | None = None,
    ) -> None:
        """Set the transition dipole moments of the response.

        Parameters
        ----------
        d1 : numpy.ndarray
            First transition dipole moment (3-vector).
        d2 : numpy.ndarray or None, optional
            Second transition dipole moment. If ``None``, defaults to ``d1``.
        d3 : numpy.ndarray or None, optional
            Third transition dipole moment. If ``None``, defaults to ``d1``.
        d4 : numpy.ndarray or None, optional
            Fourth transition dipole moment. If ``None``, defaults to ``d1``.
        """
        d = numpy.zeros((4, 3), dtype=float)

        d[0, :] = d1
        if d2 is None:
            d[1, :] = d1
        else:
            d[1, :] = d2
        if d3 is None:
            d[2, :] = d1
        else:
            d[2, :] = d3
        if d4 is None:
            d[3, :] = d1
        else:
            d[3, :] = d4

        self.F4n[0] = numpy.dot(d[3, :], d[2, :]) * numpy.dot(d[1, :], d[0, :])
        self.F4n[1] = numpy.dot(d[3, :], d[1, :]) * numpy.dot(d[2, :], d[0, :])
        self.F4n[2] = numpy.dot(d[3, :], d[0, :]) * numpy.dot(d[2, :], d[1, :])
        self.dd = d

    def set_frequencies(self, omega1: float, omega3: float) -> None:
        """Sets the frequencies of the response"""
        if self._frequencies_set:
            raise QuantarheiError("Frequencies of are already set.")

        self._omega1 = Manager().convert_energy_2_internal_u(omega1)
        self._omega3 = Manager().convert_energy_2_internal_u(omega3)

        self._frequencies_set = True

    def get_frequencies(self) -> tuple[float | numpy.ndarray, float | numpy.ndarray]:
        """Returns the two main frequencies of the response"""
        if self._frequencies_set:
            fr = (
                Manager().convert_energy_2_current_u(self._omega1),
                Manager().convert_energy_2_current_u(self._omega3),
            )

            return fr
        raise QuantarheiError("Frequencies not set.")

    def set_rwa(self, rwa: float) -> None:
        """Sets the RWA frequency"""
        if not self._frequencies_set:
            raise QuantarheiError("Frequencies must be set before setting RWA.")

        if not self._rwa_set:
            self._rwa = Manager().convert_energy_2_internal_u(rwa)
            self._omega1 = self._omega1 - self._rwa
            self._omega3 = self._omega3 - self._rwa

        else:
            raise QuantarheiError("RWA cannot be set twice. Reset first.")

        self._rwa_set = True

    def reset_rwa(self) -> None:
        """Resets the RWA setting"""
        if self._frequencies_set and self._rwa_set:
            self._omega1 = self._omega1 + self._rwa
            self._omega3 = self._omega3 + self._rwa

            self._rwa = 0.0
            self._rwa_set = False

    def get_rwa(self) -> float | numpy.ndarray:
        """Returns the RWA frequency"""
        return Manager().convert_energy_2_internal_u(self._rwa)


class ResponseFunction(LiouvillePathway):
    """Non-linear response function (deprecated).

    Parameters
    ----------
    ptype : str
        Type of the Liouville pathway (e.g. ``"R1g"``, ``"R2g"``).
    """

    def calculate_matrix(
        self, lab: Any, sys: Any, t2: float, t1s: Any, t3s: Any, rwa: float
    ) -> Any:
        """Calculates the matrix of response values in t1 and t3 times


        Parameters
        ----------
        lab :

        sys : ... None
            System for which the response is calculated. If None, we assume
            that this response is defined as a stand alone. Its dipole factors
            and frequencies are defined by its constructor.

        t2 : real
            The t2 time

        t1s :

        t2s :

        rwa: real
            Rotating wave frequency

        rmin :



        """
        om3 = self._omega3
        om1 = self._omega1

        # create the dipole factor
        dip = self.sign * numpy.dot(lab.F4eM4, self.F4n)

        # construct the phase factors in t1 and t3
        if self.rtype == "R":
            et13 = numpy.outer(numpy.exp(-1j * om3 * t3s), numpy.exp(1j * om1 * t1s))
        elif self.rtype == "NR":
            et13 = numpy.outer(numpy.exp(-1j * om3 * t3s), numpy.exp(-1j * om1 * t1s))

        val = self.func(t2, t3s, t1s, *self.args)

        return dip * val * et13

    def set_evaluation_function(self, func: Any) -> None:
        """Sets the function to be evaluated to get the response matrix"""
        self.func = func

    def set_auxliary_arguments(self, args: tuple[Any, ...]) -> None:
        """Sets additional arguments and their values for evaluation calls


        Parameters:
        -----------
        args : tuple
            A tuple of arguments that will be passed to the function after
            the expected standard arguments

        """
        self.args = args
