from __future__ import annotations

from typing import Any

import numpy as np

try:
    from qutip import (
        Qobj,
        liouvillian,
    )
    from qutip.solver.heom import (
        DrudeLorentzPadeBath,
        HEOMSolver,
    )

    _QUTIP_AVAILABLE = True
except ImportError:
    _QUTIP_AVAILABLE = False

import quantarhei as qr

_K_TO_INVCM = 0.6949
_INVFS_TO_INVCM = 1e15 / (2 * np.pi * 3e10)
_PS_TO_INVCM_TIME = 1e-12 * 2 * np.pi * 3e10

_DEFAULT_SOLVER_OPTIONS: dict[str, Any] = {
    "nsteps": 5000,
    "store_states": True,
    "rtol": 1e-12,
    "atol": 1e-12,
    "min_step": 1e-18,
    "method": "vern9",
    "progress_bar": "enhanced",
}


def prepare_simulation(Ham: Any, sbi: Any, depth: int) -> dict[str, Any]:
    if not _QUTIP_AVAILABLE:
        raise ImportError(
            "QuTiP is required for the HEOM bridge but is not installed. "
            "Install it with: pip install quantarhei[qutip]"
        )

    Nex = Ham.dim - 1

    with qr.energy_units("1/cm"):
        Hsys = Qobj(Ham.data)

    T = sbi.get_temperature() * _K_TO_INVCM
    Ltot = liouvillian(Hsys)
    baths = []

    for m in range(Nex):
        Q = Qobj(sbi.KK[m, :, :])
        cfce = sbi.get_correlation_function((m, m))
        with qr.energy_units("1/cm"):
            lam = cfce.get_reorganization_energy()
        gamma = (1.0 / cfce.get_correlation_time()) * _INVFS_TO_INVCM

        bath = DrudeLorentzPadeBath(Q, lam=lam, gamma=gamma, T=T, Nk=0, tag=str(m))
        baths.append(bath)
        _, terminator = bath.terminator()
        Ltot += terminator

    return {"depth": depth, "Ltot": Ltot, "baths": baths}


def run_simulation(
    kthprop: Any, rhoi: Any, options: dict[str, Any] | None = None
) -> Any:
    timea = kthprop.hy.sbi.get_time_axis()
    loc_options = options if options is not None else _DEFAULT_SOLVER_OPTIONS

    NC = kthprop.qutip_data["depth"]
    Ltot = kthprop.qutip_data["Ltot"]
    baths = kthprop.qutip_data["baths"]

    rho0 = Qobj(rhoi.data)
    Nt = timea.length
    end_time_ps = timea.data[Nt - 1] / 1000
    tlist = np.linspace(0, end_time_ps * _PS_TO_INVCM_TIME, Nt)

    solver = HEOMSolver(Ltot, baths, NC, options=loc_options)
    output = solver.run(rho0, tlist)

    rho_t = qr.DensityMatrixEvolution(timeaxis=timea, rhoi=rhoi)
    for t in range(Nt):
        rho_t.data[t, :, :] = output.states[t][:, :]

    return rho_t
