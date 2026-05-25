from __future__ import annotations

from typing import Any

import numpy

from .. import COMPLEX

ONE_JUMP_LINE_SHAPE_ARGUMENTS = (
    "t1",
    "t2",
    "s2",
    "t2-s2",
    "t3",
    "t1+t2",
    "t1+s2",
    "t1+t2-s2",
    "t1+t3",
    "t2+t3",
    "s2+t3",
    "t2-s2+t3",
    "t1+t2+t3",
    "t1+s2+t3",
    "t1+t2-s2+t3",
)


def _require_line_shape_arguments(gg: Any, response_name: str) -> None:
    """Checks that line-shape storage has all one-jump response arguments."""
    if not hasattr(gg, "time_mapping"):
        raise Exception(
            response_name
            + " requires FunctionStorage with one-jump line-shape arguments."
        )

    missing = [
        label for label in ONE_JUMP_LINE_SHAPE_ARGUMENTS if label not in gg.time_mapping
    ]
    if missing:
        raise Exception(
            response_name
            + " requires FunctionStorage(config=1); missing labels: "
            + ", ".join(missing)
        )


def _jump_time_metadata(evol: Any) -> dict[str, Any]:
    """Returns one-jump integration metadata from response evolution data."""
    if isinstance(evol, tuple) and len(evol) > 4 and isinstance(evol[4], dict):
        return evol[4]
    return {}


def _endpoint_t2_weight(evol: Any, left: int, right: int, default: Any) -> Any:
    """Returns an externally supplied t2 endpoint weight if available."""
    metadata = _jump_time_metadata(evol)
    endpoint = metadata.get("density_matrix_endpoint_t2", None)
    if endpoint is None:
        return default
    return endpoint[left, right]


def _single_jump_times(t2: Any, evol: Any) -> numpy.ndarray:
    """Returns jump times from the population time axis and graining."""
    t2_value = float(t2)
    if t2_value == 0.0:
        return numpy.array([0.0], dtype=numpy.float64)

    metadata = _jump_time_metadata(evol)
    graining = int(metadata.get("jump_time_graining", 1))
    axis = metadata.get("jump_time_axis", None)

    if axis is None:
        return numpy.array([0.5 * t2_value], dtype=numpy.float64)

    try:
        zero_index = axis.locate(0.0)[0]
    except Exception:
        zero_index = 0
    t2_index = int(metadata.get("jump_time_t2_index", axis.locate(t2_value)[0]))

    indices = list(range(zero_index, t2_index + 1, graining))
    if not indices or indices[-1] != t2_index:
        indices.append(t2_index)

    return numpy.asarray(
        axis.data[indices] - axis.data[zero_index], dtype=numpy.float64
    )


def _single_jump_transfer_matrix(t2: Any, s2: Any, evol: Any, KK: Any) -> Any:
    """Returns the one-jump transfer matrix density for a fixed jump time."""
    metadata = _jump_time_metadata(evol)
    axis = metadata.get("jump_time_axis", None)
    U0 = metadata.get("jump_zero_propagator", None)
    if axis is None or U0 is None or KK is None:
        return evol[3]

    zero_index = int(metadata.get("jump_time_zero_index", 0))
    t2_index = int(metadata.get("jump_time_t2_index", axis.locate(float(t2))[0]))
    s2_index = axis.locate(axis.data[zero_index] + float(s2))[0]

    rates = numpy.asarray(KK)
    if len(rates.shape) == 2:
        K_at_s2 = rates
    elif len(rates.shape) == 3:
        K_at_s2 = rates[min(s2_index, rates.shape[0] - 1)]
    else:
        return evol[3]

    transfer = numpy.zeros_like(K_at_s2)
    for final in range(K_at_s2.shape[0]):
        if U0[final, s2_index] == 0.0:
            final_survival = 0.0
        else:
            final_survival = U0[final, t2_index] / U0[final, s2_index]
        for initial in range(K_at_s2.shape[1]):
            if final != initial:
                transfer[final, initial] = (
                    K_at_s2[final, initial] * final_survival * U0[initial, s2_index]
                )

    return transfer


def _truncate_single_jump_times(
    t2: Any, s2s: numpy.ndarray, transfers: list[Any], evol: Any
) -> tuple[numpy.ndarray, list[Any]]:
    """Drop s2 points whose one-jump kernel is negligible."""
    metadata = _jump_time_metadata(evol)
    cutoff = float(metadata.get("jump_kernel_cutoff", 0.0))
    if cutoff <= 0.0 or len(s2s) <= 2:
        return s2s, transfers

    norms = numpy.asarray(
        [numpy.max(numpy.abs(transfer)) for transfer in transfers],
        dtype=numpy.float64,
    )
    max_norm = float(numpy.max(norms))
    if max_norm == 0.0:
        return numpy.array([0.0, float(t2)], dtype=numpy.float64), [
            transfers[0],
            transfers[-1],
        ]

    active = numpy.where(norms >= cutoff * max_norm)[0]
    if len(active) == 0:
        return numpy.array([0.0, float(t2)], dtype=numpy.float64), [
            transfers[0],
            transfers[-1],
        ]

    first = max(0, int(active[0]) - 1)
    last = min(len(s2s) - 1, int(active[-1]) + 1)
    return s2s[first : last + 1], transfers[first : last + 1]


def _single_jump_kernel_diagnostics(
    s2s: numpy.ndarray, transfers: list[Any], used_points: int, skipped: bool
) -> dict[str, Any]:
    """Return diagnostic information about the one-jump kernel."""
    if len(transfers) == 0:
        max_norm = 0.0
        integral_norm = 0.0
    else:
        norms = numpy.asarray(
            [numpy.max(numpy.abs(transfer)) for transfer in transfers],
            dtype=numpy.float64,
        )
        max_norm = float(numpy.max(norms))
        if len(norms) > 1:
            integral_norm = float(numpy.trapz(norms, s2s))
        else:
            integral_norm = 0.0

    return {
        "candidate_points": len(s2s),
        "used_points": int(used_points),
        "max_kernel_norm": max_norm,
        "integrated_kernel_norm": integral_norm,
        "skipped": skipped,
    }


def _set_single_jump_diagnostics(
    evol: Any, response_name: str, diagnostics: dict[str, Any]
) -> None:
    """Store one-jump diagnostics in response evolution metadata."""
    metadata = _jump_time_metadata(evol)
    all_diagnostics = metadata.get("jump_diagnostics", {})
    all_diagnostics[response_name] = diagnostics
    metadata["jump_diagnostics"] = all_diagnostics


def _replace_transfer_matrix(evol: Any, transfer_matrix: Any) -> Any:
    """Return response evolution data with a fixed transfer matrix."""
    if isinstance(evol, tuple) and len(evol) > 4:
        return evol[:3] + (transfer_matrix,) + evol[4:]
    if isinstance(evol, tuple) and len(evol) == 4:
        return evol[:3] + (transfer_matrix,)
    return evol


def _evaluate_single_jump_response(
    response_func: Any,
    t2: Any,
    s2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
    transfer_matrix: Any = None,
) -> numpy.ndarray:
    """Evaluates a placeholder one-jump response with a fixed jump time."""
    gg = system.get_lineshape_functions()
    if transfer_matrix is None:
        transfer_matrix = _single_jump_transfer_matrix(t2, s2, evol, KK)
    evol_at_s2 = _replace_transfer_matrix(evol, transfer_matrix)

    old_defaults = getattr(gg, "_reset_defaults", None)
    reset_defaults = dict(old_defaults) if isinstance(old_defaults, dict) else {}
    reset_defaults["s2"] = s2
    gg._reset_defaults = reset_defaults

    try:
        return response_func(t2, t1, t3, lab, system, evol_at_s2, KK)
    finally:
        if old_defaults is None:
            del gg._reset_defaults
        else:
            gg._reset_defaults = old_defaults


def _integrate_single_jump_response(
    response_func: Any,
    response_name: str,
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Integrates a placeholder one-jump response over the jump time ``s2``."""
    gg = system.get_lineshape_functions()
    _require_line_shape_arguments(gg, response_name)

    s2s = _single_jump_times(t2, evol)
    if t2 == 0.0:
        _set_single_jump_diagnostics(
            evol,
            response_name,
            _single_jump_kernel_diagnostics(s2s, [], 0, skipped=True),
        )
        return numpy.zeros((len(t3), len(t1)), dtype=COMPLEX)
    transfers = [_single_jump_transfer_matrix(t2, s2, evol, KK) for s2 in s2s]
    zero_cutoff = float(_jump_time_metadata(evol).get("jump_kernel_zero_cutoff", 0.0))
    original_s2s = s2s
    original_transfers = transfers
    if zero_cutoff > 0.0:
        max_norm = max(
            (float(numpy.max(numpy.abs(transfer))) for transfer in transfers),
            default=0.0,
        )
        if max_norm <= zero_cutoff:
            _set_single_jump_diagnostics(
                evol,
                response_name,
                _single_jump_kernel_diagnostics(
                    original_s2s, original_transfers, 0, skipped=True
                ),
            )
            return numpy.zeros((len(t3), len(t1)), dtype=COMPLEX)

    s2s, transfers = _truncate_single_jump_times(t2, s2s, transfers, evol)
    _set_single_jump_diagnostics(
        evol,
        response_name,
        _single_jump_kernel_diagnostics(
            original_s2s, original_transfers, len(s2s), skipped=False
        ),
    )
    if len(s2s) == 1:
        return _evaluate_single_jump_response(
            response_func, t2, s2s[0], t1, t3, lab, system, evol, KK, transfers[0]
        )

    previous_s2 = s2s[0]
    previous_value = _evaluate_single_jump_response(
        response_func, t2, previous_s2, t1, t3, lab, system, evol, KK, transfers[0]
    )
    ret = numpy.zeros_like(previous_value)
    for s2, transfer in zip(s2s[1:], transfers[1:]):
        value = _evaluate_single_jump_response(
            response_func, t2, s2, t1, t3, lab, system, evol, KK, transfer
        )
        ret += 0.5 * (previous_value + value) * (s2 - previous_s2)
        previous_s2 = s2
        previous_value = value

    return ret


def R1g(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Returns a matrix of the respose function values for given t1 and t3

    Parameters:
    -----------
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)

    t2 : float
        Value of the t2 (waiting) time of the response

    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)

    system : aggregate or molecule class
        An object storing all information about the system including
        the values of the line shape functions.


    """
    import numpy as np

    gg = system.get_lineshape_functions()
    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()
    # g = 0  # ground state index
    gg.create_data(reset={"t2": t2})

    band1 = system.get_band(1)
    band0 = system.get_band(0)

    Ut2 = evol[2]

    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: abba
    F4 = system.get_F4d("abba")
    dfac = np.einsum("i,abi->ab", lab.F4eM4, F4)

    ret: numpy.ndarray = np.zeros((len(t1), len(t3)), dtype=COMPLEX)
    for g in band0:
        for aa in band1:
            a = aa - 1
            for bb in band1:
                b = bb - 1

                # ret += Ut1[a,:][:,None]*Ut3[a,:][None,:]*\
                ret += (
                    dfac[b, a]
                    * Ut1[a, :][:, None]
                    * Ut3[a, :][None, :]
                    * _endpoint_t2_weight(evol, a, b, Ut2[a] * Ut2[b])
                    * np.exp(
                        -(np.einsum("i,ij", MM[b, a, :], gg[:, "t1"]))[:, None]
                        + (np.einsum("i,ij", MM[b, a, :], gg[:, "t1+t2"]))[:, None]
                        - (np.einsum("i,ijk", MM[a, a, :], gg[:, "t1+t2+t3"]))[:, :]
                        - np.conj(
                            (np.einsum("i,i", MM[b, b, :], gg[:, "t2"]))[None, None]
                        )
                        + np.conj(
                            (np.einsum("i,ij", MM[a, b, :], gg[:, "t2+t3"]))[None, :]
                        )
                        - np.conj(
                            (np.einsum("i,ij", MM[a, b, :], gg[:, "t3"]))[None, :]
                        )
                        - 1j * (En[aa] - En[g] - rwa) * t1[:, None]
                        - 1j * (En[aa] - En[bb]) * t2
                        - 1j * (En[aa] - En[g] - rwa) * t3[None, :]
                    )
                )

    return np.transpose(ret)


def R2g(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Returns a matrix of the respose function values for given t1 and t3

    Parameters:
    -----------
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)

    t2 : float
        Value of the t2 (waiting) time of the response

    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)

    system : aggregate or molecule class
        An object storing all information about the system including
        the values of the line shape functions.


    """
    import numpy as np

    gg = system.get_lineshape_functions()
    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()
    # g = 0  # ground state index
    gg.create_data(reset={"t2": t2})

    band1 = system.get_band(1)
    band0 = system.get_band(0)

    Ut2 = evol[2]

    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: baba
    F4 = system.get_F4d("baba")
    dfac = np.einsum("i,abi->ab", lab.F4eM4, F4)

    ret: numpy.ndarray = np.zeros((len(t1), len(t3)), dtype=COMPLEX)
    for g in band0:
        for aa in band1:
            a = aa - 1
            for bb in band1:
                b = bb - 1

                # ret += Ut1[a,:][:,None]*Ut3[b,:][None,:]*\
                ret += (
                    dfac[b, a]
                    * Ut1[a, :][:, None]
                    * Ut3[b, :][None, :]
                    * _endpoint_t2_weight(evol, a, b, Ut2[a] * Ut2[b])
                    * np.exp(
                        (np.einsum("i,i", MM[a, b, :], gg[:, "t2"]))[None, None]
                        - (np.einsum("i,ij", MM[b, b, :], gg[:, "t2+t3"]))[None, :]
                        - np.conj(
                            (np.einsum("i,ij", MM[a, a, :], gg[:, "t1+t2"]))[:, None]
                        )
                        - np.conj(
                            (np.einsum("i,ij", MM[a, b, :], gg[:, "t1"]))[:, None]
                        )
                        - np.conj(
                            (np.einsum("i,ij", MM[b, a, :], gg[:, "t3"]))[None, :]
                        )
                        + np.conj(
                            (np.einsum("i,ijk", MM[a, b, :], gg[:, "t1+t2+t3"]))[:, :]
                        )
                        - 1j * (En[g] - En[aa] + rwa) * t1[:, None]
                        - 1j * (En[bb] - En[aa]) * t2
                        - 1j * (En[bb] - En[g] - rwa) * t3[None, :]
                    )
                )

    return np.transpose(ret)


def R3g(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Returns a matrix of the respose function values for given t1 and t3

    Parameters:
    -----------
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)

    t2 : float
        Value of the t2 (waiting) time of the response

    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)

    system : aggregate or molecule class
        An object storing all information about the system including
        the values of the line shape functions.


    """
    import numpy as np

    gg = system.get_lineshape_functions()
    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()
    # g = 0  # ground state index
    gg.create_data(reset={"t2": t2})

    band1 = system.get_band(1)
    band0 = system.get_band(0)

    # Ut2 = evol[1]
    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: bbaa
    F4 = system.get_F4d("bbaa")
    dfac = np.einsum("i,abi->ab", lab.F4eM4, F4)

    ret: numpy.ndarray = np.zeros((len(t1), len(t3)), dtype=COMPLEX)
    for g in band0:
        for aa in band1:
            a = aa - 1
            for bb in band1:
                b = bb - 1

                # ret += Ut1[a,:][:,None]*Ut3[b,:][None,:]*\
                ret += (
                    dfac[b, a]
                    * Ut1[a, :][:, None]
                    * Ut3[b, :][None, :]
                    * np.exp(
                        -(np.einsum("i,ij", MM[b, b, :], gg[:, "t3"]))[None, :]
                        + np.conj(
                            (np.einsum("i,i", MM[b, a, :], gg[:, "t2"]))[None, None]
                        )
                        - np.conj(
                            (np.einsum("i,ij", MM[a, a, :], gg[:, "t1"]))[:, None]
                        )
                        - np.conj(
                            (np.einsum("i,ij", MM[a, b, :], gg[:, "t1+t2"]))[:, None]
                        )
                        - np.conj(
                            (np.einsum("i,ij", MM[b, a, :], gg[:, "t2+t3"]))[None, :]
                        )
                        + np.conj(
                            (np.einsum("i,ijk", MM[a, b, :], gg[:, "t1+t2+t3"]))[:, :]
                        )
                        - 1j * (En[g] - En[aa] + rwa) * t1[:, None]
                        - 1j * (En[g] - En[g]) * t2
                        - 1j * (En[bb] - En[g] - rwa) * t3[None, :]
                    )
                )

    return np.transpose(ret)


def R4g(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Returns a matrix of the respose function values for given t1 and t3

    Parameters:
    -----------
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)

    t2 : float
        Value of the t2 (waiting) time of the response

    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)

    system : aggregate or molecule class
        An object storing all information about the system including
        the values of the line shape functions.


    """
    import numpy as np

    gg = system.get_lineshape_functions()
    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()
    # g = 0  # ground state index
    gg.create_data(reset={"t2": t2})

    band1 = system.get_band(1)
    band0 = system.get_band(0)

    # Ut2 = evol[2]
    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: bbaa
    F4 = system.get_F4d("bbaa")
    dfac = np.einsum("i,abi->ab", lab.F4eM4, F4)

    ret: numpy.ndarray = np.zeros((len(t1), len(t3)), dtype=COMPLEX)
    for g in band0:
        for aa in band1:
            a = aa - 1
            for bb in band1:
                b = bb - 1

                # ret += Ut1[a,:][:,None]*Ut3[b,:][None,:]*\
                ret += (
                    dfac[b, a]
                    * Ut1[a, :][:, None]
                    * Ut3[b, :][None, :]
                    * np.exp(
                        -(np.einsum("i,i", MM[b, a, :], gg[:, "t2"]))[None, None]
                        - (np.einsum("i,ij", MM[a, a, :], gg[:, "t1"]))[:, None]
                        + (np.einsum("i,ij", MM[b, a, :], gg[:, "t1+t2"]))[:, None]
                        + (np.einsum("i,ij", MM[b, a, :], gg[:, "t2+t3"]))[None, :]
                        - (np.einsum("i,ij", MM[b, b, :], gg[:, "t3"]))[None, :]
                        - (np.einsum("i,ijk", MM[b, a, :], gg[:, "t1+t2+t3"]))[:, :]
                        - 1j * (En[aa] - En[g] - rwa) * t1[:, None]
                        - 1j * (En[g] - En[g]) * t2[None, None]
                        - 1j * (En[bb] - En[g] - rwa) * t3[None, :]
                    )
                )

    return np.transpose(ret)


def R1f(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Returns a matrix of the respose function values for given t1 and t3

    Parameters:
    -----------
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)

    t2 : float
        Value of the t2 (waiting) time of the response

    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)

    system : aggregate or molecule class
        An object storing all information about the system including
        the values of the line shape functions.


    """
    import numpy as np

    gg = system.get_lineshape_functions()
    gg.create_data(reset={"t2": t2})

    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()

    band0 = system.get_band(0)
    band1 = system.get_band(1)
    band2 = system.get_band(2)

    N0 = system.Nb[0]
    N1 = system.Nb[1]

    Ut2 = evol[2]

    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: fbfaba
    F4 = system.get_F4d("fbfaba")
    dfac = np.einsum("i,fabi->fab", lab.F4eM4, F4)

    ret: numpy.ndarray = np.zeros((len(t1), len(t3)), dtype=COMPLEX)
    for g in band0:
        for ff in band2:
            f = ff - N1 - N0
            p = ff - N0
            for aa in band1:
                a = aa - N0
                for bb in band1:
                    b = bb - N0

                    # ret += Ut1[a,:][:,None]*Ut3[b,:][None,:]*\
                    ret += (
                        -1.0
                        * dfac[f, b, a]
                        * Ut1[a, :][:, None]
                        * Ut3[b, :][None, :]
                        * _endpoint_t2_weight(evol, a, b, Ut2[a] * Ut2[b])
                        * np.exp(
                            -(np.einsum("i,ij", MM[a, a, :], gg[:, "t1+t2"]))[:, None]
                            - (np.einsum("i,ij", MM[b, a, :], gg[:, "t1"]))[:, None]
                            - (np.einsum("i,ij", MM[b, a, :], gg[:, "t3"]))[None, :]
                            + (np.einsum("i,ij", MM[b, p, :], gg[:, "t3"]))[None, :]
                            + (np.einsum("i,ij", MM[p, a, :], gg[:, "t1+t2"]))[:, None]
                            + (np.einsum("i,ij", MM[p, a, :], gg[:, "t3"]))[None, :]
                            - (np.einsum("i,ij", MM[p, p, :], gg[:, "t3"]))[None, :]
                            + (np.einsum("i,ijk", MM[b, a, :], gg[:, "t1+t2+t3"]))[:, :]
                            - (np.einsum("i,ijk", MM[p, a, :], gg[:, "t1+t2+t3"]))[:, :]
                            + np.conj(np.einsum("i,i", MM[a, b, :], gg[:, "t2"]))
                            - np.conj(np.einsum("i,i", MM[p, b, :], gg[:, "t2"]))
                            - np.conj(
                                (np.einsum("i,ij", MM[b, b, :], gg[:, "t2+t3"]))[
                                    None, :
                                ]
                            )
                            + np.conj(
                                (np.einsum("i,ij", MM[p, b, :], gg[:, "t2+t3"]))[
                                    None, :
                                ]
                            )
                            - 1j * (En[aa] - En[g] - rwa) * t1[:, None]
                            - 1j * (En[aa] - En[bb]) * t2[None, None]
                            - 1j * (En[ff] - En[bb] - rwa) * t3[None, :]
                        )
                    )

    return np.transpose(ret)


def R1f_wrong(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Returns a matrix of the respose function values for given t1 and t3

    Parameters:
    -----------
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)

    t2 : float
        Value of the t2 (waiting) time of the response

    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)

    system : aggregate or molecule class
        An object storing all information about the system including
        the values of the line shape functions.


    """
    import numpy as np

    gg = system.get_lineshape_functions()
    gg.create_data(reset={"t2": t2})

    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()

    band0 = system.get_band(0)
    band1 = system.get_band(1)
    band2 = system.get_band(2)

    N0 = system.Nb[0]
    N1 = system.Nb[1]

    Ut2 = evol[2]

    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: fbfaba
    F4 = system.get_F4d("fbfaba")
    dfac = np.einsum("i,fabi->fab", lab.F4eM4, F4)

    ret: numpy.ndarray = np.zeros((len(t1), len(t3)), dtype=COMPLEX)
    for g in band0:
        for ff in band2:
            f = ff - N1 - N0
            p = ff - N0
            for aa in band1:
                a = aa - N0
                for bb in band1:
                    b = bb - N0

                    # ret += Ut1[a,:][:,None]*Ut3[b,:][None,:]*\
                    ret += (
                        -1.0
                        * dfac[f, b, a]
                        * Ut1[a, :][:, None]
                        * Ut3[b, :][None, :]
                        * _endpoint_t2_weight(evol, a, b, Ut2[a] * Ut2[b])
                        * np.exp(
                            +(np.einsum("i,ij", MM[a, a, :], gg[:, "t1+t2"]))[:, None]
                            - (np.einsum("i,ij", MM[b, a, :], gg[:, "t1"]))[:, None]
                            - (np.einsum("i,ij", MM[b, a, :], gg[:, "t3"]))[None, :]
                            + (np.einsum("i,ij", MM[b, p, :], gg[:, "t3"]))[None, :]
                            - (np.einsum("i,ij", MM[p, a, :], gg[:, "t1+t2"]))[:, None]
                            + (np.einsum("i,ij", MM[p, a, :], gg[:, "t3"]))[None, :]
                            - (np.einsum("i,ij", MM[p, p, :], gg[:, "t3"]))[None, :]
                            - (np.einsum("i,ijk", MM[b, a, :], gg[:, "t1+t2+t3"]))[:, :]
                            - (np.einsum("i,ijk", MM[p, a, :], gg[:, "t1+t2+t3"]))[:, :]
                            - np.conj(np.einsum("i,i", MM[a, b, :], gg[:, "t2"]))
                            + np.conj(np.einsum("i,i", MM[p, b, :], gg[:, "t2"]))
                            - np.conj(
                                (np.einsum("i,ij", MM[b, b, :], gg[:, "t2+t3"]))[
                                    None, :
                                ]
                            )
                            - np.conj(
                                (np.einsum("i,ij", MM[p, b, :], gg[:, "t2+t3"]))[
                                    None, :
                                ]
                            )
                            - 1j * (En[aa] - En[g] - rwa) * t1[:, None]
                            - 1j * (En[aa] - En[bb]) * t2[None, None]
                            - 1j * (En[ff] - En[bb] - rwa) * t3[None, :]
                        )
                    )

    return np.transpose(ret)


def R2f(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Returns a matrix of the respose function values for given t1 and t3

    Parameters:
    -----------
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)

    t2 : float
        Value of the t2 (waiting) time of the response

    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)

    system : aggregate or molecule class
        An object storing all information about the system including
        the values of the line shape functions.


    """
    import numpy as np

    gg = system.get_lineshape_functions()
    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()

    gg.create_data(reset={"t2": t2})

    band0 = system.get_band(0)
    band1 = system.get_band(1)
    band2 = system.get_band(2)

    N0 = system.Nb[0]
    N1 = system.Nb[1]

    Ut2 = evol[2]

    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: fafbba
    F4 = system.get_F4d("fafbba")
    dfac = np.einsum("i,fabi->fab", lab.F4eM4, F4)

    ret: numpy.ndarray = np.zeros((len(t1), len(t3)), dtype=COMPLEX)

    for g in band0:
        for ff in band2:
            f = ff - N1 - N0
            p = ff - N0
            for aa in band1:
                a = aa - N0
                for bb in band1:
                    b = bb - N0

                    # ret += Ut1[a,:][:,None]*Ut3[a,:][None,:]*\
                    ret += (
                        -1.0
                        * dfac[f, b, a]
                        * Ut1[a, :][:, None]
                        * Ut3[a, :][None, :]
                        * _endpoint_t2_weight(evol, a, b, Ut2[a] * Ut2[b])
                        * np.exp(
                            -(np.einsum("i,i", MM[b, b, :], gg[:, "t2"]))[None, None]
                            + (np.einsum("i,i", MM[p, b, :], gg[:, "t2"]))[None, None]
                            + (np.einsum("i,ij", MM[a, b, :], gg[:, "t2+t3"]))[None, :]
                            - (np.einsum("i,ij", MM[a, b, :], gg[:, "t3"]))[None, :]
                            + (np.einsum("i,ij", MM[a, p, :], gg[:, "t3"]))[None, :]
                            - (np.einsum("i,ij", MM[p, b, :], gg[:, "t2+t3"]))[None, :]
                            + (np.einsum("i,ij", MM[p, b, :], gg[:, "t3"]))[None, :]
                            - (np.einsum("i,ij", MM[p, p, :], gg[:, "t3"]))[None, :]
                            - np.conj(
                                (np.einsum("i,ij", MM[a, b, :], gg[:, "t1"]))[:, None]
                            )
                            + np.conj(
                                (np.einsum("i,ij", MM[a, b, :], gg[:, "t1+t2"]))[
                                    :, None
                                ]
                            )
                            - np.conj(
                                (np.einsum("i,ij", MM[a, p, :], gg[:, "t1+t2"]))[
                                    :, None
                                ]
                            )
                            - np.conj(
                                (np.einsum("i,ijk", MM[a, a, :], gg[:, "t1+t2+t3"]))[
                                    :, :
                                ]
                            )
                            + np.conj(
                                (np.einsum("i,ijk", MM[a, p, :], gg[:, "t1+t2+t3"]))[
                                    :, :
                                ]
                            )
                            - 1j * (En[g] - En[aa] + rwa) * t1[:, None]
                            - 1j * (En[bb] - En[aa]) * t2[None, None]
                            - 1j * (En[ff] - En[aa] - rwa) * t3[None, :]
                        )
                    )

    return np.transpose(ret)


def R2f_wrong(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Returns a matrix of the respose function values for given t1 and t3

    Parameters:
    -----------
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)

    t2 : float
        Value of the t2 (waiting) time of the response

    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)

    system : aggregate or molecule class
        An object storing all information about the system including
        the values of the line shape functions.


    """
    import numpy as np

    gg = system.get_lineshape_functions()
    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()

    gg.create_data(reset={"t2": t2})

    band0 = system.get_band(0)
    band1 = system.get_band(1)
    band2 = system.get_band(2)

    N0 = system.Nb[0]
    N1 = system.Nb[1]

    Ut2 = evol[2]

    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: fafbba
    F4 = system.get_F4d("fafbba")
    dfac = np.einsum("i,fabi->fab", lab.F4eM4, F4)

    ret: numpy.ndarray = np.zeros((len(t1), len(t3)), dtype=COMPLEX)

    for g in band0:
        for ff in band2:
            f = ff - N1 - N0
            p = ff - N0
            for aa in band1:
                a = aa - N0
                for bb in band1:
                    b = bb - N0

                    # ret += Ut1[a,:][:,None]*Ut3[a,:][None,:]*\
                    ret += (
                        -1.0
                        * dfac[f, b, a]
                        * Ut1[a, :][:, None]
                        * Ut3[a, :][None, :]
                        * _endpoint_t2_weight(evol, a, b, Ut2[a] * Ut2[b])
                        * np.exp(
                            +(np.einsum("i,i", MM[b, b, :], gg[:, "t2"]))[None, None]
                            - (np.einsum("i,i", MM[p, b, :], gg[:, "t2"]))[None, None]
                            - (np.einsum("i,ij", MM[a, b, :], gg[:, "t2+t3"]))[None, :]
                            + (np.einsum("i,ij", MM[a, b, :], gg[:, "t3"]))[None, :]
                            + (np.einsum("i,ij", MM[a, p, :], gg[:, "t3"]))[None, :]
                            - (np.einsum("i,ij", MM[p, b, :], gg[:, "t2+t3"]))[None, :]
                            + (np.einsum("i,ij", MM[p, b, :], gg[:, "t3"]))[None, :]
                            - (np.einsum("i,ij", MM[p, p, :], gg[:, "t3"]))[None, :]
                            - np.conj(
                                (np.einsum("i,ij", MM[a, b, :], gg[:, "t1"]))[:, None]
                            )
                            - np.conj(
                                (np.einsum("i,ij", MM[a, b, :], gg[:, "t1+t2"]))[
                                    :, None
                                ]
                            )
                            + np.conj(
                                (np.einsum("i,ij", MM[a, p, :], gg[:, "t1+t2"]))[
                                    :, None
                                ]
                            )
                            - np.conj(
                                (np.einsum("i,ijk", MM[a, a, :], gg[:, "t1+t2+t3"]))[
                                    :, :
                                ]
                            )
                            - np.conj(
                                (np.einsum("i,ijk", MM[a, p, :], gg[:, "t1+t2+t3"]))[
                                    :, :
                                ]
                            )
                            - 1j * (En[g] - En[aa] + rwa) * t1[:, None]
                            - 1j * (En[bb] - En[aa]) * t2[None, None]
                            - 1j * (En[ff] - En[aa] - rwa) * t3[None, :]
                        )
                    )

    return np.transpose(ret)


def R1g_scM0g(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Returns a matrix of the respose function values for given t1 and t3

    Parameters:
    -----------
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)

    t2 : float
        Value of the t2 (waiting) time of the response

    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)

    system : aggregate or molecule class
        An object storing all information about the system including
        the values of the line shape functions.


    """
    import numpy as np

    gg = system.get_lineshape_functions()
    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()
    # g = 0  # ground state index
    gg.create_data(reset={"t2": t2})

    band1 = system.get_band(1)
    band0 = system.get_band(0)

    Ut2 = evol[3]

    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: abba
    F4 = system.get_F4d("bbaa")
    dfac = np.einsum("i,bai->ba", lab.F4eM4, F4)

    ret: numpy.ndarray = np.zeros((len(t1), len(t3)), dtype=COMPLEX)
    for g in band0:
        for aa in band1:
            a = aa - 1
            for bb in band1:
                b = bb - 1

                if a != b:
                    # ret += Ut1[a,:][:,None]*Ut3[b,:][None,:]*\
                    ret += (
                        dfac[b, a]
                        * Ut1[a, :][:, None]
                        * Ut3[b, :][None, :]
                        * Ut2[b, a]
                        * np.exp(
                            -(np.einsum("i,ij", MM[a, a, :], gg[:, "t1"]))[:, None]
                            - (np.einsum("i,ij", MM[b, b, :], gg[:, "t3"]))[None, :]
                            - 1j * (En[aa] - En[g] - rwa) * t1[:, None]
                            - 1j * (En[bb] - En[g] - rwa) * t3[None, :]
                        )
                    )

    return np.transpose(ret)


def R2g_scM0g(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Returns a matrix of the respose function values for given t1 and t3

    Parameters:
    -----------
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)

    t2 : float
        Value of the t2 (waiting) time of the response

    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)

    system : aggregate or molecule class
        An object storing all information about the system including
        the values of the line shape functions.


    """
    import numpy as np

    gg = system.get_lineshape_functions()
    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()
    # g = 0  # ground state index
    gg.create_data(reset={"t2": t2})

    band1 = system.get_band(1)
    band0 = system.get_band(0)

    Ut2 = evol[3]

    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: baba
    F4 = system.get_F4d("bbaa")
    dfac = np.einsum("i,bai->ba", lab.F4eM4, F4)

    ret: numpy.ndarray = np.zeros((len(t1), len(t3)), dtype=COMPLEX)
    for g in band0:
        for aa in band1:
            a = aa - 1
            for bb in band1:
                b = bb - 1

                if a != b:
                    # ret += Ut1[a,:][:,None]*Ut3[b,:][None,:]*\
                    ret += (
                        dfac[b, a]
                        * Ut1[a, :][:, None]
                        * Ut3[b, :][None, :]
                        * Ut2[b, a]
                        * np.exp(
                            -np.conj(
                                (np.einsum("i,ij", MM[a, a, :], gg[:, "t1"]))[:, None]
                            )
                            - (np.einsum("i,ij", MM[b, b, :], gg[:, "t3"]))[None, :]
                            - 1j * (En[g] - En[aa] + rwa) * t1[:, None]
                            - 1j * (En[bb] - En[g] - rwa) * t3[None, :]
                        )
                    )

    return np.transpose(ret)


def R1f_scM0g(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Returns a matrix of the respose function values for given t1 and t3

    Parameters:
    -----------
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)

    t2 : float
        Value of the t2 (waiting) time of the response

    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)

    system : aggregate or molecule class
        An object storing all information about the system including
        the values of the line shape functions.


    """
    import numpy as np

    gg = system.get_lineshape_functions()
    gg.create_data(reset={"t2": t2})

    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()

    band0 = system.get_band(0)
    band1 = system.get_band(1)
    band2 = system.get_band(2)

    N0 = system.Nb[0]
    N1 = system.Nb[1]

    Ut2 = evol[3]

    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: fbfaba
    F4 = system.get_F4d("fbfbaa")
    dfac = np.einsum("i,fbai->fba", lab.F4eM4, F4)

    ret: numpy.ndarray = np.zeros((len(t1), len(t3)), dtype=COMPLEX)
    for g in band0:
        for ff in band2:
            f = ff - N1 - N0
            p = ff - N0
            for aa in band1:
                a = aa - N0
                for bb in band1:
                    b = bb - N0

                    if a != b:
                        # ret += Ut1[a,:][:,None]*Ut3[a,:][None,:]*\
                        ret += (
                            (-1.0)
                            * dfac[f, a, b]
                            * Ut1[a, :][:, None]
                            * Ut3[b, :][None, :]
                            * Ut2[b, a]
                            * np.exp(
                                -(np.einsum("i,ij", MM[a, a, :], gg[:, "t1"]))[:, None]
                                - (
                                    np.conj(np.einsum("i,ij", MM[b, b, :], gg[:, "t3"]))
                                )[None, :]
                                - (np.einsum("i,ij", MM[p, p, :], gg[:, "t3"]))[None, :]
                                + (
                                    2.0
                                    * np.real(
                                        np.einsum("i,ij", MM[p, b, :], gg[:, "t3"])
                                    )
                                )[None, :]
                                - 1j * (En[aa] - En[g] - rwa) * t1[:, None]
                                - 1j * (En[ff] - En[bb] - rwa) * t3[None, :]
                            )
                        )

    return np.transpose(ret)


def R2f_scM0g(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Returns a matrix of the respose function values for given t1 and t3

    Parameters:
    -----------
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)

    t2 : float
        Value of the t2 (waiting) time of the response

    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)

    system : aggregate or molecule class
        An object storing all information about the system including
        the values of the line shape functions.


    """
    import numpy as np

    gg = system.get_lineshape_functions()
    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()

    gg.create_data(reset={"t2": t2})

    band0 = system.get_band(0)
    band1 = system.get_band(1)
    band2 = system.get_band(2)

    N0 = system.Nb[0]
    N1 = system.Nb[1]

    Ut2 = evol[3]

    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: fafbba
    F4 = system.get_F4d("fbfbaa")
    dfac = np.einsum("i,fbai->fba", lab.F4eM4, F4)

    ret: numpy.ndarray = np.zeros((len(t1), len(t3)), dtype=COMPLEX)

    for g in band0:
        for ff in band2:
            f = ff - N1 - N0
            p = ff - N0
            for aa in band1:
                a = aa - N0
                for bb in band1:
                    b = bb - N0

                    if b != a:
                        # ret +=  Ut1[a,:][:,None]*Ut3[a,:][None,:]*\
                        ret += (
                            (-1.0)
                            * dfac[f, a, b]
                            * Ut1[a, :][:, None]
                            * Ut3[b, :][None, :]
                            * Ut2[b, a]
                            * np.exp(
                                -np.conj(
                                    (np.einsum("i,ij", MM[a, a, :], gg[:, "t1"]))[
                                        :, None
                                    ]
                                )
                                - (
                                    np.conj(np.einsum("i,ij", MM[b, b, :], gg[:, "t3"]))
                                )[None, :]
                                - (np.einsum("i,ij", MM[p, p, :], gg[:, "t3"]))[None, :]
                                + (
                                    2.0
                                    * np.real(
                                        np.einsum("i,ij", MM[p, b, :], gg[:, "t3"])
                                    )
                                )[None, :]
                                - 1j * (En[g] - En[aa] + rwa) * t1[:, None]
                                - 1j * (En[ff] - En[bb] - rwa) * t3[None, :]
                            )
                        )

                        # print("---")
                        # if np.abs(dfac[f,a,b]) > 1.0e-12:
                        #     print(En[aa] - En[g])
                        #     print(ff, "=", system.twoex_indx[ff,0], system.twoex_indx[ff,1])
                        #     print(ff, bb, En[ff] - En[bb])

    return np.transpose(ret)


def R1f_scM0e(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Returns a matrix of the respose function values for given t1 and t3

    Parameters:
    -----------
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)

    t2 : float
        Value of the t2 (waiting) time of the response

    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)

    system : aggregate or molecule class
        An object storing all information about the system including
        the values of the line shape functions.


    """
    import numpy as np

    gg = system.get_lineshape_functions()
    gg.create_data(reset={"t2": t2})

    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()

    band0 = system.get_band(0)
    band1 = system.get_band(1)
    band2 = system.get_band(2)

    N0 = system.Nb[0]
    N1 = system.Nb[1]

    Ut2 = evol[3]

    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: fbfaba
    F4 = system.get_F4d("fbfbaa")
    dfac = np.einsum("i,fbai->fba", lab.F4eM4, F4)

    ret: numpy.ndarray = np.zeros((len(t1), len(t3)), dtype=COMPLEX)
    for g in band0:
        for ff in band2:
            f = ff - N1 - N0
            p = ff - N0
            for aa in band1:
                a = aa - N0
                for bb in band1:
                    b = bb - N0

                    if a != b:
                        # ret += Ut1[a,:][:,None]*Ut3[a,:][None,:]*\
                        ret += (
                            (-1.0)
                            * dfac[f, a, b]
                            * Ut1[a, :][:, None]
                            * Ut3[b, :][None, :]
                            * Ut2[b, a]
                            * np.exp(
                                -(np.einsum("i,ij", MM[a, a, :], gg[:, "t1"]))[:, None]
                                - (
                                    np.conj(np.einsum("i,ij", MM[b, b, :], gg[:, "t3"]))
                                )[None, :]
                                - (np.einsum("i,ij", MM[p, p, :], gg[:, "t3"]))[None, :]
                                + (
                                    2.0
                                    * np.real(
                                        np.einsum("i,ij", MM[p, b, :], gg[:, "t3"])
                                    )
                                )[None, :]
                                - 1j * (En[aa] - En[g] - rwa) * t1[:, None]
                                - 1j * (En[ff] - En[bb] - rwa) * t3[None, :]
                            )
                        )

    return np.transpose(ret)


def R2f_scM0e(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Returns a matrix of the respose function values for given t1 and t3

    Parameters:
    -----------
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)

    t2 : float
        Value of the t2 (waiting) time of the response

    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)

    system : aggregate or molecule class
        An object storing all information about the system including
        the values of the line shape functions.


    """
    import numpy as np

    gg = system.get_lineshape_functions()
    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()

    gg.create_data(reset={"t2": t2})

    band0 = system.get_band(0)
    band1 = system.get_band(1)
    band2 = system.get_band(2)

    N0 = system.Nb[0]
    N1 = system.Nb[1]

    Ut2 = evol[3]

    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: fafbba
    F4 = system.get_F4d("fbfbaa")
    dfac = np.einsum("i,fbai->fba", lab.F4eM4, F4)

    ret: numpy.ndarray = np.zeros((len(t1), len(t3)), dtype=COMPLEX)

    lam = gg.get_reorganization_energies()
    # print("lam = ", convert(lam,"int","1/cm"))

    for g in band0:
        for ff in band2:
            f = ff - N1 - N0
            p = ff - N0
            for aa in band1:
                a = aa - N0
                for bb in band1:
                    b = bb - N0

                    if b != a:
                        # ret +=  Ut1[a,:][:,None]*Ut3[a,:][None,:]*\
                        ret += (
                            (-1.0)
                            * dfac[f, a, b]
                            * Ut1[a, :][:, None]
                            * Ut3[b, :][None, :]
                            * Ut2[b, a]
                            * np.exp(
                                -np.conj(
                                    (np.einsum("i,ij", MM[a, a, :], gg[:, "t1"]))[
                                        :, None
                                    ]
                                )
                                - (
                                    np.conj(np.einsum("i,ij", MM[b, b, :], gg[:, "t3"]))
                                )[None, :]
                                - (np.einsum("i,ij", MM[p, p, :], gg[:, "t3"]))[None, :]
                                + (
                                    2.0
                                    * np.real(
                                        np.einsum("i,ij", MM[p, b, :], gg[:, "t3"])
                                    )
                                )[None, :]
                                - 1j * (En[g] - En[aa] + rwa) * t1[:, None]
                                - 1j * (En[ff] - En[bb] - rwa) * t3[None, :]
                            )
                        )

                        # print("---")
                        # if np.abs(dfac[f,a,b]) > 1.0e-12:
                        #     print(En[aa] - En[g])
                        #     print(ff, "=", system.twoex_indx[ff,0], system.twoex_indx[ff,1])
                        #     print(ff, bb, En[ff] - En[bb])

    return np.transpose(ret)


def R1g_scM1g(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Single-jump transfer version of ``R1g_scM0g``.

    This currently uses the same formula as the remainder pathway.  The
    separate implementation exists so it can be replaced by the explicit
    one-jump response formula without changing calculator bookkeeping.
    """
    return _integrate_single_jump_response(
        R1g_scM0g, "R1g_scM1g", t2, t1, t3, lab, system, evol, KK
    )


def R2g_scM1g(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Single-jump transfer version of ``R2g_scM0g``."""
    return _integrate_single_jump_response(
        R2g_scM0g, "R2g_scM1g", t2, t1, t3, lab, system, evol, KK
    )


def R1f_scM1g(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Single-jump transfer version of ``R1f_scM0g``."""
    return _integrate_single_jump_response(
        R1f_scM0g, "R1f_scM1g", t2, t1, t3, lab, system, evol, KK
    )


def R2f_scM1g(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Single-jump transfer version of ``R2f_scM0g``."""
    return _integrate_single_jump_response(
        R2f_scM0g, "R2f_scM1g", t2, t1, t3, lab, system, evol, KK
    )


def R1f_scM1e(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Single-jump transfer version of ``R1f_scM0e``."""
    return _integrate_single_jump_response(
        R1f_scM0e, "R1f_scM1e", t2, t1, t3, lab, system, evol, KK
    )


def R2f_scM1e(
    t2: Any,
    t1: numpy.ndarray,
    t3: numpy.ndarray,
    lab: Any,
    system: Any,
    evol: Any,
    KK: Any,
) -> numpy.ndarray:
    """Single-jump transfer version of ``R2f_scM0e``."""
    return _integrate_single_jump_response(
        R2f_scM0e, "R2f_scM1e", t2, t1, t3, lab, system, evol, KK
    )


#
# REGISTRATION CODE
#

dc = dict()

dc["R1g"] = R1g
dc["R2g"] = R2g
dc["R3g"] = R3g
dc["R4g"] = R4g
dc["R1f"] = R1f
dc["R2f"] = R2f
dc["R1g_scM0g"] = R1g_scM0g
dc["R2g_scM0g"] = R2g_scM0g
dc["R1f_scM0g"] = R1f_scM0g
dc["R2f_scM0g"] = R2f_scM0g
dc["R1f_scM0e"] = R1f_scM0e
dc["R2f_scM0e"] = R2f_scM0e
dc["R1g_scM1g"] = R1g_scM1g
dc["R2g_scM1g"] = R2g_scM1g
dc["R1f_scM1g"] = R1f_scM1g
dc["R2f_scM1g"] = R2f_scM1g
dc["R1f_scM1e"] = R1f_scM1e
dc["R2f_scM1e"] = R2f_scM1e


def get_implementation(name: str) -> Any:
    """Returns a dictionary of functions"""
    return dc[name]
