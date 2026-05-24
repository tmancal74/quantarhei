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


def _jump_integration_steps(evol: Any) -> int:
    """Returns requested number of explicit s2 integration points."""
    if isinstance(evol, tuple) and len(evol) > 4 and isinstance(evol[4], dict):
        return int(evol[4].get("jump_integration_steps", 1))
    return 1


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
) -> numpy.ndarray:
    """Evaluates a placeholder one-jump response with a fixed jump time."""
    gg = system.get_lineshape_functions()
    old_defaults = getattr(gg, "_reset_defaults", None)
    reset_defaults = dict(old_defaults) if isinstance(old_defaults, dict) else {}
    reset_defaults["s2"] = s2
    gg._reset_defaults = reset_defaults

    try:
        return response_func(t2, t1, t3, lab, system, evol, KK)
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

    steps = _jump_integration_steps(evol)
    if steps <= 1 or t2 == 0.0:
        return _evaluate_single_jump_response(
            response_func, t2, 0.5 * t2, t1, t3, lab, system, evol, KK
        )

    s2s = numpy.linspace(0.0, t2, steps)
    weights = numpy.ones(steps, dtype=numpy.float64)
    weights[0] = 0.5
    weights[-1] = 0.5

    ret = None
    for weight, s2 in zip(weights, s2s):
        value = _evaluate_single_jump_response(
            response_func, t2, s2, t1, t3, lab, system, evol, KK
        )
        if ret is None:
            ret = weight * value
        else:
            ret += weight * value

    assert ret is not None
    return ret / (steps - 1)


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
                    * Ut2[a]
                    * Ut2[b]
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
                    * Ut2[a]
                    * Ut2[b]
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
                        * Ut2[a]
                        * Ut2[b]
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
                        * Ut2[a]
                        * Ut2[b]
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
                        * Ut2[a]
                        * Ut2[b]
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
                        * Ut2[a]
                        * Ut2[b]
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
