"""Timing routines"""

from __future__ import annotations

import datetime
import time

from ..core.managers import Manager

# Quantarhei imports
from .logging import printlog


def timeit(
    msg: str | None = None,
    show_stamp: bool = False,
    loglevel: int = 5,
    verbose: bool = True,
) -> None:
    """Start a timing measurement and save the current time.

    Parameters
    ----------
    msg : str or None, optional
        Message to print via :func:`printlog` before recording the time.
        If ``None``, no message is printed.
    show_stamp : bool, optional
        If ``True``, also print a human-readable timestamp. Default is
        ``False``.
    loglevel : int, optional
        Log level passed to :func:`printlog`. Default is ``5``.
    verbose : bool, optional
        If ``False``, the message is suppressed. Default is ``True``.
    """
    lconf = Manager().log_conf
    if msg is not None:
        printlog(msg, loglevel=loglevel, verbose=verbose)
    if show_stamp:
        printlog(f"Time stamp: {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
    lconf.time_stamp.put(time.time())


def untimeit(show_stamp: bool = False) -> float:
    """Stop timing and return elapsed seconds since the last :func:`timeit` call.

    Parameters
    ----------
    show_stamp : bool, optional
        If ``True``, print a human-readable timestamp. Default is ``False``.

    Returns
    -------
    float
        Elapsed time in seconds.
    """
    tm2 = time.time()
    if show_stamp:
        printlog(f"Time stamp: {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
    lconf = Manager().log_conf
    return tm2 - lconf.time_stamp.get()


def finished_in(
    show_stamp: bool = False, loglevel: int = 5, verbose: bool = True
) -> None:
    """Print the elapsed time since the last :func:`timeit` call.

    Uses the phrasing ``"... finished in X sec"``.

    Parameters
    ----------
    show_stamp : bool, optional
        If ``True``, include a human-readable timestamp in the message.
        Default is ``False``.
    loglevel : int, optional
        Log level passed to :func:`printlog`. Default is ``5``.
    verbose : bool, optional
        If ``False``, the message is suppressed. Default is ``True``.
    """
    tm = untimeit()
    if show_stamp:
        printlog(
            f"... finished at {datetime.datetime.now():%Y-%m-%d %H:%M:%S} in",
            tm,
            "sec",
            loglevel=loglevel,
            verbose=verbose,
        )
    else:
        printlog("... finished in", tm, "sec", loglevel=loglevel, verbose=verbose)


def done_in(show_stamp: bool = False, loglevel: int = 5, verbose: bool = True) -> None:
    """Print the elapsed time since the last :func:`timeit` call.

    Uses the phrasing ``"... done in X sec"``.

    Parameters
    ----------
    show_stamp : bool, optional
        If ``True``, include a human-readable timestamp in the message.
        Default is ``False``.
    loglevel : int, optional
        Log level passed to :func:`printlog`. Default is ``5``.
    verbose : bool, optional
        If ``False``, the message is suppressed. Default is ``True``.
    """
    tm = untimeit()
    if show_stamp:
        printlog(
            f"... done at {datetime.datetime.now():%Y-%m-%d %H:%M:%S} in",
            tm,
            "sec",
            loglevel=loglevel,
            verbose=verbose,
        )
    else:
        printlog("... done in", tm, "sec", loglevel=loglevel, verbose=verbose)
