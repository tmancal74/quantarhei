"""Timing routines



"""
from __future__ import annotations

import datetime
import time

from ..core.managers import Manager

# Quantarhei imports
from .logging import printlog


def timeit(msg: str | None = None, show_stamp: bool = False, loglevel: int = 5, verbose: bool = True) -> None:
    """Start timing at this point and save the time

    """
    lconf = Manager().log_conf
    if msg is not None:
        printlog(msg, loglevel=loglevel, verbose=verbose)
    if show_stamp:
        printlog(f"Time stamp: {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
    lconf.time_stamp.put(time.time())


def untimeit(show_stamp: bool = False) -> float:
    """Stop timing and return current time

    """
    tm2 = time.time()
    if show_stamp:
        printlog(f"Time stamp: {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
    lconf = Manager().log_conf
    return tm2 - lconf.time_stamp.get()


def finished_in(show_stamp: bool = False, loglevel: int = 5, verbose: bool = True) -> None:
    """Print message with time past from the last timing statement

    """
    tm = untimeit()
    if show_stamp:
         printlog(f"... finished at {datetime.datetime.now():%Y-%m-%d %H:%M:%S} in",tm,"sec", loglevel=loglevel,
                 verbose=verbose)
    else:
        printlog("... finished in",tm,"sec", loglevel=loglevel,
                 verbose=verbose)


def done_in(show_stamp: bool = False, loglevel: int = 5, verbose: bool = True) -> None:
    """Print message with time past from the last timing statement

    """
    tm = untimeit()
    if show_stamp:
        printlog(f"... done at {datetime.datetime.now():%Y-%m-%d %H:%M:%S} in",tm,"sec", loglevel=loglevel,
                 verbose=verbose)
    else:
        printlog("... done in",tm,"sec", loglevel=loglevel, verbose=verbose)
