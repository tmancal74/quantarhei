# -*- coding: utf-8 -*-
"""Timing routines



"""
import time
import datetime

# Quantarhei imports 
from .logging import printlog
from ..core.managers import Manager


def timeit(msg=None,show_stamp=False, loglevel=5, verbose=True):
    """Start timing at this point and save the time
    
    """
    lconf = Manager().log_conf
    if msg is not None:
        printlog(msg, loglevel=loglevel, verbose=verbose)
    if show_stamp:
        printlog("Time stamp: {0:%Y-%m-%d %H:%M:%S}".format(
                 datetime.datetime.now()))
    lconf.time_stamp.put(time.time())


def untimeit(show_stamp=False):
    """Stop timing and return current time
    
    """
    tm2 = time.time()
    if show_stamp:
        printlog("Time stamp: {0:%Y-%m-%d %H:%M:%S}".format(
                 datetime.datetime.now()))
    lconf = Manager().log_conf
    return tm2 - lconf.time_stamp.get()


def finished_in(show_stamp=False,loglevel=5, verbose=True):
    """Print message with time past from the last timing statement
    
    """
    tm = untimeit()
    if show_stamp:
         printlog("... finished at {0:%Y-%m-%d %H:%M:%S} in".format(
                 datetime.datetime.now()),tm,"sec", loglevel=loglevel,
                 verbose=verbose)       
    else:
        printlog("... finished in",tm,"sec", loglevel=loglevel, 
                 verbose=verbose)


def done_in(show_stamp=False,loglevel=5, verbose=True):
    """Print message with time past from the last timing statement
    
    """
    tm = untimeit()
    if show_stamp:
        printlog("... done at {0:%Y-%m-%d %H:%M:%S} in".format(
                 datetime.datetime.now()),tm,"sec", loglevel=loglevel,
                 verbose=verbose)
    else:
        printlog("... done in",tm,"sec", loglevel=loglevel, verbose=verbose)
