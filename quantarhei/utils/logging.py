# -*- coding: utf-8 -*-
"""Logging module of the Quantarhei package



"""
import quantarhei as qr


def printlog(*args, verbose=True, loglevel=0, **kwargs):
    """Prints logging information


    Parameters
    ----------

    *args : arguments
        arguments like in a print function

    verbose: bool
        If True, information will be logged, otherwise
        the function leaves immediately

    loglevel : int
        Information will be logged, if loglevel
        between 0 and 10 is smaller than the value
        of verbosity set in the Manager class


    Details
    -------

    The function prints loging information on the screen or into a file
    depending on global setting the Manager class. The behaviour
    is controlled by the following Manager attributes.

    verbosity : int
        The level of verbosity required to print; 0 means nothing
        will be printed, 10 means everything is logged

    log_on_screen : bool
        If True, logging information is printed on the screen

    log_to_file : bool
        If True, logging information is saved to a file

    log_file_opened : bool
        Set to True if logging file is opened
    log_file : file
        Handle of the logging file

    log_file_name : str
        Name of the logging file



    """

    if not verbose:
        return

    if (loglevel > 10) or (loglevel < 0):
        raise Exception("Loglevel must be between 0 and 10")


    manager = qr.Manager()
    if loglevel < manager.verbosity:

        if manager.log_on_screen:
            print(*args, **kwargs)

        if manager.log_to_file:
            if not manager.log_file_opened:
                manager.log_file = open(manager.log_file_name, "w")
                manager.log_file_opened = True
            manager.log_file.write(...)


def loglevels2bool(loglevs, verbose=False):
    """Converts a list of loglevels to a list of bools

    This function converts loglevels to books (True or False values).
    Its primary usage is in logging semi-critical regions  of vode by
    printlog() function efficiently. printlog() returns fast if verbose
    argument is set to False. Otherwise the function gets Manager instance,
    asks for current verbosity, and compares it with loglevel - a lot of work!
    If you want to avoid it say inside a nested loops, you better evaluate
    the loglevels in advance.


    Parameters
    ----------

    loglevels : list of int
        List of logleves to be converted

    verbose : bool (default False)
        If False, the function is left immediatel with all values converted
        to False

    Returns
    -------

    verb : list of bool
        List of bools, one fo each loglevel

    Examples
    --------

    Standard usage:

    >>> m = qr.Manager()
    >>> m.verbosity = 5
    >>> bools = loglevels2bool([0, 2, 5, 8, 10], verbose=True)
    >>> print(bools)
    [True, True, False, False, False]

    Do not forget to set verbose to True. It is False by default
    to avoid lengthy evaluation when not required

    >>> m = qr.Manager()
    >>> m.verbosity = 5
    >>> bools = loglevels2bool([0, 2, 5, 8, 10])
    >>> print(bools)
    [False, False, False, False, False]

    """

    verb = [False]*len(loglevs)

    if verbose:
        m = qr.Manager()
        k_v = 0
        for lev in loglevs:
            if m.verbosity > lev:
                verb[k_v] = True
            k_v += 1

    return verb
