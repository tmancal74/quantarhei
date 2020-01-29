# -*- coding: utf-8 -*-
"""Logging module of the Quantarhei package



"""
import traceback
import quantarhei as qr

def init_logging():
    """Initialization of logging
    
    We test if the logging is parallel or not
    
    """
    manager = qr.Manager().log_conf
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        if size == 1:
            manager.is_serial = True
        else:
            manager.is_serial = False
            manager.log_file_appendix = "."+str(rank)
    except:
        manager.is_serial = True
    manager.initialized = True
    
    
def log_urgent(*args, **kwargs):
    printlog(*args, loglevel=qr.LOG_URGENT, **kwargs)


def log_report(*args, **kwargs):
    printlog(*args, loglevel=qr.LOG_REPORT, **kwargs)


def log_info(*args, **kwargs):
    printlog(*args, loglevel=qr.LOG_INFO, **kwargs)

        
def log_detail(*args, **kwargs):
    printlog(*args, loglevel=qr.LOG_DETAIL, **kwargs)

    
def log_quick(*args, verbose=True, **kwargs):
    if not verbose:
        return
    printlog(*args, loglevel=qr.LOG_QUICK, **kwargs)  


def printlog(*args, verbose=True, loglevel=5, 
             incr_indent=0, use_indent=True, **kwargs):
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


    manager = qr.Manager().log_conf
    if not manager.initialized:
        init_logging()
    
    if not manager.verbose:
        return
    
    manager.log_indent += incr_indent
    
    if loglevel <= manager.verbosity:

        if manager.log_on_screen:
            if use_indent:
                indent = " "*manager.log_indent
            else:
                indent = ""
            print(indent, *args, **kwargs)

    if loglevel <= manager.fverbosity:
        
        if manager.log_to_file:
            if not manager.log_file_opened:
                manager.log_file = open(manager.log_file_name
                                        +manager.log_file_appendix, "w")
                manager.log_file_opened = True
            if use_indent:
                indent = " "*manager.log_indent
            else:
                indent = ""
            if "end" in kwargs.keys():
                kwargs["end"] = None
            print(indent, *args, **kwargs, file=manager.log_file)


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
        m = qr.Manager().log_conf
        k_v = 0
        for lev in loglevs:
            if m.verbosity > lev:
                verb[k_v] = True
            k_v += 1

    return verb


def tprint(var, messg=None, default=None):
    """Test the existence of the variable with the name `var`
    
    
    If the default is specified, the non-existence of the variable
    is not a problem, and the variable is set to the default
    
    """    
    try:
        if messg is not None:
            print("#", messg)
        val = eval(globals()[var])
        if isinstance(val, str):
            val = '"'+val+'"'
        print(var, "=", val)

    except:
        if default is None:
            traceback.print_exc()
            print("Configuration file is incomplete")
            quit()
        else:
            globals()[var] = default
            val = eval(globals()[var])
            if isinstance(val, str):
                val = '"'+val+'"'
            print(var,"=", val, "# (default)")


def log_to_file(filename="qrhei.log"):
    """Set logging to file
    
    """
    manager = qr.Manager().log_conf
    #manager.log_on_screen = False
    manager.log_to_file = True
    manager.log_file_name = filename