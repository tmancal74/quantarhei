# -*- coding: utf-8 -*-

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
        
            
