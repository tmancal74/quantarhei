# -*- coding: utf-8 -*-

class LogConf:
    """Holds information of the configuration for Quantarhei logging system
    
    When initialized, it holds hard default configuration of the logging
    system.
        
    """   
    
    def __init__(self):

        from queue import LifoQueue
        
        self.verbosity = 5
        self.fverbosity = 5
        self.verbose = True
        self.log_on_screen = True
        self.log_indent = 0
        self.log_to_file = False
        self.log_file_opened = False
        self.log_file_name = "qrhei.log"
        self.log_file_appendix = ""
        self.log_file = None
        self.initialized = False
        self.time_stamp = LifoQueue()        

    def __del__(self):
        """Destructor to be used on garbage collection
        
        It closes openned log_file if it was not closed by something before
        
        """
        if self.log_file_opened:
            self.log_file.close()
            