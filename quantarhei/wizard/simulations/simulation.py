# -*- coding: utf-8 -*-

import datetime

from ...core.managers import Manager
from ...core.saveable import Saveable

class Simulation(Saveable):
    """Quantarhei simulation class
    
    
    """

    VERBOSE = 10
    WARNING = 5
    ERROR   = 1
    
    def __init__(self, loglevel=0):
        self._loglevel = loglevel
        self._indent_width = 4
        self._set_indent_string()
        self.verbosity = self.ERROR
        self.on_screen = True
        self.into_file = False
        self._indent = ""
        self._starline = "\n************************************************\n"
        
    def setup(self):
        pass
        
    def run(self):
    
        # greeting
        self._open_logfile()
        self._print_greetings()
        
        
        # setup evaluation
        self._set_indent_level(0)
        self._printlog("Evaluating setup ...", loglevel=0) 
        self._incr_indent_level()
        
        self._evaluate_setup()
        
        self._decr_indent_level()
        self._printlog("...done", loglevel=0)

        
        # Building objects for the simulation
        self._printlog("\nBuilding objects ...")
        self._incr_indent_level()
        
        self._build() 
        
        self._decr_indent_level()
        self._printlog("...done")
        
        
        # Simulation itself
        self._printlog("\nRunning simulation")
        self._incr_indent_level()

        self._implementation()

        self._decr_indent_level()
        self._printlog("...done")
        
        
        # final wrap-up and clean-up
        self._print_goodbye()
        self._close_logfile()
        
        
    def _open_logfile(self):
        pass


    def _close_logfile(self):
        pass


    def _get_timestamp(self, filename=False):
        """Returns current time stemp
        
        """
        if filename:
            return '{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now())
        return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())

    
    def _print_greetings(self):
        """Prints opening greeting of the simulation
        
        """
        time_stamp = self._get_timestamp()

        grstring = self._starline + "* Quantarhei Simulation\n*"
        grstring += "\n* Class name: "+self.__class__.__name__
        grstring += "\n*"
        grstring += "\n* Quantarhei version "+Manager().version
        grstring += "\n* Initial timestamp: "+time_stamp
        grstring += self._starline
        self._printlog(grstring, loglevel=0)

       
    def _print_goodbye(self):
        """Prints the last message before leaving
        
        """
        time_stamp = self._get_timestamp()

        grstring = self._starline + "* Simulation finished\n*"
        grstring += "\n* Final timestamp: "+time_stamp
        grstring += self._starline
        
        self._printlog(grstring, loglevel=0)

        
    def _printlog(self, *args, loglevel=0):
        """Logs output on screen and into a file
        
        """
        # define loglevel
        if loglevel < self.verbosity:
      
            if self.on_screen:
                print(self._indent, *args)
            
            if self.into_file:
                print(self._indent, args, file=self._file)
            
        
    def _set_indent_string(self, indent_character=" "):
        """Sets the form of indent
            
        """
        self._indent_string = ""
        for i in range(self._indent_width):
            self._indent_string += indent_character


    def _set_indent_level(self, level):
        """Sets the indent level and creates the indent
        
        """
        self._indent_level = level
        self._indent = ""
        for i in range(self._indent_level):
            self._indent += self._indent_string
            
    def _incr_indent_level(self):
        self._set_indent_level(self._indent_level + 1)
        
    def _decr_indent_level(self):
        self._set_indent_level(self._indent_level - 1)
            

    def _evaluate_setup(self):
        
        self.setup()
        
        
    def _build(self):
        pass
    
    
    def _implementation(self):
        pass
        
            
    # Print iterations progress
    def _printProgressBar(self, iteration, total, 
                          prefix = '', suffix = '', 
                          decimals = 1, length = 100,
                          fill='*'):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
            
        Based on: 
        https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
        """
#                          fill = '█'):
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
        # Print New Line on Complete
        if iteration == total: 
            print()        
            

