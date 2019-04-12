# -*- coding: utf-8 -*-
"""
    Module saveable

    Defines the class Saveable. Saveable objects can save themselves to and
    load themselves from a file.
    

"""

import os
from tempfile import TemporaryDirectory
import copy

from .parcel import Parcel
from .parcel import load_parcel


class Saveable:
    """Defines a class of objects that can save and load themselves
    
    
    """
    
    def save(self, filename, comment=None, test=False):
        """Saves the object with all its content into a file
        
        
        Parameters
        ----------
        
        filename : str or File
            Name of the file or a file object to which the content of
            the object will be saved 
        
        comment : str
            A comment which will be saved together with the content of
            the object
            
        test : bool
            If test is True, and file descriptor is submitted, we save into
            the file and move to its start for subsequent reading
        
        
        """
        p = Parcel()
        p.set_content(self)
        p.set_comment(comment)
        
        p.save(filename)
        
        if test:
            if not isinstance(filename, str):
                filename.seek(0)


    def load(self, filename, test=False):
        """Loads an object from a file and returns it
        
        Parameters
        ----------
        
        filename : str or File
            Filename of the file or file descriptor of the file from which
            and object should be loaded.
        
        """
        if test:
            if not isinstance(filename, str):
                filename.seek(0)
        
        return load_parcel(filename)
        
        
    def scopy(self):
        """Creates a copy of the object by saving and loading it
        
        """
        
        with TemporaryDirectory() as td:
            fname = os.path.join(td, "ssave.qrp")
            self.save(fname)
            
            no = load_parcel(fname)
            
        return no
    
    
    def deepcopy(self):
        """Returns a deep copy of the self
        
        
        """
        return copy.deepcopy(self)
    
    
    def copy(self):
        """Returns a shallow copy of the self
        
        
        """
        return copy.copy(self)