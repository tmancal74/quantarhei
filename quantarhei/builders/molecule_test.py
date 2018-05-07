# -*- coding: utf-8 -*-

from .molecules import Molecule

class TestMolecule(Molecule):
    
    def __init__(self, name=None):
        """Class to support tests on the Molecule class
        
        
        Parameters
        ----------
        
        name : str
            Name of the test molecule
            
        Names
        -----
        
        
        """
    
        if name is None:
            raise Exception("TestMolecule name not specified")


