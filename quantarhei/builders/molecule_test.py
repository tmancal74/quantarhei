# -*- coding: utf-8 -*-
"""
    Class to support tests on the Molecule class




    Class Details
    -------------    

"""
from .molecules import Molecule
from .modes import Mode

class TestMolecule(Molecule):
    """Class to support tests on the Molecule class
    
    
    Parameters
    ----------
    
    name : str
        Name of the test molecule
    
    
    Examples
    --------
    
    The TestMolecule must be created with a name. The name specifies the
    content of the TestMolecule object.
    
    >>> mol = TestMolecule()
    Traceback (most recent call last):
        ...
    Exception: TestMolecule name not specified
    
    
    Two-level molecule with one vibrational mode
    
    >>> mol = TestMolecule(name="two-levels-1-mode")
    >>> mol.get_number_of_modes()
    1
    
    """
    
    def __init__(self, name=None):

    
        if name is None:
            raise Exception("TestMolecule name not specified")


        if name == "two-levels-1-mode":
            
            super().__init__([0.0, 1.0], name=name)
            
            mod = Mode()
            self.add_Mode(mod)
            
        else:
            
            raise Exception("Unknown TestMolecule")