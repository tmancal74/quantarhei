# -*- coding: utf-8 -*-
import numpy

class MolecularModel:
    """Representation of a non-Quantarhei molecule.
    
    This class handles conversion of non-Quantarhei models of molecules
    
    """
    
    def __init__(self, model_type=None):
        
        self.model_type = model_type
        self.nstate = 2
        self.default_energies = numpy.array([0.0, 0.0], dtype=numpy.float64)
        self.default_dipole_lengths = numpy.zeros((self.nstate,self.nstate),
                                                  dtype=numpy.float64)
        
        