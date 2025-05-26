# -*- coding: utf-8 -*-
import numpy

from ..core.managers import EnergyUnitsManaged

class MolecularModel(EnergyUnitsManaged):
    """Representation of a non-Quantarhei molecule.
    
    This class handles conversion of non-Quantarhei models of molecules
    
    """
    
    def __init__(self, model_type=None):
        
        self.model_type = model_type
        self.nstate = 2
        self.default_energies = numpy.array([0.0, 0.0], dtype=numpy.float64)
        self.default_dipole_lengths = numpy.zeros((self.nstate,self.nstate),         
                                                  dtype=numpy.float64)

        
    def set_default_energies(self, elenergies):
        k = 0
        for en in elenergies:
            self.default_energies[k] = self.convert_2_internal_u(en)
            k += 1
            
        
    def set_default_dipole_length(self,transition, val):
        self.default_dipole_lengths[transition[0],transition[1]] = val
        self.default_dipole_lengths[transition[1],transition[0]] = val        