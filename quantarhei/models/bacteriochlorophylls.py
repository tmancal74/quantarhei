# -*- coding: utf-8 -*-
from ..core.ubits import cm2int
from .molecularmodel import MolecularModel

class BacterioChlorophyll(MolecularModel):
    
    def __init__(self, model_type=None):
        super().__init__(model_type=model_type)
        
        self.pdbname = "BCL"
        
        self.default_energies[1] = 12500.0*cm2int
        self.default_dipole_lengths[0,1] = 12.0
        self.default_dipole_lengths[1,0] = 12.0
        
        
    def dipole_direction(self, transition=(0,1), data_type=None, data):
        
        data_type = self._check_data_type(data_type)
        
        if data_type == "PDB":
            pass 
        
        
    def position_of_center(self, data_type=None, data):
        
        data_type = self._check_data_type(data_type)
        
        if data_type == "PDB":
            pass        
        
        
    def _check_data_type(self, data_type):
        """If non data_type is specified, the default is taken (if known)
        
        """
        if data_type is None:
            if self.model_type is None:
                raise Exception()
            else:
                return self.model_type
        else:
            return data_type
            
