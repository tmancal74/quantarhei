# -*- coding: utf-8 -*-

import numpy

from ..core.units import cm2int
from .molecularmodel import MolecularModel
from ..builders import pdb
from ..utils.vectors import normalize2

class BacterioChlorophyll(MolecularModel):
    
    def __init__(self, model_type=None):
        super().__init__(model_type=model_type)
        
        self.pdbname = "BCL"
        
        self.default_energies[1] = 12500.0*cm2int
        self.default_dipole_lengths[0,1] = 5.8
        self.default_dipole_lengths[1,0] = 5.8
        
        
    def transition_dipole(self, transition=(0,1), data_type=None, data=None):
        
        data_type = self._check_data_type(data_type)
        
        if data_type == "PDB":
            k1 = 0
            k2 = 0
            for line in data:
                if pdb.line_matches(line, by_atmName="ND"):
                    xyz1 = pdb.line_xyz(line)
                    k1 += 1
                if pdb.line_matches(line, by_atmName="NB"):
                    xyz2 = pdb.line_xyz(line)
                    k2 += 1
            # FIXME: what to do with alternate locations???
            if (k1 >= 1) and (k2 >= 1):
                d = xyz1 - xyz2
                d = normalize2(d, norm=self.default_dipole_lengths[0,1])
            else:
                print(k1,k2)
                raise Exception("No unique direction of"
                                +" a molecule's dipole found")
        else:
            raise Exception("Unknown data type")

        return d   
                                
                                
        
        return normalize2(numpy.array([1.0, 1.0, 1.0]), 5.8)
        
    def position_of_center(self, data_type=None, data=None):
        
        data_type = self._check_data_type(data_type)
        
        if data_type == "PDB":
            k1 = 0
            k2 = 0
            k3 = 0
            k4 = 0
            for line in data:
                if pdb.line_matches(line, by_atmName="NA"):
                    xyz1 = pdb.line_xyz(line)
                    k1 += 1
                if pdb.line_matches(line, by_atmName="NB"):
                    xyz2 = pdb.line_xyz(line)
                    k2 += 1
                if pdb.line_matches(line, by_atmName="NC"):
                    xyz3 = pdb.line_xyz(line)
                    k3 += 1
                if pdb.line_matches(line, by_atmName="ND"):
                    xyz4 = pdb.line_xyz(line)
                    k4 += 1
            if (k1 >= 1) and (k2 >= 1) and (k3 >= 1) and (k4 >= 1):
                pos = (xyz1 + xyz2 + xyz3 + xyz4)/4.0
            else:
                print(k1, k2, k3, k4)
                raise Exception("No unique possition of a molecule found")
        else:
            raise Exception("Unknown data type")

        return pos   
        
        
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
            
