# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 13:50:05 2017

@author: Johan
"""

# -*- coding: utf-8 -*-

from ..core.units import cm2int
from .molecularmodel import MolecularModel
from ..builders import pdb
from ..utils.vectors import normalize2

class ChlorophyllA(MolecularModel):
    
    def __init__(self, model_type=None, dp_length=4.582):
        super().__init__(model_type=model_type)
        
        self.pdbname = "CLA"
       
        self.set_default_energies([0.0, 15200.0*cm2int])
        # These values are taken from Muh, Lindorfer, et al. 
        # Physical Chemistry Chemical Physics, 2014. Chla: 4.58, Chlb 3.83:
        
        self.set_default_dipole_length((0,1), dp_length)
        
    def transition_dipole(self, transition=(0,1), data_type=None, data=None):
        """ Returns transition dipole moment vector
        
        """
        
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
                #print(k1,k2)
                raise Exception("No unique direction of"
                                +" a molecule's dipole found")
        else:
            raise Exception("Unknown data type")

        return d   
                                
        
    def position_of_center(self, data_type=None, data=None):
        """ Returns the position of the molecular center 
        
        """
        
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
                #print(k1, k2, k3, k4)
                raise Exception("No unique possition of a molecule found")
        else:
            raise Exception("Unknown data type")

        return pos   
        
        
    def pi_conjugated_system(self, data_type=None, data=None):
        """Returns the atoms and atom types in the pi-conjugated system
        
        Calculates and returns positions of all atoms in the pi-conjugated
        system of the molecule and the types of the atoms. 
        
        Parameters
        ----------
        
        data_type : string
            Type of the data; can be e.g. PDB
            
        data : 
            Data corresponding to the data type
        
        """
        data_type = self._check_data_type(data_type)
        
        if data_type == "PDB":
            
            pass

        else:
            raise Exception("Unknown data type")
        
        
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
        
        
class ChlorophyllB(MolecularModel):
    
    def __init__(self, model_type=None, dp_length=3.834):
        super().__init__(model_type=model_type)
        
        self.pdbname = "CHL"
        
        self.set_default_energies([0.0, 15700.0*cm2int])
        #These values are taken from Muh, Lindorfer, et al. Physical Chemistry Chemical Physics, 2014
        
        self.set_default_dipole_length((0,1), dp_length)
        
       
    
    def transition_dipole(self, transition=(0,1), data_type=None, data=None):
        """ Returns transition dipole moment vector
        
        """
        
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
                #print(k1,k2)
                raise Exception("No unique direction of"
                                +" a molecule's dipole found")
        else:
            raise Exception("Unknown data type")

        return d   
                                
        
    def position_of_center(self, data_type=None, data=None):
        """ Returns the position of the molecular center 
        
        """
        
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
                #print(k1, k2, k3, k4)
                raise Exception("No unique possition of a molecule found")
        else:
            raise Exception("Unknown data type")

        return pos   
        
        
    def pi_conjugated_system(self, data_type=None, data=None):
        """Returns the atoms and atom types in the pi-conjugated system
        
        Calculates and returns positions of all atoms in the pi-conjugated
        system of the molecule and the types of the atoms. 
        
        Parameters
        ----------
        
        data_type : string
            Type of the data; can be e.g. PDB
            
        data : 
            Data corresponding to the data type
        
        """
        data_type = self._check_data_type(data_type)
        
        if data_type == "PDB":
            
            pass

        else:
            raise Exception("Unknown data type")
        
        
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
            
