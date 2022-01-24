# -*- coding: utf-8 -*-

from ..core.units import cm2int
#from .molecularmodel import MolecularModel
#from ..builders import pdb
#from ..utils.vectors import normalize2
from .bacteriochlorophylls import BacterioChlorophyll

class BacterioPheophytin(BacterioChlorophyll):
    
    def __init__(self, model_type=None):
        super().__init__(model_type=model_type)
        
        self.pdbname = "BPH"
        
        self.default_energies[1] = 12500.0*cm2int
        self.default_dipole_lengths[0,1] = 5.8
        self.default_dipole_lengths[1,0] = 5.8
        

