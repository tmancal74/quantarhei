# -*- coding: utf-8 -*-

import unittest
#import numpy

"""
*******************************************************************************


    Tests of the quantarhei.models.ModelGenerator class


*******************************************************************************
"""

#legacy = False
#import tempfile
from quantarhei.models.molecularmodel import MolecularModel
from quantarhei.models.chlorophylls import ChlorophyllA
from quantarhei.models.chlorophylls import ChlorophyllB
from quantarhei.models.bacteriopheophytins import BacterioPheophytin
from quantarhei.models.bacteriochlorophylls import BacterioChlorophyll

#from quantarhei import TimeAxis
        



class TestModelGenerator(unittest.TestCase):
    """(ModelGenerator) Tests for the model generator
    
    
    """
    
    def setUp(self):
        
        pass
        
        
        
        
    def testing_molecular_mode_init(self):
        """(MolecularModel) Testing initialization of the class
        
        """

        mmod = MolecularModel()
        
        mmod.set_default_energies([0.0, 1.0])
        mmod.set_default_dipole_length((0,1), 3.5)
        
        self.assertTrue(mmod.nstate==2)
        


    def testing_chlorophyllA(self):
        """(ChlorophyllA) Testing initialization of the class
        
        """    
        chlA = ChlorophyllA()
        
        self.assertTrue(chlA.pdbname=="CLA")
        
        
    def testing_chlorophyllB(self):
        """(ChlorophyllB) Testing initialization of the class
        
        """    
        chlB = ChlorophyllB()
        
        self.assertTrue(chlB.pdbname=="CHL")        
        
    def testing_bacterioPheo(self):
        """(BacterioPheophytin) Testing initialization of the class
        
        """    
        mol = BacterioPheophytin()
        
        self.assertTrue(mol.pdbname=="BPH") 
        
    def testing_bacterioChlo(self):
        """(BacterioChlorophyll) Testing initialization of the class
        
        """    
        mol = BacterioChlorophyll()
        
        self.assertTrue(mol.pdbname=="BCL") 

        
        