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
        
        self.assertTrue(mmod.nstate==2)
        
