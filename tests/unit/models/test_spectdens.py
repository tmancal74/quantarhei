# -*- coding: utf-8 -*-

import unittest
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.models.ModelGenerator class


*******************************************************************************
"""

#legacy = False
#import tempfile
from quantarhei.models.spectral_densities.wendling_2000 import wendling_2000a
from quantarhei.models.spectral_densities.renger_2002 import renger_2002a
from quantarhei import FrequencyAxis
from quantarhei import energy_units
        



class TestModelSpectralDensities(unittest.TestCase):
    """(ModelGenerator) Tests for the model generator
    
    
    """
    
    def setUp(self):
        
        with energy_units("1/cm"):
            self.fa = FrequencyAxis(0.0, 1000, 1.0)
        self.dbW = wendling_2000a()
        self.dbR = renger_2002a()
        
        
        
    def testing_wendling_2000(self):
        """(LITERATURE SPECT. DENSITIES) Spectral density of Wendling at al.
        
        """

        sp = self.dbW.get_SpectralDensity(self.fa)
        
        with energy_units("1/cm"):
            val = sp.at(100.0)
        
            print(val)
        
        self.assertAlmostEqual(val,1.52848396142e-08)
        
        
    def testing_renger_2002(self):
        """(LITERATURE SPECT. DENSITIES) Spectral density of Renger at al.
        
        """

        sp = self.dbR.get_SpectralDensity(self.fa)
        
        with energy_units("1/cm"):
            val = sp.at(100.0)
        
            print(val)
        
        self.assertAlmostEqual(val,0.0228204533609)        