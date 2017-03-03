# -*- coding: utf-8 -*-
import unittest

"""
*******************************************************************************


    Tests of the quantarhei.spectroscopy.twod package


*******************************************************************************
"""
import numpy

from quantarhei.spectroscopy.twod2 import TwoDSpectrumContainer
from quantarhei.spectroscopy.twod2 import TwoDSpectrumCalculator

from quantarhei import TimeAxis


class TestTwod(unittest.TestCase):
    """Tests for the twod package
    
    
    """
    
    def setUp(self,verbose=False):
        pass
        

        
    def test_TwoDSpectrumBase(self):
        """Testing basic functions of the TwoDSpectrumBase class
        
        """
        pass
    
    
    def test_TwoDSpectrumCalculator(self):
        """Testing basic functions of the TwoDSpectrumCalculator class
        
        """
        t1 = TimeAxis(0.0, 1000, 1.0)
        t3 = TimeAxis(0.0, 1000, 1.0)
        
        t2 = TimeAxis(30, 10, 10.0)
        
        twod_calc = TwoDSpectrumCalculator(t1, t2, t3)
        
        
            
    