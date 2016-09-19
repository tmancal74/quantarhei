# -*- coding: utf-8 -*-

import unittest



"""
*******************************************************************************


    Tests of the quantarhei.core.units package


*******************************************************************************
"""

from quantarhei import FrequencyAxis

class TestFrequencyAxis(unittest.TestCase):
    """Tests for the units package
    
    
    """
    
    def test_of_frequency_axis_creation(self):
        """Testing FrequencyAxis creation """
        
        wa = FrequencyAxis(0.0, 1000, 0.1)
        
        self.assertEqual(wa.omin,wa.frequency[0])
        self.assertEqual(wa.omax,wa.frequency[len(wa.data)-1])
        