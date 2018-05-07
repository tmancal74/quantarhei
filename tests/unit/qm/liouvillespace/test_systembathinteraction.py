# -*- coding: utf-8 -*-

import unittest
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.qm.SystemBathInteraction class


*******************************************************************************
"""

from quantarhei.qm.liouvillespace.systembathinteraction_test \
import TestSystemBathInteraction

class SBITest(unittest.TestCase):
    """Tests for the SystemBathInteraction class
    
    
    """
    
    def setUp(self,verbose=False):
        
        self.verbose = verbose
        
        
    def test_system_bath_interaction_creation(self):
        """(SystemBathInteraction) Testing creation
        
        """
        
        sbi = TestSystemBathInteraction(name="dimer-2-env")
        
        self.assertEqual(sbi.sbitype, "Linear_Coupling")
        
        
    def test_get_temperature_and_has_temperature(self):
        """(SystemBathInteraction) Testing has_temperature() and get_temperature()
        
        """
        
        sbi = TestSystemBathInteraction(name="dimer-2-env")
        
        self.assertTrue(sbi.has_temperature())
        
        T = sbi.get_temperature()
        self.assertAlmostEqual(300.0, T)
        
        
        sbi = TestSystemBathInteraction(name="dimer-2-lind")
        self.assertFalse(sbi.has_temperature())
        
        
    
        