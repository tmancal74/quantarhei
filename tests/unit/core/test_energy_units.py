# -*- coding: utf-8 -*-

import unittest

"""
*******************************************************************************


    Tests of the quantarhei.energy_units class


*******************************************************************************
"""  

from quantarhei import energy_units
from quantarhei.core.units import cm2int
from quantarhei import set_current_units

# let us reuse a class from previous test
from .test_UnitsManaged import UnitsManagedObject        

class TestEnergyUnits(unittest.TestCase):
    
    def setUp(self):
        # reusining class from previous test
        self.u = UnitsManagedObject()
        
    def test_using_different_units(self):
        """Testing that 'energy_units' context manager works
        
        """
        # set value in 1/cm
        with energy_units("1/cm"):
            val = 100.0
            val2 = 50.0
            self.u.set_value(val)
            
            w = UnitsManagedObject()
            w.set_value(val2)
            
        # compare it in internal units
        self.assertEqual(self.u.value,val*cm2int)
        self.assertEqual(w.value,val2*cm2int)
        
    def test_using_global_units_testing(self):
        """Testing that 'set_current_units' works
        
        """
        a = UnitsManagedObject()
        a.set_value(1.0)
        
        #
        # Test of energy units
        #
        units = dict(energy="1/cm")
        
        set_current_units(units)
        
        self.assertEqual(a.get_value(),1.0/cm2int)
        
        val1 = 100
        val2 = 50.0
        w = UnitsManagedObject()
        w.set_value(val1)
        self.u.set_value(val2)
        
        self.assertEqual(self.u.get_value(), val2)
        self.assertEqual(w.get_value(), val1)
        
        set_current_units()
        
        self.assertEqual(self.u.value, val2*cm2int)
        self.assertEqual(w.value, val1*cm2int)
        
        self.assertEqual(self.u.get_value(), val2*cm2int)
        self.assertEqual(w.get_value(), val1*cm2int)
        
        #
        # Test of frequency units
        #
        units = dict(frequency="1/cm")        
        
        set_current_units(units)
        
        self.assertEqual(a.get_value(),1.0/cm2int)
        
        val1 = 100
        val2 = 50.0
        w = UnitsManagedObject()
        w.set_value(val1)
        self.u.set_value(val2)
        
        self.assertEqual(self.u.get_value(), val2)
        self.assertEqual(w.get_value(), val1)
        
        set_current_units()
        
        self.assertEqual(self.u.value, val2*cm2int)
        self.assertEqual(w.value, val1*cm2int)
        
        self.assertEqual(self.u.get_value(), val2*cm2int)
        self.assertEqual(w.get_value(), val1*cm2int)
          
        
