# -*- coding: utf-8 -*-
import unittest


"""
*******************************************************************************


    Tests of the quantarhei.Manager class


*******************************************************************************
"""

from quantarhei import Manager

class TestManager(unittest.TestCase):
    """Tests for the Manager class
    
    
    """
    
    def setUp(self):
        pass
     
    
    def test_that_Manager_is_a_singleton(self):
        """Testing that Manager object is a singleton
        
        
        """
        m = Manager()
        n = Manager()
    
        if m is not n:
            raise Exception()
        
        
"""
*******************************************************************************


    Tests of the quantarhei.core.UnitsManaged class


*******************************************************************************
"""        
        
from quantarhei.core import UnitsManaged
    
    
class UnitsManagedObject(UnitsManaged):
    """Test object

    This class should have methods for conversion of units

    """
    def __init__(self):
        self.value = 0.0
        
    def set_value(self,val):
        self.value = self.convert_energy_2_internal_u(val)
    

class TestUnitsManaged(unittest.TestCase):
    
    def setUp(self):
        self.u = UnitsManagedObject()
        
    def test_trivial_units_conversion(self):
        """Test of the inheritance from the UnitsManaged class
        
        """       
        val = 3.141568
        self.u.set_value(val)
        
        self.assertEqual(val,self.u.value)
        
    def test_internal_units_representation(self):
        """Testing that default internal units are used
        
        """
        self.assertEqual("2pi/fs",self.u.unit_repr("energy"))

    
"""
*******************************************************************************


    Tests of the quantarhei.energy_units class


*******************************************************************************
"""  

from quantarhei import energy_units
from quantarhei.core.units import cm2int


class TestEnergyUnits(unittest.TestCase):
    
    def setUp(self):
        # Using a class defined for the previous test
        self.u = UnitsManagedObject()
        
    def test_using_different_units(self):
        """Testing that 'energy_units' context manager works
        
        """
        # set value in 1/cm
        with energy_units("1/cm"):
            val = 100.0
            self.u.set_value(val)
            
        # compare it in internal units
        self.assertEqual(self.u.value,val*cm2int)
        
        
        