# -*- coding: utf-8 -*-

import unittest

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
        
    def get_value(self):
        return self.convert_energy_2_current_u(self.value)
    

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
        self.assertEqual("1/fs",self.u.unit_repr("energy"))

