# -*- coding: utf-8 -*-


import unittest
import numpy


"""
*******************************************************************************


    Tests of the quantarhei.core.valueaxis package


*******************************************************************************
"""

from quantarhei.core.valueaxis import ValueAxis

class TestValueAxis(unittest.TestCase):
    """Tests for the units package
    
    
    """
    
    def test_of_value_axis_creation(self):
        """Testing ValueAxis creation """
        
        ta = ValueAxis(0.0, 1000, 0.1)
        
        self.assertEqual(ta.min,ta.data[0])
        self.assertEqual(ta.max,ta.data[ta.length-1])
        
        
    def test_if_value_axis_is_saveable(self):
        """Testing the Saveability of ValueAxis
        
        """
        
        ta = ValueAxis(0.0, 1000, 0.1)
        
        ta.save("value.hdf5")
        
        tb = ValueAxis()
        tb.load("value.hdf5")
        
        numpy.testing.assert_arrays_equal(ta.data,tb.data)

