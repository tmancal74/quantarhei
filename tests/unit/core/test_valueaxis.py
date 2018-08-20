# -*- coding: utf-8 -*-


import unittest
import numpy
#import h5py
import tempfile


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
        
        #with h5py.File("test_file_ValueAxes",driver="core", 
        #                   backing_store=False) as f:
        with tempfile.TemporaryFile() as f:
        
            ta = ValueAxis(0.0, 1000, 0.1)
        
            ta.save(f) #, test=True)
            f.seek(0)
        
            tb = ValueAxis()
            tb = tb.load(f) #, test=True)
        
        numpy.testing.assert_array_equal(ta.data,tb.data)
        self.assertEqual(ta.length,tb.length)
        self.assertEqual(ta.start,tb.start)
        self.assertEqual(ta.step,tb.step)
