# -*- coding: utf-8 -*-


import unittest
#import h5py
import tempfile
import numpy



"""
*******************************************************************************


    Tests of the quantarhei.core.time package


*******************************************************************************
"""

from quantarhei import TimeAxis

class TestTimeAxis(unittest.TestCase):
    """Tests for the units package
    
    
    """
    
    def test_of_time_axis_creation(self):
        """Testing TimeAxis creation """
        
        ta = TimeAxis(0.0, 1000, 0.1)
        
        self.assertEqual(ta.min,ta.data[0])
        self.assertEqual(ta.max,ta.data[ta.length-1])
        
        
    def test_if_time_axis_is_saveable(self):
        """Testing the Saveability of TimeAxis
        
        """
        
        #with h5py.File("test_file_ValueAxes",driver="core", 
        #                   backing_store=False) as f:   
        with tempfile.TemporaryFile() as f:
            ta = TimeAxis(0.0, 1000, 0.1)
        
            ta.save(f) #, test=True)
            f.seek(0)
            tb = TimeAxis()
            tb = tb.load(f) #, test=True)
            
        
        numpy.testing.assert_array_equal(ta.data,tb.data)
        
        
            
        
        