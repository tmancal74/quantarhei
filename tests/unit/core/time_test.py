# -*- coding: utf-8 -*-


import unittest



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
        
        self.assertEqual(ta.tmin,ta.time[0])
        self.assertEqual(ta.tmax,ta.time[len(ta.data)-1])
        
        