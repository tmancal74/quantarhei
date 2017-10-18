# -*- coding: utf-8 -*-

import unittest
import numpy

from quantarhei import DFunction, FrequencyAxis


class TestDFunction(unittest.TestCase):
    """Tests for the units package
    
    
    """
    
    def test_dfunction_saveable(self):
        """Testing if DFunction is saveable """
        
        wa = FrequencyAxis(0.0, 1000, 0.1)
        
        fw = numpy.exp(-wa.data)
        
        f = DFunction(wa,fw)
        
        f.plot()
        
        f.save("dfunction.hdf5")
        
        f2 = DFunction()
        f2.load("dfunction.hdf5")
        
        f2.plot()
        
        #self.assertEqual(wa.min,wa.data[0])
        #self.assertEqual(wa.max,wa.data[len(wa.data)-1])