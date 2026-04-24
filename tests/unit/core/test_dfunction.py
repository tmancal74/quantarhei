# -*- coding: utf-8 -*-

import unittest
import numpy
#import h5py
import tempfile

from quantarhei import DFunction, FrequencyAxis


class TestDFunction(unittest.TestCase):
    """Tests for the units package
    
    
    """
    
    def test_dfunction_saveable(self):
        """Testing if DFunction is saveable """
        
        wa = FrequencyAxis(0.0, 1000, 0.1)
        
        fw = numpy.exp(-wa.data)
        
        fce = DFunction(wa,fw)
        
        #fce.plot()

        #with h5py.File("test_file_1",driver="core", 
        #                   backing_store=False) as f:
        with tempfile.TemporaryFile() as f:
        
            fce.save(f, test=True)
        
            fce2 = DFunction()
            fce2 = fce2.load(f, test=True)
        
        #fce2.plot()
        
        numpy.testing.assert_array_equal(fce.data, fce2.data)
        
        
    def test_dfunction_saveable_2(self):
        """Testing if DFunction can save spline initiated object correctly
        
        """
        
        wa = FrequencyAxis(0.0, 100, 1.0)
        
        fw = numpy.exp(-wa.data/30.0)
        
        fce = DFunction(wa,fw)
        
        val_mez_linear = fce.at(20.5)
        val_mez_spline = fce.at(20.5, approx="spline")
        
#        print("init", fce._splines_initialized)
#        print(val_mez_linear)
#        print(val_mez_spline)

        #with h5py.File("test_file_1",driver="core", 
        #                   backing_store=False) as f:
        with tempfile.TemporaryFile() as f:
        
            fce.save(f, test=True)
        
            fce2 = DFunction()
            fce2 = fce2.load(f, test=True)

        self.assertEqual(fce._splines_initialized, True)
        self.assertEqual(fce2._splines_initialized, True)
            
        new_mez_linear = fce2.at(20.5)
        new_mez_spline = fce2.at(20.5, approx="spline")
#        print(new_mez_linear)
#        print(new_mez_spline)

        self.assertEqual(val_mez_spline, new_mez_linear)
        self.assertEqual(val_mez_spline, new_mez_spline)
        self.assertEqual(fce._splines_initialized, fce2._splines_initialized)
        numpy.testing.assert_array_equal(fce.data, fce2.data)
