# -*- coding: utf-8 -*-

import unittest
import numpy
#import h5py
import tempfile

from quantarhei import energy_units

from quantarhei import TimeAxis


"""
*******************************************************************************


    Tests of the quantarhei.core.units package


*******************************************************************************
"""

from quantarhei import FrequencyAxis

class TestFrequencyAxis(unittest.TestCase):
    """Tests for the units package
    
    
    """
    
    def test_of_frequency_axis_creation(self):
        """Testing FrequencyAxis creation """
        
        wa = FrequencyAxis(0.0, 1000, 0.1)
        
        self.assertEqual(wa.min,wa.data[0])
        self.assertEqual(wa.max,wa.data[len(wa.data)-1])
        
        
    def test_(self):
        """Testing FrequencyAxis values """
        
        ta = TimeAxis(0.0, 1000, 0.1)
        
        ta.atype='complete'
        wa = ta.get_FrequencyAxis()
        
        # in internal units
        wadata =  2.0*numpy.pi*\
        numpy.fft.fftshift(numpy.fft.fftfreq(ta.length,d=ta.step))
        
        numpy.testing.assert_allclose(wa.data, wadata, atol=1.0e-12)
        
        
    def test_of_fa_location_in_different_units(self):
        """Testing location of value in different units """
        ta = TimeAxis(0.0, 1000, 0.1)
        
        ta.atype='complete'
        wa = ta.get_FrequencyAxis()
        
        with energy_units("1/cm"):

            val = 10000.0
            nsni = numpy.int(numpy.floor((val-wa.start)/wa.step))

            i1 = wa.locate(10000.0)
            self.assertEqual(nsni, i1[0])
            
        #with h5py.File("test_file_ValueAxes",driver="core", 
        #                   backing_store=False) as f:    
        with tempfile.TemporaryFile() as f:
        
            wa.save(f, test=True)
            
            tb = FrequencyAxis()
            tb = tb.load(f, test=True)
            
            
        
        numpy.testing.assert_array_equal(wa.data,tb.data)
            #
            
            
            
            
    