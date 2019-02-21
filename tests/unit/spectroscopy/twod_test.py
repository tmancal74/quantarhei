# -*- coding: utf-8 -*-
import unittest

import tempfile
import os

import numpy

from nose.tools import assert_raises

import quantarhei as qr
from quantarhei.spectroscopy.twod2 import TwoDSpectrumBase

"""
*******************************************************************************


    Tests of the quantarhei.spectroscopy.twod package


*******************************************************************************
"""

import quantarhei as qr


class TestTwod(unittest.TestCase):
    """Tests for the twod package
    
    
    """
    
    def setUp(self,verbose=False):
        pass
        

        
    def test_TwoDSpectrumBase(self):
        """Testing basic functions of the TwoDSpectrumBase class
        
        """
        t1 = qr.TimeAxis(0.0, 1000, 1.0)
        t3 = qr.TimeAxis(0.0, 1000, 1.0)
        
        twodB = TwoDSpectrumBase()
        twodB.set_axis_1(t1)
        twodB.set_axis_3(t3)
        data1 = numpy.zeros((t1.length, t3.length), dtype=qr.COMPLEX)
        data2 = numpy.zeros((t1.length, t3.length), dtype=qr.COMPLEX)
        
        data1[:,:] = 1.1
        
        twodB._add_data(data1, resolution="off", dtype="total")

        numpy.testing.assert_equal(twodB.data, data1)
        
        data2[:,:] = 2.2
        
        with assert_raises(Exception):      
            twodB.data = data2

        twodB._allow_data_writing = True
        twodB.data = data2

        numpy.testing.assert_equal(twodB.data, data2)
        
        with tempfile.TemporaryDirectory() as tdir:
            
            names = ["data.mat", "data.npy"]
            for name in names:
                fname = os.path.join(tdir,name)
                twodB.save_data(fname)
                
                twodB2 = TwoDSpectrumBase()
                twodB2.set_axis_1(t1)
                twodB2.set_axis_3(t3)
    
                twodB2.load_data(fname)
                
                numpy.testing.assert_equal(twodB.data, twodB2.data)


    
    def test_TwoDSpectrumCalculator(self):
        """Testing basic functions of the TwoDSpectrumCalculator class
        
        """
        t1 = qr.TimeAxis(0.0, 1000, 1.0)
        t3 = qr.TimeAxis(0.0, 1000, 1.0)
        
        t2 = qr.TimeAxis(30, 10, 10.0)
        
        twod_calc = qr.TwoDSpectrumCalculator(t1, t2, t3)
        
        
            
    