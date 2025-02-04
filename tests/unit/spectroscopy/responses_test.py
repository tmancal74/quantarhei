# -*- coding: utf-8 -*-
import unittest

import tempfile
import os

import numpy

from nose.tools import assert_raises

import quantarhei as qr

from quantarhei.utils.vectors import X 

#import matplotlib.pyplot as plt


"""
*******************************************************************************


    Tests of the quantarhei.spectroscopy.responses package


*******************************************************************************
"""


class TestResponses(unittest.TestCase):
    """Tests for the response package
    
    
    """
    
    def setUp(self,verbose=False):
        pass

        
    def test_Responses(self):
        """Testing basic functions of the TwoDSpectrumBase class
        
        """
        t1 = qr.TimeAxis(0.0, 1000, 1.0)
        t3 = qr.TimeAxis(0.0, 1000, 1.0)
        
        resp = qr.ResponseFunction("R2g")
        
        
        def gfce(t):
            
            return 1.0*(t**2)
        
        
        def efce(t2, t1s, t3s, gg):
            
            N1 = t1s.shape[0]
            N3 = t3s.shape[0]
            ret = numpy.zeros((N1,N3), dtype=complex)
            it1 = 0
            for t1 in t1s:
                it2 = 0
                for t3 in t3s:
                    ret[it1,it2] = numpy.exp(-gg(t1) -gg(t3))
                    it2 += 1
                it2 += 1
        
            return ret
        
        
        resp.set_evaluation_function(efce)
        
        
        d1 = [1.0, 0.0, 0.0]
        ome1 = 10100.0
        ome3 = 10100.0
        rwa = 1.0 
        
        # Fixme: set correct units
        with qr.energy_units("1/cm"):
            resp.set_frequencies(ome1, ome3)
            resp.set_rwa(rwa)          
        
        resp.set_dipoles(d1)
        
        args = (gfce,)
        resp.set_auxliary_arguments(args)
        
        
        # FIXME: create a laboratory object
        lab = qr.LabSetup()
        lab.set_pulse_polarizations(pulse_polarizations=(X,X,X), detection_polarization=X)
        
        # FIXME: create a reasonable molecule with everything needed
        sys = qr.Molecule([0.0, 1.0])
        sys.set_dipole(0,1,[1.0, 0.8, 0.8])

        
        mx = resp.calculate_matrix(lab, sys, 10, t1.data, t3.data, rwa)


        #print("shape =", mx.shape)
        
        #with assert_raises(Exception):      
        #    twodB.data = data2


        #numpy.testing.assert_equal(twodB.data, data2)
        
        #with tempfile.TemporaryDirectory() as tdir:
        
        #    numpy.testing.assert_equal(twodB.data, twodB2.data)


    
    def test_LouvillePathway(self):
        """Testing basic functions of the TwoDSpectrumCalculator class
        
        """
        t1 = qr.TimeAxis(0.0, 1000, 1.0)
        t3 = qr.TimeAxis(0.0, 1000, 1.0)
        
        t2 = qr.TimeAxis(30, 10, 10.0)
        
        #twod_calc = qr.TwoDResponseCalculator(t1, t2, t3)
        
        
            
    
