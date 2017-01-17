# -*- coding: utf-8 -*-
import unittest

"""
*******************************************************************************


    Tests of the quantarhei.spectroscopy.abs package


*******************************************************************************
"""
import numpy

from quantarhei.spectroscopy.abs import AbsSpectrumBase, AbsSpectrumDifference
from quantarhei import FrequencyAxis
from quantarhei import energy_units

class TestAbs(unittest.TestCase):
    """Tests for the abs package
    
    
    """
    
    def setUp(self,verbose=False):
        
        abss = AbsSpectrumBase()
        with energy_units("1/cm"):
            f = FrequencyAxis(10000.0,2000, 1.0)
            a = self._spectral_shape(f, [1.0, 11000.0, 100.0])
            abss.axis = f
            abss.data = a
        self.abss = abss
        self.axis = abss.axis
        
    def _spectral_shape(self, f, par):
        a = par[0]
        fd = par[1]
        sig = par[2]
        with energy_units("1/cm"):
            return a*numpy.exp(-((f.data-fd)/sig)**2)
        
    def _opt_spectral_shape(self, par):
        f = self.axis
        return self._spectral_shape(f, par)
        
            
    def test_abs_spectrum_base(self):
        """Testing basic function of AbsSpectrumBase class
        
        """
        pass


    def test_abs_difference_equal_zero(self):
        """Testing calculation of zero difference between abs spectra
        
        """
        
        #
        #  Test difference = 0
        # 
        par = [1.0, 11000.0, 100.0]
        with energy_units("1/cm"):
            f2 = self._spectral_shape(self.abss.axis, par) 
        
        abss2 = AbsSpectrumBase(axis=self.abss.axis, data=f2)
        ad = AbsSpectrumDifference(target=self.abss, optfce=abss2)
        d = ad.difference()
        
        self.assertAlmostEqual(d, 0.0)
        
    def test_abs_difference_formula(self):
        """Testing formula for calculation of difference between abs spectra
        
        """

        #
        #  Test difference formula
        #
        par = [1.0, 10900.0, 100.0]
        #with energy_units("1/cm"):
        f2 = self._spectral_shape(self.abss.axis, par) 
        
        abss2 = AbsSpectrumBase(axis=self.abss.axis, data=f2)
        ad = AbsSpectrumDifference(target=self.abss, optfce=abss2)
        d = ad.difference()            
        
        target = self.abss.data
        secabs = abss2.data
        x = self.abss.axis.data
        d2 = numpy.sum(numpy.abs((target-secabs)/
                                   (x[len(x)-1]-x[0])))
        
        self.assertAlmostEqual(d, d2)        
        
    def test_minimize(self):
        """Testing difference minimization using scipy.optimize package
        
        """
        ad = AbsSpectrumDifference(target=self.abss, 
                                   optfce=self._opt_spectral_shape)
        
        method = "Nelder-Mead"
        ini = [0.5, 10900, 80.0]
        p = ad.minimize(ini, method=method)
        
        numpy.testing.assert_array_almost_equal(p, [1.0, 11000.0, 100.0])

        
        