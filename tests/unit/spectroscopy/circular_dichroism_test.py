# -*- coding: utf-8 -*-
import unittest
import copy

"""
*******************************************************************************


    Tests of the quantarhei.spectroscopy.circular_dichroism package


*******************************************************************************
"""
import numpy
#import h5py

from quantarhei.spectroscopy.circular_dichroism import CircDichSpectrumBase, CircDichSpectrumContainer
from quantarhei import FrequencyAxis
from quantarhei import energy_units

from quantarhei import CircDichSpectrum, CircDichSpectrumCalculator
from quantarhei import Molecule, Aggregate, CorrelationFunction, TimeAxis

class TestCircDich(unittest.TestCase):
    """Tests for the circular dichroism package
    
    
    """
    
    def setUp(self,verbose=False):
        
        circdichs = CircDichSpectrumBase()
        with energy_units("1/cm"):
            f = FrequencyAxis(10000.0,2000, 1.0)
            a = self._spectral_shape(f, [1.0, 11000.0, 100.0])
            circdichs.axis = f
            circdichs.data = a
        self.circdichs = circdichs
        self.axis = circdichs.axis
        
        #make a single-molecule system
        time = TimeAxis(0.0,1000,1.0)
        self.ta = time
        with energy_units("1/cm"):
            mol1 = Molecule(elenergies=[0.0, 12000.0])
            mol2 = Molecule(elenergies=[0.0, 12000.0])

            params = dict(ftype="OverdampedBrownian", reorg=20, cortime=100,
                          T=300)
            mol1.set_dipole(0,1,[0.0, 1.0, 0.0])
            mol2.set_dipole(0,1,[0.0, 1.0, 0.0])

            cf = CorrelationFunction(time, params)
            mol1.set_transition_environment((0,1),cf)
            mol2.set_transition_environment((0,1),cf)
            
            mol1.position = [0, 0, 0]
            mol2.position = [10, 10, 10]
            
            agg = Aggregate(name="tester_dim", molecules=[mol1, mol2])
            agg.build()
            
            circdich_calc = CircDichSpectrumCalculator(time, system=agg)
            circdich_calc.bootstrap(rwa=12000)
            
        circdich1 = circdich_calc.calculate()
        
        self.circdich1 = circdich1
        
    def _spectral_shape(self, f, par):
        a = par[0]
        fd = par[1]
        sig = par[2]
        with energy_units("1/cm"):
            return a*numpy.exp(-((f.data-fd)/sig)**2)
        
    def _opt_spectral_shape(self, par):
        f = self.axis
        a = CircDichSpectrumBase(axis=f,data=self._spectral_shape(f, par))
        return a
    
    def test_CircDichSpectrumBase(self):
        """Testing the functions of the CircDichSpectrumBase class
        
        """
        def reset_tester():
            circdichs2 = CircDichSpectrumBase()
            circdichs2.set_axis(copy.deepcopy(self.circdichs.axis))
            circdichs2.set_data(copy.deepcopy(self.circdichs.data))
            return circdichs2
        
        def gaussian(x, height, center, fwhm, offset=0.0):
            return height*numpy.exp(-(((x - center)**2)*4.0*numpy.log(2.0))/
                            (fwhm**2)) + offset 
                                    
        circdichs2 = reset_tester()
        numpy.testing.assert_equal(circdichs2.axis.data, self.circdichs.axis.data)
        numpy.testing.assert_equal(circdichs2.data, self.circdichs.data)
        
        
        circdichs2.set_by_interpolation(self.circdichs.axis.data,\
                                                    self.circdichs.data)
        numpy.testing.assert_allclose(circdichs2.axis.data,\
                                        self.circdichs.axis.data, 1e-03)
        
        circdichs2.normalize2(1)
        a = numpy.copy(circdichs2.data)
        circdichs2.normalize()
        b = circdichs2.data
        numpy.testing.assert_equal(a, b)
        
        circdichs2.subtract(1)
        numpy.testing.assert_equal(numpy.max(circdichs2.data), 0)

        circdichs2 = reset_tester()
        a = numpy.copy(circdichs2.data)
        circdichs2.add_to_data(self.circdichs)
        numpy.testing.assert_equal(circdichs2.data, 2*self.circdichs.data)
        
#        circdichs2 = reset_tester()
#        axis = circdichs2.axis.data
#        circdichs2.data = gaussian(axis, 1.0, 11000.0, 100.0)
#        popt,_ = circdichs2.gaussian_fit(guess = [1.0001, 11000.0001, 100.0001])
#        gauss = gaussian(axis, popt[0], popt[1], popt[2])
#        numpy.testing.assert_allclose(gauss, self.circdichs.data)
      
        circdichs2.clear_data()
        numpy.testing.assert_equal(circdichs2.data, 0)
    
    def test_container_get_set(self):
        """Testing the getters and setters of the CircDichSpectrumContainer class
        
        """ 
        circdich1 = self.circdich1
        
        cntnr = CircDichSpectrumContainer()
        cntnr.set_axis(circdich1.axis)
        cntnr.set_spectrum(circdich1, tag='tester')
        spec = cntnr.get_spectrum(tag='tester')
        
        numpy.testing.assert_array_equal(circdich1.data, spec.data)
    
        
    def test_circdich_spectrum_saveablity(self):
        """Testing if CircDichSpectrumContainer is saveable
        """
        import tempfile
        circdich1 = self.circdich1
        cntnr = CircDichSpectrumContainer()
        cntnr.set_spectrum(circdich1, tag='tester')
        
        with tempfile.TemporaryFile() as f:
            cntnr.save(f)
            f.seek(0)
            circdich2 = CircDichSpectrumContainer()
            circdich2 = circdich2.load(f)

        data2 = circdich2.get_spectrum(tag='tester').data
        
        numpy.testing.assert_array_equal(circdich1.data, data2)

