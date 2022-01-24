# -*- coding: utf-8 -*-
import unittest
import copy
"""
*******************************************************************************


    Tests of the quantarhei.spectroscopy.linear_dichroism package


*******************************************************************************
"""
import numpy
#import h5py

from quantarhei.spectroscopy.linear_dichroism import LinDichSpectrumBase, LinDichSpectrumContainer
from quantarhei import FrequencyAxis
from quantarhei import energy_units

from quantarhei import LinDichSpectrum, LinDichSpectrumCalculator
from quantarhei import Molecule, Aggregate, CorrelationFunction, TimeAxis

class TestLinDich(unittest.TestCase):
    """Tests for the linear dichroism package
    
    
    """
    
    def setUp(self,verbose=False):
        
        lindichs = LinDichSpectrumBase()
        with energy_units("1/cm"):
            f = FrequencyAxis(10000.0,2000, 1.0)
            a = self._spectral_shape(f, [1.0, 11000.0, 100.0])
            lindichs.axis = f
            lindichs.data = a
        self.lindichs = lindichs
        self.axis = lindichs.axis
        
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
            
            lindich_calc = LinDichSpectrumCalculator(time, system=agg, \
                                vector_perp_to_membrane=numpy.array((0,1,0)))
            lindich_calc.bootstrap(rwa=12000)
            
        lindich1 = lindich_calc.calculate()
        
        self.lindich1 = lindich1
        
    def _spectral_shape(self, f, par):
        a = par[0]
        fd = par[1]
        sig = par[2]
        with energy_units("1/cm"):
            return a*numpy.exp(-((f.data-fd)/sig)**2)
        
    def _opt_spectral_shape(self, par):
        f = self.axis
        a = LinDichSpectrumBase(axis=f,data=self._spectral_shape(f, par))
        return a
    
    def test_LinDichSpectrumBase(self):
        """Testing the functions of the LinDichSpectrumBase class
        
        """
        def reset_tester():
            lindichs2 = LinDichSpectrumBase()
            lindichs2.set_axis(copy.deepcopy(self.lindichs.axis))
            lindichs2.set_data(copy.deepcopy(self.lindichs.data))
            return lindichs2
        
        def gaussian(x, height, center, fwhm, offset=0.0):
            return height*numpy.exp(-(((x - center)**2)*4.0*numpy.log(2.0))/
                            (fwhm**2)) + offset 
                                    
        lindichs2 = reset_tester()
        numpy.testing.assert_equal(lindichs2.axis.data, self.lindichs.axis.data)
        numpy.testing.assert_equal(lindichs2.data, self.lindichs.data)
        
        
        lindichs2.set_by_interpolation(self.lindichs.axis.data,\
                                                    self.lindichs.data)
        numpy.testing.assert_allclose(lindichs2.axis.data,\
                                        self.lindichs.axis.data, 1e-03)
        
        lindichs2.normalize2(1)
        a = numpy.copy(lindichs2.data)
        lindichs2.normalize()
        b = lindichs2.data
        numpy.testing.assert_equal(a, b)
        
        lindichs2.subtract(1)
        numpy.testing.assert_equal(numpy.max(lindichs2.data), 0)

        lindichs2 = reset_tester()
        a = numpy.copy(lindichs2.data)
        lindichs2.add_to_data(self.lindichs)
        numpy.testing.assert_equal(lindichs2.data, 2*self.lindichs.data)
        
#        lindichs2 = reset_tester()
#        axis = lindichs2.axis.data
#        lindichs2.data = gaussian(axis, 1.0, 11000.0, 100.0)
#        popt,_ = lindichs2.gaussian_fit(guess = [1.0001, 11000.0001, 100.0001])
#        gauss = gaussian(axis, popt[0], popt[1], popt[2])
#        numpy.testing.assert_allclose(gauss, self.lindichs.data)
      
        lindichs2.clear_data()
        numpy.testing.assert_equal(lindichs2.data, 0)
    
    def test_container_get_set(self):
        """Testing the getters and setters of the LinDichSpectrumContainer class
        
        """ 
        lindich1 = self.lindich1
        
        cntnr = LinDichSpectrumContainer()
        cntnr.set_axis(lindich1.axis)
        cntnr.set_spectrum(lindich1, tag='tester')
        spec = cntnr.get_spectrum(tag='tester')
        
        numpy.testing.assert_array_equal(lindich1.data, spec.data)
    
        
    def test_lindich_spectrum_saveablity(self):
        """Testing if LinDichSpectrumContainer is saveable
        """
        
        import tempfile
        lindich1 = self.lindich1
        cntnr = LinDichSpectrumContainer()
        cntnr.set_spectrum(lindich1, tag='tester')
        
        with tempfile.TemporaryFile() as f:
            cntnr.save(f)
            f.seek(0)
            lindich2 = LinDichSpectrumContainer()
            lindich2 = lindich2.load(f)

        data2 = lindich2.get_spectrum(tag='tester').data
        
        numpy.testing.assert_array_equal(lindich1.data, data2)


