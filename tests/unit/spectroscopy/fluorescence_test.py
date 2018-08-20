# -*- coding: utf-8 -*-
import unittest

"""
*******************************************************************************


    Tests of the quantarhei.spectroscopy.fluorescence package


*******************************************************************************
"""
import numpy
#import h5py
import copy

from quantarhei.spectroscopy.fluorescence import FluorSpectrumBase, FluorSpectrumContainer
from quantarhei import FrequencyAxis
from quantarhei import energy_units

from quantarhei import FluorSpectrum, FluorSpectrumCalculator
from quantarhei import Molecule, Aggregate, CorrelationFunction, TimeAxis

class TestFluor(unittest.TestCase):
    """Tests for the fluorescence package
    
    
    """
    
    def setUp(self,verbose=False):
        
        fluors = FluorSpectrumBase()
        with energy_units("1/cm"):
            f = FrequencyAxis(10000.0,2000, 1.0)
            a = self._spectral_shape(f, [1.0, 11000.0, 100.0])
            fluors.axis = f
            fluors.data = a
        self.fluors = fluors
        self.axis = fluors.axis
        
        #make a single-molecule system
        time = TimeAxis(0.0,1000,1.0)
        self.ta = time
        with energy_units("1/cm"):
            mol1 = Molecule(elenergies=[0.0, 12000.0])
            params = dict(ftype="OverdampedBrownian", reorg=20, cortime=100,
                          T=300)
            mol1.set_dipole(0,1,[0.0, 1.0, 0.0])
            cf = CorrelationFunction(time, params)
            mol1.set_transition_environment((0,1),cf)
            
            fluor_calc = FluorSpectrumCalculator(time, system=mol1)
            fluor_calc.bootstrap(rwa=12000)
            
        fluor1 = fluor_calc.calculate()
        
        self.fluor1 = fluor1
        
    def _spectral_shape(self, f, par):
        a = par[0]
        fd = par[1]
        sig = par[2]
        with energy_units("1/cm"):
            return a*numpy.exp(-((f.data-fd)/sig)**2)
        
    def _opt_spectral_shape(self, par):
        f = self.axis
        a = FluorSpectrumBase(axis=f,data=self._spectral_shape(f, par))
        return a
    
    def test_FluorSpectrumBase(self):
        """Testing the functions of the FluorSpectrumBase class
        
        """
        def reset_tester():
            fluors2 = FluorSpectrumBase()
            fluors2.set_axis(copy.deepcopy(self.fluors.axis))
            fluors2.set_data(copy.deepcopy(self.fluors.data))
            return fluors2
        
        def gaussian(x, height, center, fwhm, offset=0.0):
            return height*numpy.exp(-(((x - center)**2)*4.0*numpy.log(2.0))/
                            (fwhm**2)) + offset 
                                    
        fluors2 = reset_tester()
        numpy.testing.assert_equal(fluors2.axis.data, self.fluors.axis.data)
        numpy.testing.assert_equal(fluors2.data, self.fluors.data)
        
        
        fluors2.set_by_interpolation(self.fluors.axis.data,\
                                                    self.fluors.data)
        numpy.testing.assert_allclose(fluors2.axis.data,\
                                        self.fluors.axis.data, 1e-03)
        
        fluors2.normalize2(1)
        a = numpy.copy(fluors2.data)
        fluors2.normalize()
        b = fluors2.data
        numpy.testing.assert_equal(a, b)
        
        fluors2.subtract(1)
        numpy.testing.assert_equal(numpy.max(fluors2.data), 0)

        fluors2 = reset_tester()
        a = numpy.copy(fluors2.data)
        fluors2.add_to_data(self.fluors)
        numpy.testing.assert_equal(fluors2.data, 2*self.fluors.data)
        
#        fluors2 = reset_tester()
#        axis = fluors2.axis.data
#        fluors2.data = gaussian(axis, 1.0, 11000.0, 100.0)
#        popt,_ = fluors2.gaussian_fit(guess = [1.0001, 11000.0001, 100.0001])
#        gauss = gaussian(axis, popt[0], popt[1], popt[2])
#        numpy.testing.assert_allclose(gauss, self.fluors.data)
      
        fluors2.clear_data()
        numpy.testing.assert_equal(fluors2.data, 0)
              
    
    def test_container_get_set(self):
        """Testing the getters and setters of the FluorSpectrumContainer class
        
        """ 
        fluor1 = self.fluor1
        
        cntnr = FluorSpectrumContainer()
        cntnr.set_axis(fluor1.axis)
        cntnr.set_spectrum(fluor1, tag='tester')
        spec = cntnr.get_spectrum(tag='tester')
        
        numpy.testing.assert_array_equal(fluor1.data, spec.data)
    
        
    def test_fluor_spectrum_saveablity(self):
        """Testing if FluorSpectrumContainer is saveable
        
        """
        import tempfile
        fluor1 = self.fluor1
        cntnr = FluorSpectrumContainer()
        cntnr.set_spectrum(fluor1, tag='tester')
        
        with tempfile.TemporaryFile() as f:
            cntnr.save(f)
            f.seek(0)
            
            fluor2 = FluorSpectrumContainer()
            fluor2 = fluor2.load(f)

        data2 = fluor2.get_spectrum(tag='tester').data
        
        numpy.testing.assert_array_equal(fluor1.data, data2)
    
    def test_fluor_calculator(self):
        """Testing some basic methods of the FluorSpectrumCalculator class
        
        """
        
        #Calculations for monomer are already performed. Now testing that
        #fluorescence spectrum can be calculated for an aggregate.
        try:
            mol1 = Molecule(elenergies=[0.0, 12000.0])
            mol2 = Molecule(elenergies=[0.0, 12000.0])

            params = dict(ftype="OverdampedBrownian", reorg=20, cortime=100,
                          T=300)
            mol1.set_dipole(0,1,[0.0, 1.0, 0.0])
            mol2.set_dipole(0,1,[0.0, 1.0, 0.0])

            cf = CorrelationFunction(self.ta, params)
            mol1.set_transition_environment((0,1),cf)
            mol2.set_transition_environment((0,1),cf)

            agg = Aggregate(name="tester_dim", molecules=[mol1, mol2])
            agg.build()
            HH = agg.get_Hamiltonian()
            
        
            fluor_calc = FluorSpectrumCalculator(self.ta,
                 system=agg,
                 dynamics="secular",
                 relaxation_tensor=None,
                 rate_matrix=None,
                 effective_hamiltonian=HH,
                 temperature=300)
            fluor_calc.bootstrap(rwa=12000)           
            fluor_calc = fluor_calc.calculate()
        except:
            raise Exception('Fluorescence not calculatable for aggregate')
        

        
