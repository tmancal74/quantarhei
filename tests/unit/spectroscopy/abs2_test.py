# -*- coding: utf-8 -*-
import unittest

"""
*******************************************************************************


    Tests of the quantarhei.spectroscopy.abs package


*******************************************************************************
"""
import numpy
#import h5py
import copy

#from quantarhei.spectroscopy.absbase import AbsSpectrumBase
from quantarhei import FrequencyAxis
from quantarhei import energy_units

from quantarhei import AbsSpectrum, AbsSpectrumCalculator, AbsSpectrumContainer
from quantarhei import Molecule, Aggregate, CorrelationFunction, TimeAxis

class TestAbs(unittest.TestCase):
    """Tests for the abs package
    
    
    """
    
    def setUp(self,verbose=False):
        
        abss = AbsSpectrum()
        with energy_units("1/cm"):
            f = FrequencyAxis(10000.0,2000, 1.0)
            a = self._spectral_shape(f, [1.0, 11000.0, 100.0])
            abss.axis = f
            abss.data = a
        self.abss = abss
        self.axis = abss.axis
        
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
            
            abs_calc = AbsSpectrumCalculator(time, system=mol1)
            abs_calc.bootstrap(rwa=12000)
            
        abs1 = abs_calc.calculate()
        
        self.abs1 = abs1
        
    def _spectral_shape(self, f, par):
        a = par[0]
        fd = par[1]
        sig = par[2]
        with energy_units("1/cm"):
            return a*numpy.exp(-((f.data-fd)/sig)**2)
        
    def _opt_spectral_shape(self, par):
        f = self.axis
        a = AbsSpectrum(axis=f,data=self._spectral_shape(f, par))
        return a
    
    
    def test_AbsSpectrumBase(self):
        """Testing the functions of the AbsSpectrumBase class
        
        """
        def reset_tester():
            abss2 = AbsSpectrum()
            abss2.set_axis(copy.deepcopy(self.abss.axis))
            abss2.set_data(copy.deepcopy(self.abss.data))
            return abss2
        
        def gaussian(x, height, center, fwhm, offset=0.0):
            return height*numpy.exp(-(((x - center)**2)*4.0*numpy.log(2.0))/
                            (fwhm**2)) + offset 
                                    
        abss2 = reset_tester()
        numpy.testing.assert_equal(abss2.axis.data, self.abss.axis.data)
        numpy.testing.assert_equal(abss2.data, self.abss.data)
        
        
        abss2.set_by_interpolation(self.abss.axis.data,\
                                                    self.abss.data)
        numpy.testing.assert_allclose(abss2.axis.data,\
                                        self.abss.axis.data, 1e-03)
        
        abss2.normalize2(1)
        a = numpy.copy(abss2.data)
        abss2.normalize()
        b = abss2.data
        numpy.testing.assert_equal(a, b)
        
        abss2.subtract(1)
        numpy.testing.assert_equal(numpy.max(abss2.data), 0)

        abss2 = reset_tester()
        a = numpy.copy(abss2.data)
        abss2.add_to_data(self.abss)
        numpy.testing.assert_equal(abss2.data, 2*self.abss.data)
        
#        abss2 = reset_tester()
#        axis = abss2.axis.data
#        abss2.data = gaussian(axis, 1.0, 11000.0, 100.0)
#        popt,_ = abss2.gaussian_fit(guess = [1.0001, 11000.0001, 100.0001])
#        gauss = gaussian(axis, popt[0], popt[1], popt[2])
#        numpy.testing.assert_allclose(gauss, self.abss.data)
      
        abss2.clear_data()
        numpy.testing.assert_equal(abss2.data, 0)
              
    
    def test_container_get_set(self):
        """Testing the getters and setters of the AbsSpectrumContainer class
        
        """ 
        abs1 = self.abs1
        
        cntnr = AbsSpectrumContainer()
        cntnr.set_axis(abs1.axis)
        cntnr.set_spectrum(abs1, tag='tester')
        spec = cntnr.get_spectrum(tag='tester')
        
        numpy.testing.assert_array_equal(abs1.data, spec.data)
    
        
    def test_abs_spectrum_saveablity(self):
        """Testing if AbsSpectrumContainer is saveable
        
        """
        import tempfile
        abs1 = self.abs1
        cntnr = AbsSpectrumContainer()
        cntnr.set_spectrum(abs1, tag='tester')
        
        with tempfile.TemporaryFile() as f:
            cntnr.save(f) #'tempfile', test=True)
            f.seek(0)
            abs2 = AbsSpectrumContainer()
            abs2 = abs2.load(f) #'tempfile')

        data2 = abs2.get_spectrum(tag='tester').data
        
        numpy.testing.assert_array_equal(abs1.data, data2)
        
        
    def test_abs_calculator(self):
        """Testing some basic methods of the AbsSpectrumCalculator class
        
        """
        
        #Calculations for monomer are already performed. Now testing that
        #absorption spectrum can be calculated for an aggregate.
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
            
        
            abs_calc = AbsSpectrumCalculator(self.ta,
                 system=agg,
                 dynamics="secular",
                 relaxation_tensor=None,
                 rate_matrix=None)
            abs_calc.bootstrap(rwa=12000)           
            abs_calc = abs_calc.calculate()
        except:
            raise Exception('Absorption not calculatable for aggregate')

