# -*- coding: utf-8 -*-
import unittest

"""
*******************************************************************************


    Tests of the quantarhei.spectroscopy.abs package


*******************************************************************************
"""
import numpy
#import h5py
import tempfile

from quantarhei.spectroscopy.abs2 import AbsSpectrumBase #, AbsSpectrumDifference
from quantarhei import FrequencyAxis
from quantarhei import energy_units, eigenbasis_of
#from quantarhei import convert

from quantarhei import AbsSpectrum, AbsSpectrumCalculator
from quantarhei import Molecule, CorrelationFunction, TimeAxis
from quantarhei import Aggregate, Mode

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
        
        
        
        
        time = TimeAxis(0.0,1000,1.0)
        with energy_units("1/cm"):
            mol1 = Molecule(elenergies=[0.0, 10000.0])
            params = dict(ftype="OverdampedBrownian", reorg=30, cortime=100,
                          T=300)
            mol1.set_dipole(0,1,[0.0, 1.0, 0.0])
            cf = CorrelationFunction(time, params)
            mol1.set_transition_environment((0,1),cf)
            
            self.mol1 = mol1
            
            abs_calc = AbsSpectrumCalculator(time, system=mol1)
            abs_calc.bootstrap(rwa=10000)

            mol3 = Molecule(elenergies=[0.0, 10000.0])
            params = dict(ftype="OverdampedBrownian", reorg=30, cortime=100,
                          T=300)
            mol3.set_dipole(0,1,[0.0, 1.0, 0.0])
            cf = CorrelationFunction(time, params)
            mol3.set_transition_environment((0,1),cf)

            self.mol3 = mol3
            
            if True:
                mod = Mode(frequency=1000.0)
                mol3.add_Mode(mod)
                mod.set_nmax(0, 1)
                mod.set_nmax(1, 4)
                mod.set_HR(1, 0.1)
            
        abs1 = abs_calc.calculate()
        
        self.abs1 = abs1
        
        with energy_units("1/cm"):
            mol2 = Molecule(elenergies=[0.0, 10000.0])
            mol2.set_dipole(0,1,[1.0, 0.0, 0.0])
            mol2.set_transition_environment((0,1), cf)
            
        agg = Aggregate(molecules=[mol1, mol2])
        agc = Aggregate(molecules=[mol1, mol2])
        
        with energy_units("1/cm"):
            agg.set_resonance_coupling(0,1,0.0)
            agc.set_resonance_coupling(0,1,300.0)
        
        agg.build()
        agc.build()
        self.agg0 = agg
        self.agg = agc
        
        
    def _spectral_shape(self, f, par):
        a = par[0]
        fd = par[1]
        sig = par[2]
        with energy_units("1/cm"):
            return a*numpy.exp(-((f.data-fd)/sig)**2)
        
    def _opt_spectral_shape(self, par):
        f = self.axis
        a = AbsSpectrumBase(axis=f,data=self._spectral_shape(f, par))
        return a
        
            
    def test_abs_spectrum_saveablity(self):
        """Testing if AbsSpectrum is saveble
        
        """
        abs1 = self.abs1
        
        #drv = "core"
        #bcs = False
        
        with tempfile.TemporaryFile() as f:
        #with h5py.File('tempfile.hdf5', 
        #               driver=drv, 
        #               backing_store=bcs) as f:
    
            abs1.save(f) #, test=True)
            f.seek(0)
            abs2 = AbsSpectrum()
            abs2 = abs2.load(f) #, test=True)
            
        numpy.testing.assert_array_equal(abs1.data, abs2.data)


    # def test_abs_difference_equal_zero(self):
    #     """Testing calculation of zero difference between abs spectra
        
    #     """
        
    #     #
    #     #  Test difference = 0
    #     # 
    #     par = [1.0, 11000.0, 100.0]
    #     with energy_units("1/cm"):
    #         f2 = self._spectral_shape(self.abss.axis, par) 
        
    #         abss2 = AbsSpectrumBase(axis=self.abss.axis, data=f2)
    #         ad = AbsSpectrumDifference(target=self.abss, optfce=abss2, 
    #                                bounds=(10100.0, 11900.0))
    #     d = ad.difference()
        
    #     self.assertAlmostEqual(d, 0.0)
        
        
        
    # def test_abs_difference_formula(self):
    #     """Testing formula for calculation of difference between abs spectra
        
    #     """

    #     #
    #     #  Test difference formula
    #     #
    #     par = [1.0, 10900.0, 100.0]
    #     with energy_units("1/cm"):
    #         f2 = self._spectral_shape(self.abss.axis, par) 
        
    #         abss2 = AbsSpectrumBase(axis=self.abss.axis, data=f2)
    #         ad = AbsSpectrumDifference(target=self.abss, optfce=abss2, 
    #                                bounds=(10100.0, 11900.0))
    #     d = ad.difference()  
          
    #     #with energy_units("1/cm"):   
    #     if True:
    #         target = self.abss.data[ad.nl:ad.nu]
    #         secabs = abss2.data[ad.nl:ad.nu]
    #         x = self.abss.axis.data[ad.nl:ad.nu]
    #         d2 = 1000.0*numpy.sum(numpy.abs((target-secabs)**2/
    #                                (x[len(x)-1]-x[0])))
        
    #     self.assertAlmostEqual(d, d2)        
        
#    def test_minimize(self):
#        """Testing difference minimization using scipy.optimize package
#        
#        """
#        with energy_units("1/cm"):
#            ad = AbsSpectrumDifference(target=self.abss, 
#                                   optfce=self._opt_spectral_shape, 
#                                   bounds=(10100.0, 11900.0))
#        
#        method = "Nelder-Mead"
#        ini = [0.5, 10900, 80.0]
#        p = ad.minimize(ini, method=method)
#        
#        numpy.testing.assert_array_almost_equal(p, [1.0, 11000.0, 100.0])


    def test_of_molecular_absorption(self):
        """(AbsSpectrum) Testing absorption spectrum of a molecule
        
        
        """
        mol1 = self.mol1
        mol1.set_electronic_rwa([0, 1])
        time = mol1.get_transition_environment((0,1)).axis
                      
        with energy_units("1/cm"):
            # data for comparison
            x = self.abs1.axis.data
            y = self.abs1.data/3.0        
                                               
        abs_calc = AbsSpectrumCalculator(time, system=mol1)
        

        prop = mol1.get_ReducedDensityMatrixPropagator(time, 
                                                relaxation_theory="stR",
                                                time_dependent=True)
                
        abs_calc.bootstrap(prop=prop)
        abs1 = abs_calc.calculate(from_dynamics=True)   
        
        with energy_units("1/cm"):
            x1 = abs1.axis.data
            y1 = abs1.data 
        
        diff = numpy.max(numpy.abs(y1-y))
        rdiff = diff/numpy.max(numpy.abs(y))
        
        self.assertTrue(rdiff < 0.01)
        
        _plot_ = False
        if _plot_:
            import matplotlib.pyplot as plt
            plt.plot(x,y,"-b")
            plt.plot(x1,y1,"--r")
            plt.show() 

        mol2 = self.mol3
        mol2.set_electronic_rwa([0, 1])        
        abs_calc = AbsSpectrumCalculator(time, system=mol2)
         
        
        prop = mol2.get_ReducedDensityMatrixPropagator(time, 
                                                 relaxation_theory="stR",
                                                 time_dependent=True)
                 
        abs_calc.bootstrap(prop=prop)
        ham = mol2.get_Hamiltonian()
        #with energy_units("1/cm"):
        #    with eigenbasis_of(ham):
        #        print("ham.dat")
        #        print(numpy.diag(ham.data))
                
        #print("****************")

        #
        # we test that calculate cannot be called within a basis context   
        #
        with self.assertRaises(Exception) as exp:
            with eigenbasis_of(ham):
                abs2 = abs_calc.calculate(from_dynamics=True)  
             
        abs2 = abs_calc.calculate(from_dynamics=True)
        abs3 = abs_calc.calculate(from_dynamics=True, alt=True)
        #print("****************")


        with energy_units("1/cm"):
            x1 = abs2.axis.data
            y1 = abs2.data
            y = abs3.data

        diff = numpy.max(numpy.abs(y1-y))
        rdiff = diff/numpy.max(numpy.abs(y))
        self.assertTrue(rdiff < 0.01)
        
        _plot_ = False
        if _plot_:
            import matplotlib.pyplot as plt
            #plt.plot(x,y,"-b")
            plt.plot(x1,y1,"--r")
            plt.show() 
            
            #raise Exception()
            
    
    def test_of_aggregate_absorption(self):
        """(AbsSpectrum) Testing absorption spectrum of an aggregate
        
        
        """
        mol1 = self.mol1
        agg = self.agg0
        agc = self.agg
        
        time = mol1.get_transition_environment((0,1)).axis
                      
        with energy_units("1/cm"):
            # data for comparison
            x = self.abs1.axis.data
            y = self.abs1.data/3.0
                                                
        abs_calc = AbsSpectrumCalculator(time, system=agg)
        
        prop = agg.get_ReducedDensityMatrixPropagator(time, 
                                                relaxation_theory="stR",
                                                time_dependent=True)

        abs_calc.bootstrap(prop=prop)
        
        abs1 = abs_calc.calculate(from_dynamics=True)   
        
        with energy_units("1/cm"):
            x1 = abs1.axis.data
            y1 = abs1.data/2.0   # we divide by 2, becuase there are 2
                                 # molecules in the aggregate
        
        diff = numpy.max(numpy.abs(y1-y))
        rdiff = diff/numpy.max(numpy.abs(y))
        
        self.assertTrue(rdiff < 0.01)
        
        _plot_ = False
        if _plot_:
            import matplotlib.pyplot as plt
            plt.plot(x,y,"-b")
            plt.plot(x1,y1,"--r")
            plt.show() 
        
        abs_calc = AbsSpectrumCalculator(time, system=agc)
        
        prop = agc.get_ReducedDensityMatrixPropagator(time, 
                                                relaxation_theory="stR",
                                                time_dependent=True)
        
        abs_calc.bootstrap(prop=prop)
        abs1 = abs_calc.calculate(from_dynamics=True)   
        
        with energy_units("1/cm"):
            x1 = abs1.axis.data
            y1 = abs1.data   

        _plot_ = False
        if _plot_:
            import matplotlib.pyplot as plt
            plt.plot(x,y,"-b")
            plt.plot(x1,y1,"--r")
            plt.show()     
            


if __name__ == '__main__':
    unittest.main()        
        