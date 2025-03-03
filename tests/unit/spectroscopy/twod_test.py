# -*- coding: utf-8 -*-
import unittest

from pathlib import Path

import tempfile
import os

import numpy
import numpy.testing as npt

from nose.tools import assert_raises

import quantarhei as qr
from quantarhei.spectroscopy.twod2 import TwoDSpectrumBase
from quantarhei.utils.vectors import X 

import matplotlib.pyplot as plt


"""
*******************************************************************************


    Tests of the quantarhei.spectroscopy.twod package


*******************************************************************************
"""
TEST_DIR = Path(__file__).parent


class TestTwod(unittest.TestCase):
    """Tests for the twod package
    
    
    """
    
    def setUp(self,verbose=False):
 
        
        #######################################################################
        #
        # Create a system to calculate 2D spectrum with Mock calculator
        #
        #######################################################################
        
        #
        # System parameters
        #
        E0 = 12000.0
        E1 = E0 + 100.0
        
        
        self.E0 = E0
        self.E1 = E1
        

        with qr.energy_units("1/cm"):
            m1 = qr.Molecule([0.0, E0])
            m2 = qr.Molecule([0.0, E1])
            
            m1.set_dipole(0,1,[2.0, 0.0, 0.0])
            m2.set_dipole(0,1,[0.5, 0.5, 0.0])
            
            m1.position = [0.0, 0.0, 0.0]
            m2.position = [0.0, 0.0, 2.0]
            
            
        agg = qr.Aggregate(molecules=[m1, m2])
        agg.set_coupling_by_dipole_dipole()
        
        with qr.energy_units("1/cm"):
            m1.set_transition_width((0,1), 200)
            m2.set_transition_width((0,1), 200)
            
        self.agg2D = agg.deepcopy()
        
        agg.build(mult=1)
        
        self.H = agg.get_Hamiltonian()
        
        gammas = [1.0/200.0]
        with qr.eigenbasis_of(self.H):
            
            K = qr.qm.ProjectionOperator(1, 2, dim=self.H.dim)
            
        sbi = qr.qm.SystemBathInteraction(sys_operators=[K], rates=gammas,
                                          system=agg)
        
        self.LL = qr.qm.LindbladForm(self.H, sbi)




    def test_MockTwoD(self):
        """Testing calculation with artificial lineshapes
        
        """
        
        agg2D = self.agg2D
        H = self.H
        LL = self.LL            
        
        t2axis = qr.TimeAxis(0.0, 200, 5.0)
        t1axis = qr.TimeAxis(0.0, 100, 10.0)
        t3axis = qr.TimeAxis(0.0, 100, 10.0)
        
        eUt = qr.EvolutionSuperOperator(time=t2axis, ham=H, relt=LL)
        eUt.set_dense_dt(10)
        eUt.calculate()

        msc = qr.MockTwoDResponseCalculator(t1axis, t2axis, t3axis)
        msc.bootstrap(rwa=qr.convert((self.E0+self.E1)/2.0,"1/cm","int"), 
                      shape="Gaussian")
    
        agg2D.build(mult=2)
        agg2D.diagonalize()

        lab = qr.LabSetup()
        lab.set_pulse_polarizations(pulse_polarizations=[X,X,X],
                              detection_polarization=X)   
        
        #
        # Calculation one by one
        # 
        cont = qr.TwoDSpectrumContainer(t2axis=t2axis)
        
        for t2 in t2axis.data:
            twod = msc.calculate_one_system(t2, agg2D, eUt, lab) 
        
            cont.set_spectrum(twod)        
        
        cont.set_data_flag(qr.signal_TOTL)
        
        sp = cont.get_spectrum(800.0)
        
        file_path_1 = TEST_DIR / "twod_test_sp_800.dat"
        
        #
        # This is how you save the data for future comparison
        #
        save_data = False
        if save_data:
            sp.save_data(file_path_1)
                         
        
        #
        # Here we load spectrum for comparison
        #
        sp800 = qr.TwoDSpectrum()
        sp800.load_data(file_path_1)
        
        #
        # At the moment, axes are not saved with the spectrum and we have to
        # set them manually.
        #
        sp800.set_axis_1(t1axis)
        sp800.set_axis_3(t3axis)
        
        """
        with qr.energy_units("1/cm"):
            spl.plot(show=True)
            sp.savefig("twod.png")
        """
            
        #
        # Calculation of all at once
        #
        cont2 = msc.calculate_all_system(agg2D, eUt, lab)
        cont2.set_data_flag(qr.signal_TOTL)
        sp2 = cont.get_spectrum(800.0)
        
        
        #
        # comparison of 2D spectra at 800 fs
        #
        npt.assert_allclose(sp800.data, sp2.data)
        npt.assert_allclose(sp800.data, sp.data)


        """ 
        with qr.energy_units("1/cm"):
            sp2.plot(show=True)
            sp2.savefig("twod.png")
            
        
            
       
#            
#            cont.make_movie("movie.mp4")
          
        pevol = cont.get_point_evolution(qr.convert(12600,"1/cm","int"),
                                         qr.convert(11500,"1/cm","int"),
                                         t2axis)
        
        plt.plot(t2axis.data, pevol)
        plt.show()
        

        pars = cont.global_fit_exponential()
        
        print("Optimum: ", pars.x)
        
        raise Exception()
        """
        
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
        
        twodB._add_data(data1, resolution="off", dtype=qr.signal_TOTL)

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


    
    def test_TwoDResponseCalculator(self):
        """Testing basic functions of the TwoDSpectrumCalculator class
        
        """
        t1 = qr.TimeAxis(0.0, 1000, 1.0)
        t3 = qr.TimeAxis(0.0, 1000, 1.0)
        
        t2 = qr.TimeAxis(30, 10, 10.0)
        
        twod_calc = qr.TwoDResponseCalculator(t1, t2, t3)
        
        
            
    