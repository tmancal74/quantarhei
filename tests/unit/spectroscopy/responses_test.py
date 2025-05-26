# -*- coding: utf-8 -*-
import unittest
from pathlib import Path
import numpy.testing as npt

import tempfile
import os

import numpy

import quantarhei as qr

from quantarhei.utils.vectors import X 
from quantarhei.symbolic.cumulant import GFInitiator

#import matplotlib.pyplot as plt


"""
*******************************************************************************


    Tests of the quantarhei.spectroscopy.responses package


*******************************************************************************
"""

TEST_DIR = Path(__file__).parent


def efce_R1g_vec(t2, t1s, t3s, aa):
    """Our R1g response: vectorial version
    
    
    """
    import numpy as np

    aa.set_t2(t2)
    
    return np.exp(-aa.gg_t1["ab"] + aa.gg_t1_t2["ab"] - aa.gg_t1_t2_t3["aa"] - np.conj(aa.gg_t2["bb"])
                  + np.conj(aa.gg_t2_t3["ab"]) - np.conj(aa.gg_t3["ab"]))

def efce_R2g_vec(t2, t1s, t3s, aa):
    """Our R2g response: vectorial version
    
    
    """
    import numpy as np

    aa.set_t2(t2)
    
    return np.exp(aa.gg_t2["ab"] - aa.gg_t2_t3["bb"] - np.conj(aa.gg_t1["ab"]) 
            - np.conj(aa.gg_t1_t2["aa"])  + np.conj(aa.gg_t1_t2_t3["ab"]) - np.conj(aa.gg_t3["ab"]))

def efce_R3g_vec(t2, t1s, t3s, aa):
    """Our R23g response: vectorial version
    
    
    """
    import numpy as np

    aa.set_t2(t2)
    
    return np.exp(-aa.gg_t3["bb"] - np.conj(aa.gg_t1["aa"]) - np.conj(aa.gg_t1_t2["ab"]) 
                  + np.conj(aa.gg_t1_t2_t3["ab"]) + np.conj(aa.gg_t2["ab"]) - np.conj(aa.gg_t2_t3["ab"]))

def efce_R4g_vec(t2, t1s, t3s, aa):
    """Our R4g response: vectorial version
    
    
    """
    import numpy as np

    aa.set_t2(t2)
    
    return np.exp(-aa.gg_t1["aa"] + aa.gg_t1_t2["ab"] - aa.gg_t1_t2_t3["ab"] - aa.gg_t2["ab"]
                  + aa.gg_t2_t3["ab"] - aa.gg_t3["bb"])



class TestResponses(unittest.TestCase):
    """Tests for the response package
    
    
    """
    
    def setUp(self,verbose=False):
        pass

        
    def test_Responses(self):
        """Testing basic functions of the TwoDSpectrumBase class
        
        """
        #
        # Lineshape function
        #
        Nt1 = 100
        dt1 = 11
        Nt3 = 100
        dt3 = 12
        
        Nt2 = 50
        dt2 = 10.0
        
        t1_axis = qr.TimeAxis(0.0, Nt1, dt1)
        t3_axis = qr.TimeAxis(0.0, Nt3, dt3)
        t2_axis = qr.TimeAxis(0.0, Nt2, dt2)
        
        temperature = 300.0
        # parameters of the correlation function
        params = {"ftype":    "OverdampedBrownian",
                "reorg":    80.0,
                "cortime":  150.0,
                "T":        temperature,
                "matsubara":20}
        
        deg = [1.0, 0.0, 0.0]
        dfe = [1.0, 0.0, 0.0]
        
        d1 = deg
        d2 = deg
        d3 = deg
        d4 = deg
        d3f = dfe
        d4f = dfe
        
        omeg = 12100.0
        omfe = 12100.0
        
        #
        # we get the g-functions as spline interpolated functions
        #
        tt = qr.TimeAxis(0.0, 3*t1_axis.length, t1_axis.step)
        with qr.energy_units("1/cm"):
            cf = qr.LineshapeFunction(tt, params)
        
        sfce = cf.as_spline_function()

        responses = []
        
        response_types = ["R1g", "R2g", "R3g", "R4g"]
        response_functions = [efce_R1g_vec, efce_R2g_vec, efce_R3g_vec, efce_R4g_vec]
        
        gdict = dict(aa=sfce, bb=sfce, ab=sfce, ba=sfce)
        
        gf = GFInitiator(t1_axis.data, t3_axis.data, gdict)
        
        for rftype, rffce in zip(response_types, response_functions):
            rfce = qr.ResponseFunction(rftype)
            rfce.set_evaluation_function(rffce)
            responses.append(rfce)
        
        #
        # This part is the same for all two-level system responses
        #
        for rfce in responses:
            rfce.set_dipoles(d1,d2,d3,d4)
            args = (gf,)
            rfce.set_auxliary_arguments(args)
            
            
        
        #
        # Set transition frequencies for all responses (beware of units)
        #
        with qr.energy_units("1/cm"):
            for rsp in responses:
                rsp.set_frequencies(omeg, omeg)

    

        #
        # Calculation of the 2D spectra from response functions
        #
        
        rwa = 12000.0
        Npad = 0
        
        # laboratory settings
        lab = qr.LabSetup()
        lab.set_pulse_polarizations(pulse_polarizations=(X,X,X), detection_polarization=X)
        
        
        calc = qr.TwoDResponseCalculator(t1_axis, t2_axis, t3_axis, responses=responses)
        with qr.energy_units("1/cm"):
            calc.bootstrap(rwa=rwa, pad=Npad, lab=lab, verbose=True) 
            
        one = False
        T2 = 0.0
        
        if one:
            print("Calculating one spectrum")
            twod = calc.calculate_one(0) 
        
        else:
            print("Calculating", Nt2,"spectra")
            tcont = calc.calculate()
            tcont = tcont.get_TwoDSpectrumContainer()
            twod = tcont.get_spectrum(T2)            


        file_path_1 = TEST_DIR / "responses_test_twod_0.dat"
        save_data = False
        if save_data:
            twod.save_data(file_path_1)

        #
        # Here we load spectrum for comparison
        #
        twod0 = qr.TwoDSpectrum()
        twod0.load_data(file_path_1)
        twod0.set_axis_1(t1_axis)
        twod0.set_axis_3(t3_axis)


        npt.assert_allclose(twod0.data, twod.data)


    def test_LouvillePathway(self):
        """Testing basic functions of the TwoDSpectrumCalculator class
        
        """
        t1 = qr.TimeAxis(0.0, 1000, 1.0)
        t3 = qr.TimeAxis(0.0, 1000, 1.0)
        
        t2 = qr.TimeAxis(30, 10, 10.0)
        
        #twod_calc = qr.TwoDResponseCalculator(t1, t2, t3)
        
        
            
    
