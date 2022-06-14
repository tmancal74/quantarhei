# -*- coding: utf-8 -*-
import unittest
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.qm.TDRedfieldRelaxationTensor against 
                 quantarhei.qm.RefieldRateMatrix class


*******************************************************************************
"""


from quantarhei.qm import RedfieldRateMatrix
from quantarhei.qm import TDRedfieldRelaxationTensor
from quantarhei import Hamiltonian
from quantarhei import CorrelationFunction
from quantarhei.qm.corfunctions import CorrelationFunctionMatrix 
from quantarhei.qm.corfunctions import SpectralDensity
from quantarhei import TimeAxis, eigenbasis_of
from quantarhei.qm import Operator
from quantarhei.qm import SystemBathInteraction

from quantarhei import energy_units

import quantarhei as qr

class TDTestRedfield(unittest.TestCase):
    """Tests for the TDRedfieldRelaxationTensor class
    
    
    """
    
    def setUp(self,verbose=False):
        
        self.verbose = verbose
        
        time = TimeAxis(0.0,1000,1.0)
        with energy_units("1/cm"):
            params = {"ftype":"OverdampedBrownian",
                      "reorg":30.0,
                      "T":300.0,
                      "cortime":100.0}
                      
            cf1 = CorrelationFunction(time,params)
            cf2 = CorrelationFunction(time,params)
            sd = SpectralDensity(time,params)
                    
        cm1 = CorrelationFunctionMatrix(time,2,1)
        cm1.set_correlation_function(cf1,[(0,0),(1,1)])
        cm2 = CorrelationFunctionMatrix(time,2,1)
        cm2.set_correlation_function(cf2,[(0,0),(1,1)])  
          
        K11 = numpy.array([[1.0, 0.0],[0.0, 0.0]],dtype=qr.REAL)
        K21 = numpy.array([[1.0, 0.0],[0.0, 0.0]],dtype=qr.REAL)
        K12 = K11.copy()
        K22 = K21.copy()
            
        KK11 = Operator(data=K11)
        KK21 = Operator(data=K21)
        KK12 = Operator(data=K12)
        KK22 = Operator(data=K22) 
        
        self.sbi1 = SystemBathInteraction([KK11,KK21],cm1)
        self.sbi2 = SystemBathInteraction([KK12,KK22],cm2)
        
        with energy_units("1/cm"):
            h1 = [[0.0, 100.0],[100.0, 0.0]]
            h2 = [[0.0, 100.0],[100.0, 0.0]]
            self.H1 = Hamiltonian(data=h1)
            self.H2 = Hamiltonian(data=h2)
            
        #sd.convert_2_spectral_density()
        with eigenbasis_of(self.H1):
            de = self.H1.data[1,1]-self.H1.data[0,0]
        self.c_omega_p = sd.at(de,approx="spline") #interp_data(de)
        self.c_omega_m = sd.at(-de,approx="spline") #interp_data(-de)
        

    
    def test_comparison_of_rates(self):
        """Testing that Time-dependent Redfield tensor and rate matrix are compatible
        
        
        """ 
        tensor = True
        matrix = True
        
        dim = self.H1.dim
        KT = numpy.zeros((dim,dim), dtype=numpy.float64)
        KM = numpy.zeros((dim,dim), dtype=numpy.float64)
        
        if tensor:
            #print(self.H1)
            RT = TDRedfieldRelaxationTensor(self.H1, self.sbi1)
            
            timeaxis = self.sbi1.CC.timeAxis
            Nt = timeaxis.length - 1
            
            for n in range(2):
                for m in range(2):
                    #print(n,m,numpy.real(RT.data[n,n,m,m]))
                    KT[n,m] = numpy.real(RT.data[Nt,n,n,m,m])
                    
                    
        if matrix:
            #print(self.H2)
            RR = RedfieldRateMatrix(self.H2,self.sbi2)
            
            for n in range(2):
                for m in range(2):
                    #print(n,m,numpy.real(RR.data[n,m]))    
                    KM[n,m] = numpy.real(RR.data[n,m])
                    
        #print(self.c_omega_p)
        #print(self.c_omega_m)
        
        numpy.testing.assert_allclose(KT,KM, rtol=1.0e-2)
                
            

