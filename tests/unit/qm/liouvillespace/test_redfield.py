# -*- coding: utf-8 -*-

import unittest
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.qm.RedfieldRelaxationTensor and 
                 quantarhei.qm.RefieldRateMatrix classes


*******************************************************************************
"""


from quantarhei.qm import RedfieldRateMatrix
from quantarhei.qm import RedfieldRelaxationTensor
from quantarhei import Hamiltonian
from quantarhei import CorrelationFunction
from quantarhei.qm.corfunctions import CorrelationFunctionMatrix 
from quantarhei.qm.corfunctions import SpectralDensity
from quantarhei import TimeAxis, eigenbasis_of
from quantarhei.qm import Operator
from quantarhei.qm import SystemBathInteraction
from quantarhei.qm import ReducedDensityMatrixPropagator
from quantarhei.qm import ReducedDensityMatrix
from quantarhei import energy_units

class TestRedfield(unittest.TestCase):
    """Tests for the RedfieldRelaxationTensor class
    
    
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
          
        K11 = numpy.array([[1.0, 0.0],[0.0, 0.0]],dtype=numpy.float)
        K21 = numpy.array([[1.0, 0.0],[0.0, 0.0]],dtype=numpy.float)
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
        """Testing that Redfield tensor and rate matrix are compatible
        
        
        """ 
        tensor = True
        matrix = True
        
        dim = self.H1.dim
        KT = numpy.zeros((dim,dim), dtype=numpy.float64)
        KM = numpy.zeros((dim,dim), dtype=numpy.float64)
        
        if tensor:
            #print(self.H1)
            RT = RedfieldRelaxationTensor(self.H1,self.sbi1)
            
            for n in range(2):
                for m in range(2):
                    #print(n,m,numpy.real(RT.data[n,n,m,m]))
                    KT[n,m] = numpy.real(RT.data[n,n,m,m])
                    
                    
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
                

    def test_propagation_in_different_basis(self):
        """(REDFIELD) Testing comparison of propagations in different bases

        """

        LT1 = RedfieldRelaxationTensor(self.H1, self.sbi1, as_operators=True)
        LT2 = RedfieldRelaxationTensor(self.H1, self.sbi1, as_operators=False)
        
        time = TimeAxis(0.0, 1000, 1.0)
        
        prop1 = ReducedDensityMatrixPropagator(time, self.H1, LT1)
        prop2 = ReducedDensityMatrixPropagator(time, self.H1, LT2)
        
        rho0 = ReducedDensityMatrix(dim=self.H1.dim)
        rho0.data[1,1] = 1.0
          
        with eigenbasis_of(self.H1):
        #if True:
            rhot1_e = prop1.propagate(rho0)
            
        with eigenbasis_of(self.H1):
            rhot2_e = prop2.propagate(rho0)

        rhot1_l = prop1.propagate(rho0)
        rhot2_l = prop2.propagate(rho0)
            
        numpy.testing.assert_allclose(rhot1_l.data, rhot1_e.data,
                                      rtol=1.0e-5, atol=1.0e-12)
        numpy.testing.assert_allclose(rhot2_l.data, rhot1_e.data,
                                      rtol=1.0e-5, atol=1.0e-12)
        numpy.testing.assert_allclose(rhot1_e.data, rhot2_e.data,
                                      rtol=1.0e-5, atol=1.0e-12)              
