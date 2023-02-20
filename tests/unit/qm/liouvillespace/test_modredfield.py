# -*- coding: utf-8 -*-


import unittest
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.qm.RedfieldRelaxationTensor and 
                 quantarhei.qm.RefieldRateMatrix classes


*******************************************************************************
"""


from quantarhei.qm import ModifiedRedfieldRateMatrix
from quantarhei.qm import ModRedfieldRelaxationTensor
from quantarhei import Hamiltonian
from quantarhei import CorrelationFunction
from quantarhei import Molecule
from quantarhei import Aggregate
from quantarhei.qm.corfunctions import CorrelationFunctionMatrix 
from quantarhei.qm.corfunctions import SpectralDensity
from quantarhei import TimeAxis, eigenbasis_of
from quantarhei.qm import Operator
from quantarhei.qm import SystemBathInteraction
from quantarhei.qm import ReducedDensityMatrixPropagator
from quantarhei.qm import ReducedDensityMatrix
from quantarhei import energy_units
from quantarhei import REAL

class TestModRedfield(unittest.TestCase):
    """Tests for the RedfieldRelaxationTensor class
    
    
    """
    
    def setUp(self,verbose=False):
        
        self.verbose = verbose
        
        time = TimeAxis(0.0,1000,1.0)
        self.time = time
        with energy_units("1/cm"):
            params = {"ftype":"OverdampedBrownian",
                      "reorg":30.0,
                      "T":300.0,
                      "cortime":100.0}
                      
            cf1 = CorrelationFunction(time,params)
            cf2 = CorrelationFunction(time,params)
            sd = SpectralDensity(time,params)
           
            
            m1 = Molecule([0.0, 10000.0])
            m1.set_transition_environment((0,1), cf1)
            m2 = Molecule([0.0, 10000.0])
            m2.set_transition_environment((0,1), cf2)
            
            agg = Aggregate(molecules=[m1, m2])
            
            agg.set_resonance_coupling(0,1, 100.0)
            
        agg.build()
        
        self.H1 = agg.get_Hamiltonian()
        self.sbi1 = agg.get_SystemBathInteraction()
            
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
        KT = numpy.zeros((dim,dim), dtype=REAL)
        KM = numpy.zeros((dim,dim), dtype=REAL)
        
        if tensor:
            
            RT = ModRedfieldRelaxationTensor(self.H1,self.sbi1)
            
            # rates have to be compared in eigenbasis of the Hamiltonian
            with eigenbasis_of(self.H1):
                for n in range(dim):
                    for m in range(dim):
                        #print(n,m,numpy.real(RT.data[n,n,m,m]))
                        KT[n,m] = numpy.real(RT.data[n,n,m,m])
                    
                    
        if matrix:
            
            RR = ModifiedRedfieldRateMatrix(self.H1,self.sbi1) #,self.time)
            
            for n in range(dim):
                for m in range(dim):
                    #print(n,m,numpy.real(RR.data[n,m]))    
                    KM[n,m] = numpy.real(RR.data[n,m])
                    
        
        numpy.testing.assert_allclose(KT,KM, rtol=1.0e-10)
                

    def test_propagation_in_different_basis(self):
        """(REDFIELD) Testing comparison of propagations in different bases

        """

        #LT1 = ModRedfieldRelaxationTensor(self.H1, self.sbi1, as_operators=True)
        LT2 = ModRedfieldRelaxationTensor(self.H1, self.sbi1, as_operators=False)
        
        time = TimeAxis(0.0, 1000, 1.0)
        
        #prop1 = ReducedDensityMatrixPropagator(time, self.H1, LT1)
        prop2 = ReducedDensityMatrixPropagator(time, self.H1, LT2)
        
        rho0 = ReducedDensityMatrix(dim=self.H1.dim)
        rho0.data[1,1] = 1.0
          
        #with eigenbasis_of(self.H1):
        ##if True:
        #    rhot1_e = prop1.propagate(rho0)
            
        with eigenbasis_of(self.H1):
            rhot2_e = prop2.propagate(rho0)

        #rhot1_l = prop1.propagate(rho0)
        rhot2_l = prop2.propagate(rho0)
            
        #numpy.testing.assert_allclose(rhot1_l.data, rhot1_e.data,
        #                              rtol=1.0e-5, atol=1.0e-12)
        #numpy.testing.assert_allclose(rhot2_l.data, rhot1_e.data,
        #                              rtol=1.0e-5, atol=1.0e-12)
        #numpy.testing.assert_allclose(rhot1_e.data, rhot2_e.data,
        #                              rtol=1.0e-5, atol=1.0e-12)              
        numpy.testing.assert_allclose(rhot2_l.data, rhot2_e.data,
                                      rtol=1.0e-5, atol=1.0e-12)

              
