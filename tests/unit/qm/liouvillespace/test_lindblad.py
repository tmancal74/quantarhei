# -*- coding: utf-8 -*-

import unittest
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.qm.LindbladForm class


*******************************************************************************
"""

from quantarhei.qm import LindbladForm
from quantarhei.qm import ElectronicLindbladForm
from quantarhei import Hamiltonian
from quantarhei.qm import Operator
from quantarhei.qm import SystemBathInteraction

from quantarhei import energy_units

class TestLindblad(unittest.TestCase):
    """Tests for the LindbladForm class
    
    
    """
    
    def setUp(self,verbose=False):
        
        self.verbose = verbose
                  
        K12 = numpy.array([[0.0, 1.0],[0.0, 0.0]],dtype=numpy.float)
        K21 = numpy.array([[0.0, 0.0],[1.0, 0.0]],dtype=numpy.float)

            
        KK12 = Operator(data=K12)
        KK21 = Operator(data=K21)

        self.KK12 = KK12
        self.KK21 = KK21
        self.rates = (1.0/100.0, 1.0/200.0)
        
        self.sbi1 = SystemBathInteraction([KK12,KK21],
                                          rates=self.rates)
        self.sbi2 = SystemBathInteraction([KK12,KK21], 
                                          rates=self.rates)
        
        with energy_units("1/cm"):
            h1 = [[100.0, 0.0],[0.0, 0.0]]
            h2 = [[100.0, 0.0],[0.0, 0.0]]
            self.H1 = Hamiltonian(data=h1)
            self.H2 = Hamiltonian(data=h2)
            

    def test_comparison_of_rates(self):
        """Testing that Lindblad tensor and rate matrix are compatible
        
        
        """ 
        tensor = True
#        matrix = True
        
        dim = self.H1.dim
        KT = numpy.zeros((dim,dim), dtype=numpy.float64)
        KM = numpy.zeros((dim,dim), dtype=numpy.float64)
        
        if tensor:
            #print(self.H1)
            LT = LindbladForm(self.H1, self.sbi1, as_operators=False)
            
            for n in range(2):
                for m in range(2):
                    #print(n,m,numpy.real(RT.data[n,n,m,m]))
                    KT[n,m] = numpy.real(LT.data[n,n,m,m])
                    
        KM = numpy.zeros((dim,dim))
        KM[0,0] = -self.rates[1]
        KM[1,1] = -self.rates[0]
        KM[0,1] = self.rates[0]
        KM[1,0] = self.rates[1]
                
        print("KT = ", KT)
        print("KM = ", KM)
        print("LT = ", LT.data)
        
        numpy.testing.assert_allclose(KT,KM, rtol=1.0e-2)
           
        
class TestElectronicLindblad(unittest.TestCase):
    """Tests for the ElectronicLindbladForm class
    
    
    """
    
    def setUp(self,verbose=False):
        
        self.verbose = verbose
                  
        K12 = numpy.array([[0.0, 1.0],[0.0, 0.0]],dtype=numpy.float)
        K21 = numpy.array([[0.0, 0.0],[1.0, 0.0]],dtype=numpy.float)

            
        KK12 = Operator(data=K12)
        KK21 = Operator(data=K21)

        self.KK12 = KK12
        self.KK21 = KK21
        self.rates = (1.0/100.0, 1.0/200.0)
        
        self.sbi1 = SystemBathInteraction([KK12,KK21],
                                          rates=self.rates)
        self.sbi2 = SystemBathInteraction([KK12,KK21], 
                                          rates=self.rates)
        
        with energy_units("1/cm"):
            h1 = [[100.0, 0.0],[0.0, 0.0]]
            h2 = [[100.0, 0.0],[0.0, 0.0]]
            self.H1 = Hamiltonian(data=h1)
            self.H2 = Hamiltonian(data=h2)
            

    def test_comparison_of_rates_in_electronic(self):
        """Testing that ElectronicLindblad tensor and rate matrix are compatible
        
        This is for the case that everything is electronic
        
        
        """ 
        tensor = True
#        matrix = True
        
        dim = self.H1.dim
        KT = numpy.zeros((dim,dim), dtype=numpy.float64)
        KM = numpy.zeros((dim,dim), dtype=numpy.float64)
        
        if tensor:
            #print(self.H1)
            LT = ElectronicLindbladForm(self.H1, self.sbi1, as_operators=False)
            
            for n in range(2):
                for m in range(2):
                    #print(n,m,numpy.real(RT.data[n,n,m,m]))
                    KT[n,m] = numpy.real(LT.data[n,n,m,m])
                    
        KM = numpy.zeros((dim,dim))
        KM[0,0] = -self.rates[1]
        KM[1,1] = -self.rates[0]
        KM[0,1] = self.rates[0]
        KM[1,0] = self.rates[1]
                
        print("KT = ", KT)
        print("KM = ", KM)
        print("LT = ", LT.data)
        
        numpy.testing.assert_allclose(KT,KM, rtol=1.0e-2)
                
