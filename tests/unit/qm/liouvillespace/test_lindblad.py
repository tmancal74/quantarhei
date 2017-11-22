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

from quantarhei.qm import Operator
from quantarhei.qm import SystemBathInteraction
from quantarhei.qm import ReducedDensityMatrixPropagator
from quantarhei.qm import ReducedDensityMatrix

from quantarhei import Hamiltonian
from quantarhei import energy_units
from quantarhei import TimeAxis


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
                
        #print("KT = ", KT)
        #print("KM = ", KM)
        #print("LT = ", LT.data)
        
        numpy.testing.assert_allclose(KT,KM, rtol=1.0e-2)
           
        
from quantarhei import Molecule, Aggregate, Mode
from quantarhei.qm import ProjectionOperator

class TestElectronicLindblad(unittest.TestCase):
    """Tests for the ElectronicLindbladForm class
    
    
    """
    
    def setUp(self,verbose=False):
        
        self.verbose = verbose
        
        #
        # PURE ELECTRONIC AGGREGATE
        #
    
        m1 = Molecule([0.0, 1.0])
        m2 = Molecule([0.0, 1.0])
        
        agg = Aggregate(molecules=[m1, m2])
        
        agg.set_resonance_coupling(0,1, 0.1)
        
        agg.build()
        
        self.ham = agg.get_Hamiltonian()
                
        KK12 = ProjectionOperator(1, 2, dim=self.ham.dim)
        KK21 = ProjectionOperator(2, 1, dim=self.ham.dim)
        
        self.rates = (1.0/100.0, 1.0/200.0)
        
        self.sbi = SystemBathInteraction([KK12,KK21],
                                          rates=self.rates)
        self.sbi.set_system(agg)

        #
        #  VIBRONIC AGGREGATE
        #
        
        vm1 = Molecule([0.0, 1.0])
        vm2 = Molecule([0.0, 1.0])
        mod1 = Mode(0.01)
        mod2 = Mode(0.01)
        vm1.add_Mode(mod1)
        vm2.add_Mode(mod2)
                
        mod1.set_nmax(0, 3)
        mod1.set_nmax(1, 3)
        
        mod2.set_nmax(0, 3)
        mod2.set_nmax(1, 3)
        
        vagg = Aggregate(molecules=[vm1, vm2])
        vagg.set_resonance_coupling(0, 1, 0.1)
        
        vagg.build()
        
        self.vham = vagg.get_Hamiltonian()
        
        self.vsbi = SystemBathInteraction([KK12, KK21], rates=self.rates)
        self.vsbi.set_system(vagg)       
        
        
    def test_comparison_of_rates_in_electronic(self):
        """Testing that ElectronicLindblad and rate matrix are compatible
        
        This is for the case that everything is electronic
        
        
        """ 
        
        dim = self.ham.dim
        KT = numpy.zeros((dim,dim), dtype=numpy.float64)
        KM = numpy.zeros((dim,dim), dtype=numpy.float64)
        
        LT = ElectronicLindbladForm(self.ham, self.sbi, as_operators=False)
        
        for n in range(dim):
            for m in range(dim):
                KT[n,m] = numpy.real(LT.data[n,n,m,m])
                    
        KM = numpy.zeros((dim,dim))
        KM[1,1] = -self.rates[1]
        KM[2,2] = -self.rates[0]
        KM[1,2] = self.rates[0]
        KM[2,1] = self.rates[1]
                
        #print("KT = ", KT)
        #print("KM = ", KM)
        #print("LT = ", LT.data)
        
        numpy.testing.assert_allclose(KT,KM, rtol=1.0e-2)               


    def test_construction_electronic_lindblad(self):
        """Testing construction of electronic Lindblad form for vibronic system
         
        """
         
        LT = ElectronicLindbladForm(self.vham, self.vsbi, as_operators=False)
        
        dim = self.vham.dim
        zer = numpy.zeros((dim, dim, dim, dim), dtype=numpy.float64)
        numpy.testing.assert_array_equal(LT._data-LT._data, zer)


    def test_comparison_with_and_without_vib(self):
        """Testing ElectronicLindbladForm in propagation
        
        """
        
        # Lindblad forms
        LT = ElectronicLindbladForm(self.ham, self.sbi, as_operators=False)
        vLT = ElectronicLindbladForm(self.vham, self.vsbi, as_operators=False)
       
        # Propagators
        time = TimeAxis(0.0, 1000, 1.0)
        el_prop = ReducedDensityMatrixPropagator(time, self.ham, LT)
        vib_prop = ReducedDensityMatrixPropagator(time, self.vham, vLT)
        
        # electronic initial condition
        rho0 = ReducedDensityMatrix(dim=self.ham.dim)
        rho0.data[2,2] = 1.0
        
        # vibronic intial condition
        vrho0 = ReducedDensityMatrix(dim=self.vham.dim)
        agg = self.vsbi.system
        Nin = agg.vibindices[2][0]
        vrho0.data[Nin, Nin] = 1.0
        
        # propagations
        rhot = el_prop.propagate(rho0)
        vrhot = vib_prop.propagate(vrho0)
        
        # averaging over vibrational level at t=0
        Nel = agg.Nel
        
        aver_vrhot0 = agg.trace_over_vibrations(vrho0) 
                
        # Test
        numpy.testing.assert_array_equal(rhot._data[0,:,:], rho0._data)
        numpy.testing.assert_array_equal(rhot._data[0,:,:], aver_vrhot0._data)
        
        aver_vrhot10 = agg.trace_over_vibrations(vrhot, 10)
        numpy.testing.assert_array_equal(rhot._data[10,:,:], 
                                         aver_vrhot10._data)        
     
        aver_vrhot800 = agg.trace_over_vibrations(vrhot, 800)
        numpy.testing.assert_array_equal(rhot._data[800,:,:], 
                                         aver_vrhot800._data)        
        