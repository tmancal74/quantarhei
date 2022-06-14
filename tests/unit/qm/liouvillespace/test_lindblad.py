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
from quantarhei.qm import ProjectionOperator

from quantarhei import Hamiltonian
from quantarhei import energy_units
from quantarhei import TimeAxis

from quantarhei import eigenbasis_of, Manager


class TestLindblad(unittest.TestCase):
    """Tests for the LindbladForm class
    
    
    """
    
    def setUp(self,verbose=False):
        
        self.verbose = verbose
           
        #
        # Lindblad projection operators
        #
        K12 = numpy.array([[0.0, 1.0],[0.0, 0.0]],dtype=numpy.float)
        K21 = numpy.array([[0.0, 0.0],[1.0, 0.0]],dtype=numpy.float)

        KK12 = Operator(data=K12)
        KK21 = Operator(data=K21)

        self.KK12 = KK12
        self.KK21 = KK21
        
        #
        # Linbdlad rates
        #
        self.rates = (1.0/100.0, 1.0/200.0)
        
        #
        # System-bath interaction using operators and rates in site basis
        #
        self.sbi1 = SystemBathInteraction([KK12,KK21],
                                          rates=self.rates)
        self.sbi2 = SystemBathInteraction([KK12,KK21], 
                                          rates=self.rates)
        
        #
        # Test Hamiltonians
        #
        with energy_units("1/cm"):
            h1 = [[100.0, 0.0],[0.0, 0.0]]
            h2 = [[100.0, 0.0],[0.0, 0.0]]
            self.H1 = Hamiltonian(data=h1)
            self.H2 = Hamiltonian(data=h2)
            
            h3 = [[100.0, 20.0],[20.0, 0.0]]
            self.H3 = Hamiltonian(data=h3)
            
            # less trivial Hamiltonian
            h4 = [[100.0, 200.0, 30.0  ],
                  [200.0, 50.0,  -100.0],
                  [30.0, -100.0,  0.0 ]]
            self.H4 = Hamiltonian(data=h4)
            
            h4s = [[100.0, 0.0, 0.0  ],
                  [0.0, 50.0,  0.0],
                  [0.0, 0.0,  0.0 ]]
            
            self.H4s = Hamiltonian(data=h4s)
            
            
        #
        # Projection operators in eigenstate basis
        #
        with eigenbasis_of(self.H3):
            K_12 = ProjectionOperator(0, 1, dim=2)
            K_21 = ProjectionOperator(1, 0, dim=2)
            self.K_12 = K_12
            self.K_21 = K_21
            
        with eigenbasis_of(self.H4):
            Ke_12 = ProjectionOperator(0, 1, dim=3)
            Ke_21 = ProjectionOperator(1, 0, dim=3)
            Ke_23 = ProjectionOperator(1, 2, dim=3)
            Ke_32 = ProjectionOperator(2, 1, dim=3)

        Ks_12 = ProjectionOperator(0, 1, dim=3)
        Ks_21 = ProjectionOperator(1, 0, dim=3)
        Ks_23 = ProjectionOperator(1, 2, dim=3)
        Ks_32 = ProjectionOperator(2, 1, dim=3) 
           
        self.rates4 = [1.0/100, 1.0/200, 1.0/150, 1.0/300]
            
        # 
        # System-bath operators defined in exciton basis
        #
        self.sbi3 = SystemBathInteraction([K_12, K_21],
                                          rates=self.rates)
        
        self.sbi4e = SystemBathInteraction([Ke_12, Ke_21, Ke_23, Ke_32],
                                          rates=self.rates4)
        self.sbi4s = SystemBathInteraction([Ks_12, Ks_21, Ks_23, Ks_32],
                                          rates=self.rates4)
        
        

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
                
        
        numpy.testing.assert_allclose(KT,KM, rtol=1.0e-2)
     
        
    def test_comparison_of_dynamics(self):
        """Testing site basis dynamics by Lindblad

        """

        LT1 = LindbladForm(self.H1, self.sbi1, as_operators=True)
        LT2 = LindbladForm(self.H1, self.sbi1, as_operators=False)
        
        time = TimeAxis(0.0, 1000, 1.0)
        
        prop1 = ReducedDensityMatrixPropagator(time, self.H1, LT1)
        prop2 = ReducedDensityMatrixPropagator(time, self.H1, LT2)
        
        rho0 = ReducedDensityMatrix(dim=self.H1.dim)
        rho0.data[1,1] = 1.0
          
        rhot1 = prop1.propagate(rho0)
        rhot2 = prop2.propagate(rho0)
        
        numpy.testing.assert_allclose(rhot1.data,rhot2.data) #, rtol=1.0e-2) 


    def test_propagation_in_different_basis(self):
        """(LINDBLAD) Testing comparison of propagations in different bases

        """

        LT1 = LindbladForm(self.H1, self.sbi1, as_operators=True)
        LT2 = LindbladForm(self.H1, self.sbi1, as_operators=False)
        
        time = TimeAxis(0.0, 1000, 1.0)
        
        prop1 = ReducedDensityMatrixPropagator(time, self.H1, LT1)
        prop2 = ReducedDensityMatrixPropagator(time, self.H1, LT2)
        
        rho0 = ReducedDensityMatrix(dim=self.H1.dim)
        rho0.data[1,1] = 1.0
          
        with eigenbasis_of(self.H1):
            rhot1_e = prop1.propagate(rho0)
            
        with eigenbasis_of(self.H1):
            rhot2_e = prop2.propagate(rho0)

        rhot1_l = prop1.propagate(rho0)
        rhot2_l = prop2.propagate(rho0)
            
        numpy.testing.assert_allclose(rhot1_l.data, rhot1_e.data)
        numpy.testing.assert_allclose(rhot2_l.data, rhot2_e.data)
        numpy.testing.assert_allclose(rhot1_e.data, rhot2_e.data) #, rtol=1.0e-2) 
        
    
    def test_transformation_in_different_basis(self):
        """(LINDBLAD) Testing transformations into different bases

        """
        #Manager().warn_about_basis_change = True
        #Manager().warn_about_basis_changing_objects = True
        
        LT1 = LindbladForm(self.H1, self.sbi1, as_operators=True, name="LT1")
        LT2 = LindbladForm(self.H1, self.sbi1, as_operators=False, name="LT2")

        rho0 = ReducedDensityMatrix(dim=self.H1.dim, name="ahoj")
        with eigenbasis_of(self.H1):
            rho0.data[1,1] = 0.7
            rho0.data[0,0] = 0.3
                
          
        with eigenbasis_of(self.H1):
            rhot1_e = LT1.apply(rho0, copy=True)
            
        with eigenbasis_of(self.H1):
            rhot2_e = LT2.apply(rho0, copy=True)

        rhot1_l = LT1.apply(rho0, copy=True)
        rhot2_l = LT2.apply(rho0, copy=True)
            
        numpy.testing.assert_allclose(rhot1_l.data, rhot1_e.data)
        numpy.testing.assert_allclose(rhot2_l.data, rhot2_e.data)
        numpy.testing.assert_allclose(rhot1_e.data, rhot2_e.data) #, rtol=1.0e-2) 

    

    def test_comparison_of_exciton_dynamics(self):
        """Testing exciton basis dynamics by Lindblad

        """
        
        # site basis form to be compared with
        LT1 = LindbladForm(self.H1, self.sbi1, as_operators=True)
        
        # exciton basis forms
        LT13 = LindbladForm(self.H3, self.sbi3, as_operators=True)
        LT23 = LindbladForm(self.H3, self.sbi3, as_operators=False)

        LT4e = LindbladForm(self.H4, self.sbi4e, as_operators=True)
        LT4s = LindbladForm(self.H4s, self.sbi4s, as_operators=True)
        
        time = TimeAxis(0.0, 1000, 1.0)
        
        #
        # Propagators
        #
        prop0 = ReducedDensityMatrixPropagator(time, self.H1, LT1)
        prop1 = ReducedDensityMatrixPropagator(time, self.H3, LT13)
        prop2 = ReducedDensityMatrixPropagator(time, self.H3, LT23)
        prop4e = ReducedDensityMatrixPropagator(time, self.H4, LT4e)
        prop4s = ReducedDensityMatrixPropagator(time, self.H4s, LT4s)        

        # 
        # Initial conditions
        #
        rho0 = ReducedDensityMatrix(dim=self.H3.dim)
        rho0c = ReducedDensityMatrix(dim=self.H1.dim) # excitonic
        with eigenbasis_of(self.H3):
            rho0c.data[1,1] = 1.0
        rho0.data[1,1] = 1.0
        
        rho04e = ReducedDensityMatrix(dim=self.H4.dim)
        rho04s = ReducedDensityMatrix(dim=self.H4.dim)
        with eigenbasis_of(self.H4):
            rho04e.data[2,2] = 1.0
        rho04s.data[2,2] = 1.0        
          
        #
        # Propagations
        #
        rhotc = prop0.propagate(rho0c)
        rhot1 = prop1.propagate(rho0)
        rhot2 = prop2.propagate(rho0)
        
        rhot4e = prop4e.propagate(rho04e)
        rhot4s = prop4s.propagate(rho04s)
        
        # propagation with operator- and tensor forms should be the same
        numpy.testing.assert_allclose(rhot1.data,rhot2.data) #, rtol=1.0e-2) 

        #
        # Population time evolution by Lindblad is independent 
        # of the level structure and basis, as long as I compare
        # populations in basis in which the Lindblad form was defined
        #

        P = numpy.zeros((2, time.length))
        Pc = numpy.zeros((2, time.length))
        
        P4e = numpy.zeros((3, time.length))
        P4s = numpy.zeros((3, time.length))
        
        with eigenbasis_of(self.H3):
            for i in range(time.length):
                P[0,i] = numpy.real(rhot1.data[i,0,0])  # population of exciton 0
                P[1,i] = numpy.real(rhot1.data[i,1,1])  # population of exciton 1
                

        for i in range(time.length):
            Pc[0,i] = numpy.real(rhotc.data[i,0,0])  # population of exciton 0
            Pc[1,i] = numpy.real(rhotc.data[i,1,1])  # population of exciton 1
        
        # we compare populations
        numpy.testing.assert_allclose(Pc,P) #, rtol=1.0e-2) 

        with eigenbasis_of(self.H4):
            for i in range(time.length):
                P4e[0,i] = numpy.real(rhot4e.data[i,0,0])  # population of exciton 0
                P4e[1,i] = numpy.real(rhot4e.data[i,1,1])  # population of exciton 1
                P4e[2,i] = numpy.real(rhot4e.data[i,2,2])  # population of exciton 1

        for i in range(time.length):
            P4s[0,i] = numpy.real(rhot4s.data[i,0,0])  # population of exciton 0
            P4s[1,i] = numpy.real(rhot4s.data[i,1,1])  # population of exciton 1
            P4s[2,i] = numpy.real(rhot4s.data[i,2,2])  # population of exciton 1

        
#        import matplotlib.pyplot as plt
        #
#        plt.plot(time.data,P4s[0,:], "-r")
#        plt.plot(time.data,P4s[1,:], "-r")
#        plt.plot(time.data,P4s[2,:], "-r")
#        plt.plot(time.data,P4e[0,:], "-g")
#        plt.plot(time.data,P4e[1,:], "-g")
#        plt.plot(time.data,P4e[2,:], "-g")
        
#        plt.show()

        numpy.testing.assert_allclose(P4e, P4s, atol=1.0e-8) 
        


from quantarhei import Molecule, Aggregate, Mode
#from quantarhei.qm import ProjectionOperator

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
        numpy.testing.assert_allclose(rhot._data[0,:,:], rho0._data)
        numpy.testing.assert_allclose(rhot._data[0,:,:], aver_vrhot0._data)
        
        aver_vrhot10 = agg.trace_over_vibrations(vrhot, 10)
        numpy.testing.assert_allclose(rhot._data[10,:,:], 
                                         aver_vrhot10._data)        
     
        aver_vrhot800 = agg.trace_over_vibrations(vrhot, 800)
        numpy.testing.assert_allclose(rhot._data[800,:,:], 
                                         aver_vrhot800._data)        
        