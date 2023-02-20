# -*- coding: utf-8 -*-


import unittest
import numpy
#import h5py


"""
*******************************************************************************


    Tests of the quantarhei.qm.corfunctions.cfmatrix module


*******************************************************************************
"""


import quantarhei as qr
import quantarhei.qm.corfunctions as cors

class TestCFMatrix(unittest.TestCase):
    """Tests of matrix of corelation functions module
    
    
    """
 
    def setUp(self, verbose=False):
        """Initializes the calculation
    
        """
        
        self.verbose = verbose
        
        self.time = qr.TimeAxis(0, 1000, 1.0)
        self.temperature = 300
        pars1 = dict(ftype="OverdampedBrownian", reorg=30, cortime=100,
                     T=self.temperature)
        pars2 = dict(ftype="OverdampedBrownian", reorg=80, cortime=200,
                     T=self.temperature)

        with qr.energy_units("1/cm"):
            self.cf1 = qr.CorrelationFunction(self.time, pars1)
            self.cf2 = qr.CorrelationFunction(self.time, pars2)

        #
        # We explicitely create a correlation function matrix
        #
        cfm = cors.CorrelationFunctionMatrix(self.time, nob=5)
        
        cfm.set_correlation_function(self.cf1, [(1,1),(2,2),(3,3)])
        cfm.set_correlation_function(self.cf2, [(4,4)])

        self.cfm = cfm
        
        #
        # We create an aggregate and look at its correlation function matrix
        #
        m1 = qr.Molecule([0.0, 1.0])
        m2 = qr.Molecule([0.0, 1.01])
        m1.set_transition_environment((0,1), self.cf1)
        m2.set_transition_environment((0,1), self.cf2)
        
        agg = qr.Aggregate(molecules=[m1, m2])
        agg2 = agg.deepcopy()
        agg.build()
        
        self.agg = agg

        agg2.set_resonance_coupling(0, 1, 0.02)
        agg2.build()
        
        self.agg2 = agg2
        
        
    
    def test_of_creation(self):
        """(CorrelationFunctionMatrix) Test of creation 
        
        """
        cfm = cors.CorrelationFunctionMatrix(self.time, nob=5)
        
        cfm.set_correlation_function(self.cf1, [(1,1),(2,2),(3,3)])
        cfm.set_correlation_function(self.cf2, [(4,4)])
    
    
    def test_of_temperature_retrieval(self):
        """(CorrelationFunctionMatrix) Test of temperature retrieval 
        """
        cfm = self.cfm
        
        T = cfm.get_temperature()
        
        numpy.testing.assert_almost_equal(T, self.temperature)

                
    def test_reorganization_energy(self):
        """(CorrelationFunctionMatrix) Test of reorganization energy retrieval
        """
        
        cfm = self.cfm
        
        lam = cfm.get_reorganization_energy(1, 1)
        self.assertTrue(lam==self.cf1.get_reorganization_energy())
 
        lam = cfm.get_reorganization_energy(4, 4)
        self.assertTrue(lam==self.cf2.get_reorganization_energy())
        
        lam = cfm.get_reorganization_energy(4, 2)
        self.assertTrue(lam==0.0)
 
        with qr.energy_units("1/cm"):
            self.assertTrue(30.0==self.cf1.get_reorganization_energy())
            
            
    def test_correlation_function4(self):
        """(CorrelationFunctionMatrix) Test of generalized correlation functions
        
        """
        
        agg = self.agg
        sbi = agg.get_SystemBathInteraction()
        
        cfm = sbi.CC
        cfm.init_site_mapping()
        
        
        lam = cfm.get_reorganization_energy4(0,0,0,0)
        self.assertTrue(lam==0.0)
        
        lam = cfm.get_reorganization_energy4(1,1,1,1) 
        self.assertTrue(qr.convert(lam,"int","1/cm")==30.0)
        
        lam = cfm.get_reorganization_energy4(2,2,2,2)
        self.assertTrue(qr.convert(lam,"int","1/cm")==80.0)
        
        agg2 = self.agg2
        
        sbi = agg2.get_SystemBathInteraction()
        cfm = sbi.CC
        
        HH = agg2.get_Hamiltonian()
        
        ee, SS = numpy.linalg.eigh(HH._data)
        
        cfm.transform(SS)
        
        lam = cfm.get_reorganization_energy4(2,2,2,2)
        self.assertTrue(qr.convert(lam,"int","1/cm")!=80.0)
        
        
        
