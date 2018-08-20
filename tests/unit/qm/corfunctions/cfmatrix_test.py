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

        self.cf1 = qr.CorrelationFunction(self.time, pars1)
        self.cf2 = qr.CorrelationFunction(self.time, pars2)

    
    def test_of_creation(self):
        """(CorrelationFunctionMatrix) Test of creation 
        
        """
        cfm = cors.CorrelationFunctionMatrix(self.time, nob=5)
        
        cfm.set_correlation_function(self.cf1, [(1,1),(2,2),(3,3)])
        cfm.set_correlation_function(self.cf2, [(4,4)])
    
    
    def test_of_temperature_retrieval(self):
        """(CorrelationFunctionMatrix) Test of temperature retrieval 
        """
        cfm = cors.CorrelationFunctionMatrix(self.time, nob=5)
        
        cfm.set_correlation_function(self.cf1, [(1,1),(2,2),(3,3)])
        cfm.set_correlation_function(self.cf2, [(4,4)])
        
        T = cfm.get_temperature()
        
        numpy.testing.assert_almost_equal(T, self.temperature)

                

