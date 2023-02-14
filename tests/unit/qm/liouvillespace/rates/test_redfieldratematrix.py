# -*- coding: utf-8 -*-

import unittest
import numpy

from quantarhei import TimeAxis
from quantarhei import energy_units
from quantarhei import CorrelationFunction
from quantarhei.qm.corfunctions import CorrelationFunctionMatrix 
from quantarhei.qm import SystemBathInteraction
from quantarhei import Hamiltonian
from quantarhei.qm import Operator
from quantarhei.qm import RedfieldRateMatrix
from quantarhei import REAL

class TestRedfieldRateMatrix(unittest.TestCase):
    """Tests for the RateMatrix class
    
    
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

                    
        cm1 = CorrelationFunctionMatrix(time,2,1)
        cm1.set_correlation_function(cf1,[(0,0),(1,1)])
 
          
        K11 = numpy.array([[1.0, 0.0],[0.0, 0.0]],dtype=REAL)
        K21 = numpy.array([[1.0, 0.0],[0.0, 0.0]],dtype=REAL)

            
        KK11 = Operator(data=K11)
        KK21 = Operator(data=K21)
 
        
        self.sbi1 = SystemBathInteraction([KK11,KK21],cm1) 
        
        with energy_units("1/cm"):
            h1 = [[0.0, 100.0],[100.0, 0.0]]
            self.H1 = Hamiltonian(data=h1)
            
        self.sbi2 = SystemBathInteraction()
       
        
        
    def test_create_of_RedfieldRateMatrix(self):
        """(RedfieldRateMatrix) Testing creation 
        
        """
        op = Operator(dim=self.H1.dim)
        
        # testing wrong Hamiltonian
        with self.assertRaises(Exception) as context:
            RR = RedfieldRateMatrix(op,self.sbi1)
            
        self.assertTrue("First argument must be a Hamiltonian"
                        in str(context.exception))
        
        # testing wrong SBI
        with self.assertRaises(Exception) as context:
            RR = RedfieldRateMatrix(self.H1,op)

        self.assertTrue("Second argument must be a SystemBathInteraction"
                        in str(context.exception))
        
        # testing cutoff time
        RR = RedfieldRateMatrix(self.H1, self.sbi1, cutoff_time=100.0)
        
        self.assertTrue(RR._has_cutoff_time)
        self.assertTrue(RR.cutoff_time == 100.0)
        
        # testing empty SBI
        with self.assertRaises(Exception) as context:
            RR = RedfieldRateMatrix(self.H1, self.sbi2)
            
        self.assertTrue("No system bath intraction components present"
                        in str(context.exception))
        
        
        