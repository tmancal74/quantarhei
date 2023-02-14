# -*- coding: utf-8 -*-


import unittest
#import numpy

from quantarhei import TimeAxis
from quantarhei import energy_units
from quantarhei import eigenbasis_of
from quantarhei import CorrelationFunction
#from quantarhei.qm.corfunctions import CorrelationFunctionMatrix 
from quantarhei.qm import SystemBathInteraction
from quantarhei import Hamiltonian
from quantarhei.qm import Operator
from quantarhei.qm import ModifiedRedfieldRateMatrix as RedfieldRateMatrix
from quantarhei import Molecule
from quantarhei import Aggregate
from quantarhei import REAL

class TestModifiedRedfieldRateMatrix(unittest.TestCase):
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

 
        m1 = Molecule([0.0, 1.0])
        m2 = Molecule([0.0, 1.0])
          
        agg = Aggregate(molecules=[m1, m2])
        m1.set_transition_environment((0,1), cf1)
        m2.set_transition_environment((0,1), cf1)

        with energy_units("1/cm"):        
            agg.set_resonance_coupling(0,1, 100.0)
            
        agg.build(mult=1)
        
        self.agg = agg

        self.sbi1 = agg.get_SystemBathInteraction() 
        
        
        with energy_units("1/cm"):
            h1 = [[0.0, 100.0],[100.0, 0.0]]
            self.H1 = Hamiltonian(data=h1)
            
        self.sbi2 = SystemBathInteraction()
        
        
       
        
        
    def test_create_ModifiedRedfieldRateMatrix(self):
        """(ModifiedRedfieldRateMatrix) Testing creation 
        
        """
        op = Operator(dim=self.H1.dim)
        
        # testing wrong Hamiltonian
        with self.assertRaises(Exception) as context:
            RR = RedfieldRateMatrix(op,self.sbi1,self.sbi1.TimeAxis)
            
        self.assertTrue("First argument must be a Hamiltonian"
                        in str(context.exception))
        
        
        # testing wrong SBI
        with self.assertRaises(Exception) as context:
            RR = RedfieldRateMatrix(self.H1,op,self.sbi1.TimeAxis)

        #self.assertTrue("Second argument must be a SystemBathInteraction"
        #                in str(context.exception))
        
        # testing initialization
        RR = RedfieldRateMatrix(self.H1, self.sbi1, self.sbi1.TimeAxis,
                                initialize=False, cutoff_time=100.0)
        
        self.assertTrue(RR._has_cutoff_time)
        self.assertTrue(RR.cutoff_time == 100.0)
        
        # testing that calculation works
        HH = self.agg.get_Hamiltonian()
        sbi = self.agg.get_SystemBathInteraction()
        HH.protect_basis()
        with eigenbasis_of(HH):
            RR = RedfieldRateMatrix(HH, sbi, sbi.TimeAxis)  
        HH.unprotect_basis()
        
        