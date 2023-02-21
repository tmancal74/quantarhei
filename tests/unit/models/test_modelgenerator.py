# -*- coding: utf-8 -*-

import unittest
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.models.ModelGenerator class


*******************************************************************************
"""

#legacy = False
#import tempfile
from quantarhei.models.modelgenerator import ModelGenerator
from quantarhei import TimeAxis
        



class TestModelGenerator(unittest.TestCase):
    """(ModelGenerator) Tests for the model generator
    
    
    """
    
    def setUp(self):
        
        self.mgen = ModelGenerator()
        
        
        
        
    def testing_model_system_creation(self):
        """(ModelGenerator) Testing generation of test aggregates
        
        """

        mgen = self.mgen
        
        agg1 = mgen.get_Aggregate("dimer-1")
        agg1.build(mult=1)
        H = agg1.get_Hamiltonian()
        self.assertEqual(H.dim, 3)
        
        agg1 = mgen.get_Aggregate("dimer-1")
        agg1.build(mult=2)
        H = agg1.get_Hamiltonian()
        self.assertEqual(H.dim, 4) 
        
        agg1 = mgen.get_Aggregate("trimer-1")
        agg1.build(mult=1)
        H = agg1.get_Hamiltonian()
        self.assertEqual(H.dim, 4)
        
        agg1 = mgen.get_Aggregate("trimer-1")
        agg1.build(mult=2)
        H = agg1.get_Hamiltonian()
        self.assertEqual(H.dim, 7)    
        
        agg1 = mgen.get_Aggregate("pentamer-1")
        agg1.build(mult=1)
        H = agg1.get_Hamiltonian()
        self.assertEqual(H.dim, 6)
        
        agg1 = mgen.get_Aggregate("pentamer-1")
        agg1.build(mult=2)
        H = agg1.get_Hamiltonian()
        self.assertEqual(H.dim, 16)         
        

        with self.assertRaises(Exception) as context:
        #if True:
            agg1 = mgen.get_Aggregate("pentamer")
            
        self.assertTrue("Unknown model name pentamer" 
                        in str(context.exception))
        

    def testing_model_system_creation_with_env(self):
        """(ModelGenerator) Testing generation of aggregates with environment
        
        """        
        
        mgen = self.mgen
        
        agg1 = mgen.get_Aggregate_with_environment("dimer-1_env")
        agg1.build(mult=1)
        H = agg1.get_Hamiltonian()
        self.assertEqual(H.dim, 3)
        sbi = agg1.get_SystemBathInteraction()
        T = sbi.get_temperature()
        self.assertEqual(T,300)
        
        ta = TimeAxis(0.0, 1000, 1.0)
        agg1 = mgen.get_Aggregate_with_environment("dimer-1_env", ta)
        agg1.build(mult=1)
        H = agg1.get_Hamiltonian()
        self.assertEqual(H.dim, 3)
        sbi = agg1.get_SystemBathInteraction()
        T = sbi.get_temperature()
        self.assertEqual(T,300)        
        
        agg1 = mgen.get_Aggregate_with_environment("dimer-1_env")
        agg1.build(mult=2)
        H = agg1.get_Hamiltonian()
        self.assertEqual(H.dim, 4) 
        sbi = agg1.get_SystemBathInteraction()
        T = sbi.get_temperature()
        self.assertEqual(T,300)
        
        agg1 = mgen.get_Aggregate_with_environment("trimer-1_env")
        agg1.build(mult=1)
        H = agg1.get_Hamiltonian()
        self.assertEqual(H.dim, 4)
        sbi = agg1.get_SystemBathInteraction()
        T = sbi.get_temperature()
        self.assertEqual(T,300)
        
        agg1 = mgen.get_Aggregate_with_environment("trimer-1_env")
        agg1.build(mult=2)
        H = agg1.get_Hamiltonian()
        self.assertEqual(H.dim, 7)    
        sbi = agg1.get_SystemBathInteraction()
        T = sbi.get_temperature()
        self.assertEqual(T,300)
        
        agg1 = mgen.get_Aggregate_with_environment("pentamer-1_env")
        agg1.build(mult=1)
        H = agg1.get_Hamiltonian()
        self.assertEqual(H.dim, 6)
        sbi = agg1.get_SystemBathInteraction()
        T = sbi.get_temperature()
        self.assertEqual(T,300)
        
        agg1 = mgen.get_Aggregate_with_environment("pentamer-1_env")
        agg1.build(mult=2)
        H = agg1.get_Hamiltonian()
        self.assertEqual(H.dim, 16)    
        sbi = agg1.get_SystemBathInteraction()
        T = sbi.get_temperature()
        self.assertEqual(T,300)

        agg1 = mgen.get_Aggregate_with_environment("pentamer-1_env", ta)
        agg1.build(mult=2)
        H = agg1.get_Hamiltonian()
        self.assertEqual(H.dim, 16)    
        sbi = agg1.get_SystemBathInteraction()
        T = sbi.get_temperature()
        self.assertEqual(T,300)
        
        with self.assertRaises(Exception) as context:
        #if True:
            agg1 = mgen.get_Aggregate_with_environment("pentamer")
            
        self.assertTrue("Unknown model name pentamer" 
                        in str(context.exception))        
        
        