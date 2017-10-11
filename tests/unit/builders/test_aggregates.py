# -*- coding: utf-8 -*-

import unittest
import h5py
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.Aggregate class


*******************************************************************************
"""

from quantarhei import Aggregate
from quantarhei import Molecule


class TestAggregate(unittest.TestCase):
    """Tests for the Manager class
    
    
    """
    
    def setUp(self):
        m1 = Molecule(name="Molecule 1", elenergies=[0.0, 1.0])
        m2 = Molecule(name="Molecule 2", elenergies=[0.0, 1.0])
        self.agg = Aggregate(name="TestAgg", molecules=[m1,m2])
        
        self.agg.init_coupling_matrix()
    
    def test_saving_of_aggregate(self):
        """Testing the saving capability of the Aggregate class
        
        
        """
        use_temporary_file = True
        
        if use_temporary_file: 
            
            drv = "core"
            bcs = False
        
            with h5py.File('tempfile.hdf5', 
                           driver=drv, 
                           backing_store=bcs) as f:
                                             
                self.agg.save_as(f,"Aggregate")
                
                # reread it
                agg = Aggregate()
                agg.load_as(f,"Aggregate")

        else:

            with h5py.File('tempfile.hdf5') as f:                                           
                self.agg.save_as(f,"Aggregate")
            
            with h5py.File('tempfile.hdf5') as f:
                agg = Aggregate()
                agg.load_as(f,"Aggregate")
            
            
        
        self.assertEqual(self.agg.nmono, agg.nmono)
        self.assertEqual(self.agg.mult, agg.mult)
        self.assertEqual(self.agg.sbi_mult, agg.sbi_mult)
        self.assertEqual(self.agg.name, agg.name)
        self.assertEqual(self.agg._has_egcf_matrix, agg._has_egcf_matrix)
        self.assertEqual(self.agg._has_system_bath_interaction,
                         agg._has_system_bath_interaction)
        self.assertEqual(self.agg.coupling_initiated, agg.coupling_initiated)
        self.assertEqual(self.agg._has_relaxation_tensor,
                         agg._has_relaxation_tensor)
        self.assertEqual(self.agg._relaxation_theory, agg._relaxation_theory)
        self.assertEqual(self.agg._built, agg._built)
        self.assertEqual(self.agg.Nel,agg.Nel)
        self.assertEqual(self.agg.Ntot,agg.Ntot)
        self.assertEqual(self.agg.Nb,agg.Nb)
        for key in agg.mnames.keys():
            self.assertEqual(self.agg.mnames[key],agg.mnames[key])
        numpy.testing.assert_array_equal(self.agg.resonance_coupling, 
                                         agg.resonance_coupling) 

        