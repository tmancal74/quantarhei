# -*- coding: utf-8 -*-

import unittest
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.Molecule class


*******************************************************************************
"""

from quantarhei import Molecule

class TestMolecule(unittest.TestCase):
    """Tests for the Manager class
    
    
    """
    
    def setUp(self):
        self.en = [0.0, 1.0, 2.0]
        self.m = Molecule(name="Molecule",elenergies=self.en)     
    
    def test_Molecule_instantiation(self):
        """Testing Molecule instantiation
        
        
        """

        
        self.assertEqual(self.m.name,"Molecule")
        for i in range(2):
            self.assertEqual(self.m.elenergies[i],self.en[i])
            
            
    def test_get_Hamiltonian(self):
        """Testing that Molecule returns correct Hamiltonian 
        
        """
        
        H = self.m.get_Hamiltonian()

        h = numpy.zeros((3,3),dtype=numpy.float)
        h[1,1] = 1.0
        h[2,2] = 2.0
        
        self.assertTrue(numpy.allclose(H.data,h))
        
        