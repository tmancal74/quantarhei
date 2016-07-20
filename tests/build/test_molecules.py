# -*- coding: utf-8 -*-

import unittest


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
        pass
     
    
    def test_Molecule_instantiation(self):
        """Testing Molecule instantiation
        
        
        """
        en = [0.0, 1.0, 2.0]
        m = Molecule(name="Molecule",elenergies=en)
        
        self.assertEqual(m.name,"Molecule")
        for i in range(2):
            self.assertEqual(m.elenergies[i],en[i])