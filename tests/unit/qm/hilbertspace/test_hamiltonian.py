# -*- coding: utf-8 -*-

import unittest
import numpy
import scipy.constants as const

"""
*******************************************************************************


    Tests of the quantarhei.Hamiltonian class


*******************************************************************************
"""


from quantarhei import Hamiltonian
from quantarhei import energy_units

class TestHamiltonian(unittest.TestCase):
    """Tests for the Hamiltonian class
    
    
    """
    
    def setUp(self,verbose=False):
        self.h = [[0.0, 1.0],[1.0, 2.0]]
        self.H = Hamiltonian(data=self.h)
        
        self.verbose = verbose

    
    def test_units_management(self):
        """Testing that Hamiltonian is EnergyUnitsManaged
        
        
        """        
        cm2int = 2.0*const.pi*const.c*1.0e-13
        self.assertEquals(self.H.unit_repr(),"1/fs")
        
        if self.verbose:
            print("In internal:")
            print(self.H.data)
        
        with energy_units("1/cm"):
            self.assertEquals(self.H.unit_repr(),"1/cm")
            self.assertEquals(self.H.data[0,1],1.0/cm2int)
            if self.verbose:
                print("In 1/cm:")
                print(self.H.data)
                
            H2 = Hamiltonian(data=[[1000.0, 0.0],[0.0, 2000.0]])
            
            self.assertAlmostEqual(H2.data[0,0],1000.0)
            self.assertAlmostEqual(H2.data[1,1],2000.0)
            
        self.assertAlmostEqual(H2.data[0,0],1000.0*cm2int)
        self.assertAlmostEqual(H2.data[1,1],2000.0*cm2int)
            
        if self.verbose:
            print("In internal:")
            print(self.H.data)        
        