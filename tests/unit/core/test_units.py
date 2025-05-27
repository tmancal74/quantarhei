# -*- coding: utf-8 -*-

import unittest

"""
*******************************************************************************


    Tests of the quantarhei.core.units package


*******************************************************************************
"""

from quantarhei.core.units import conversion_facs_energy
import scipy.constants as const
from quantarhei.core.units import eps0_int

class TestUnits(unittest.TestCase):
    """Tests for the units package
    
    
    """
    
    def setUp(self,verbose=False):
        self.verbose = verbose
        self.J2eV = 1.0/const.e
        self.J2THz = (1.0/(1.0e12*const.h))
        self.eV2cm = 8065.54413 
        self.THz2cm = 33.35641
     
    
    def test_energy_units_conversions(self):
        """Testing energy units conversion
        
        
        """
        

        
        # J to eV
        eV2int = conversion_facs_energy["eV"]
        J2int = conversion_facs_energy["J"]
        
        J2eV = J2int/eV2int
        
        if self.verbose:
            print("J to eV conversion factor: %s - %s " %
              (J2eV,self.J2eV))
              
        self.assertAlmostEqual((J2int/eV2int)*1.0e-18,self.J2eV*1.0e-18)
        
        # J to THz
        THz2int = conversion_facs_energy["THz"]
        
        if self.verbose:
            print("J to THz conversion factor: %s - %s" %
              ((J2int/THz2int),self.J2THz))
              
        self.assertAlmostEqual((J2int/THz2int)*1.0e-21,self.J2THz*1.0e-21)
        
        # eV to cm-1
        cm2int = conversion_facs_energy["1/cm"]
        
        if self.verbose:
            print("cm to eV conversion factor: %s - %s " %
              (eV2int/cm2int,self.eV2cm))        
                      
        self.assertAlmostEqual((eV2int/cm2int)/8000.0,self.eV2cm/8000.0)

        # THz to cm-1
        THz2cm = THz2int/cm2int
        if self.verbose:
            print("THz to cm conversion factor: %s - %s " %
              (THz2cm,self.THz2cm))        
                      
        self.assertAlmostEqual(THz2cm/33.0,self.THz2cm/33.0)

            
        E_int = (1.0/(4.0*const.pi*eps0_int))
        E_eV = E_int/eV2int
        
        oneDebye = 1.0e-21/const.c
        R3 = (1.0e-10)**3
        E_J = (1.0/(4.0*const.pi*const.epsilon_0))*(oneDebye**2)/R3
        
        if self.verbose:
            print("epsilon_0 [int] = ", eps0_int)
            
        self.assertAlmostEqual(E_eV,E_J*J2eV)
        
        print(100.0*self.THz2cm)
        
       
        