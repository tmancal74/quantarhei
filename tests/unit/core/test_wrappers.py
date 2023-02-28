# -*- coding: utf-8 -*-
import unittest
import numpy

from quantarhei.core.wrappers import prevent_basis_context
from quantarhei.core.wrappers import enforce_basis_context
from quantarhei.core.wrappers import prevent_energy_units_context
from quantarhei.core.wrappers import enforce_energy_units_context
from quantarhei.core.managers import eigenbasis_of
from quantarhei.core.managers import energy_units
from quantarhei.core.managers import Manager

from quantarhei import Hamiltonian
from quantarhei import CorrelationFunction
from quantarhei import TimeAxis

class TestWrappers(unittest.TestCase):
    """Test for function wrappers
    
    
    """
    
    def test_eigenbasis_of_prevention(self):
        """(Wrappers) Testing eigenbasis_of context prevention """
        
        
        if not Manager()._enforce_contexts:
            return
        
        @prevent_basis_context
        def func(text):
            return "Ahoj"
        
        self.assertEqual(func("Ciao"), "Ahoj")
        
        HH = Hamiltonian(data=[[0.0, 0.0],[0.0, 1.0]])
        
        with self.assertRaises(Exception) as context:
            with eigenbasis_of(HH):
                func("Ham")
                
        self.assertTrue("This function MUST NOT be called "
                        in str(context.exception))
            
        
    def test_eigenbasis_of_enforcement(self):
        """(Wrappers) Testing eigenbasis_of context enforcement """
        
        if not Manager()._enforce_contexts:
            return
   
        @enforce_basis_context
        def func(text):
            return "Ahoj"
        
        
        HH = Hamiltonian(data=[[0.0, 0.0],[0.0, 1.0]])
        
        with eigenbasis_of(HH):
            self.assertEqual(func("Ciao"), "Ahoj")
        
        with self.assertRaises(Exception) as context:
            func("Ham")
                
        self.assertTrue("This function MUST be called "
                        in str(context.exception))    


    def test_units_context_prevention(self):
        """(Wrappers) Testing energy_units context prevention """
        
        
        if not Manager()._enforce_contexts:
            return
        
        @prevent_energy_units_context
        def func(text):
            return "Ahoj"
        
        self.assertEqual(func("Ciao"), "Ahoj")
        
        with self.assertRaises(Exception) as context:
            with energy_units("1/cm"):
                func("Ham")
                
        self.assertTrue("This function MUST NOT be called "
                        in str(context.exception))  
        
    def test_units_context_enforcement(self):
        """(Wrappers) Testing energy_units context enforcement """
        
        if not Manager()._enforce_contexts:
            return
   
        @enforce_energy_units_context
        def func(text):
            return "Ahoj"
 
        with energy_units("1/cm"):
            self.assertEqual(func("Ciao"), "Ahoj")
        
        with self.assertRaises(Exception) as context:
            func("Ham")
                
        self.assertTrue("This function MUST be called "
                        in str(context.exception)) 
        
        
    def test_units_enforcement_corrfce(self):
        """(Wrappers) Testing enforcement of units in CorrelationFunction
        
        """
        params = dict(ftype="OverdampedBrownian",
                       reorg = 30.0,
                       cortime = 100.0,
                       T = 300.0)
        
        time = TimeAxis(0.0, 1000, 1.0)
        
        with self.assertRaises(Exception) as context:
            cf = CorrelationFunction(time, params)       
            
        self.assertTrue("This function MUST be called "
                        in str(context.exception)) 

        with energy_units("1/cm"):
            cf = CorrelationFunction(time, params) 

        self.assertTrue(numpy.abs(cf.at(100.0)) > 0.0)

        