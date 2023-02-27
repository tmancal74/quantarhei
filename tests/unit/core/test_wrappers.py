# -*- coding: utf-8 -*-
import unittest

from quantarhei.core.wrappers import prevent_basis_context
from quantarhei.core.wrappers import enforce_basis_context
from quantarhei.core.managers import eigenbasis_of
from quantarhei.core.managers import Manager

from quantarhei import Hamiltonian

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