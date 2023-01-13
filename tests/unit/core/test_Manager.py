# -*- coding: utf-8 -*-

import unittest

"""
*******************************************************************************


    Tests of the quantarhei.Manager class


*******************************************************************************
"""

from quantarhei import Manager

class TestManager(unittest.TestCase):
    """Tests for the Manager class
    
    
    """
    
    def setUp(self):
        pass
     
    
    def test_that_Manager_is_a_singleton(self):
        """Testing that Manager object is a singleton
        
        
        """
        m = Manager()
        n = Manager()
    
        if m is not n:
            raise Exception()
        
        
    
        