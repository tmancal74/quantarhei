# -*- coding: utf-8 -*-

import unittest
import numpy


from quantarhei.qm import RateMatrix
from quantarhei import REAL

class TestRateMatrix(unittest.TestCase):
    """Tests for the RateMatrix class
    
    
    """
    
    def setUp(self,verbose=False):
        
        self.verbose = verbose
        
        RM = RateMatrix(dim=5)
        
        self.RM = RM
        
        
    def test_rate_matrix_creation(self):
        """(RATEMATRIX) Testing rate matrix creation

        """        
        
        tst = numpy.zeros((5,5), dtype=REAL)
        RM = self.RM.data
        
        numpy.testing.assert_allclose(tst, RM,
                                      rtol=1.0e-5, atol=1.0e-12)
        
        OK = False
        try:
            RM2 = RateMatrix()
        except:
            OK = True
            
        assert OK == True
        
        RM2 = RateMatrix(data=numpy.ones((6,6),dtype=REAL))
        
        assert RM2.data[3,4] == 1.0e0
        
    def test_rate_value_setting(self):
        """(RATEMATRIX) Testing that a rate can be set
        
        """
        RM = self.RM
        
        RM.set_rate((2,3), 12.0)
        
        assert RM.data[3,3] == -12.0
        assert RM.data[2,3] == 12.0
        assert RM.data[3,2] == 0.0
        
        
        
        
        