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
        
        
        
    def test_rate_matrix_creation(self):
        """(RATEMATRIX) Testing rate matrix creation

        """        
        
        # testing that rate matrix was correctly created from dimension
        tst = numpy.zeros((5,5), dtype=REAL)
        RM = RateMatrix(dim=5)
        RM = RM.data
        
        numpy.testing.assert_allclose(tst, RM,
                                      rtol=1.0e-5, atol=1.0e-12)
        
        # testing that creation without arguments raises exception
        with self.assertRaises(Exception) as context:
            RM2 = RateMatrix()
        self.assertTrue("One of the arguments has to be specified" 
                        in str(context.exception))
        
        # test some data values
        RM2 = RateMatrix(data=numpy.ones((6,6),dtype=REAL))
        
        assert RM2.data[3,4] == 1.0e0
        assert RM2.data[0,0] == 1.0e0
        
        # testing for non-square data refusal
        tst2 = numpy.zeros((2,3), dtype=REAL)
        with self.assertRaises(Exception) as context:
            RM = RateMatrix(data=tst2)
        self.assertTrue("Expecting rectangular matrix" 
                        in str(context.exception))    

        # testing for wrong dimension specification
        with self.assertRaises(Exception) as context:
            RM = RateMatrix(data=tst, dim=6)
        self.assertTrue("Inconsistent data dimensions" 
                        in str(context.exception)) 
        
        # testing for correct dimension specification
        RM = RateMatrix(data=tst, dim=5)
        numpy.testing.assert_allclose(tst, RM.data,
                                      rtol=1.0e-5, atol=1.0e-12)   
        
        
    def test_rate_value_setting(self):
        """(RATEMATRIX) Testing that a rate can be set
        
        """
        
        # testing setting the rate values with consistent depopulation rate
        RM = RateMatrix(dim=5)
        
        RM.set_rate((2,3), 12.0)
        
        assert RM.data[3,3] == -12.0
        assert RM.data[2,3] == 12.0
        assert RM.data[3,2] == 0.0
        
        RM.set_rate((1,3), 10.0)
        
        assert RM.data[3,3] == -22.0
        assert RM.data[1,3] == 10.0
        assert RM.data[2,3] == 12.0
        assert RM.data[3,1] == 0.0
        
        with self.assertRaises(Exception) as context:
            RM.set_rate((3,3), 10.0)
        self.assertTrue("Diagonal (depopulation) rates cannot be set"
                        in str(context.exception))
        
        
        
        
        