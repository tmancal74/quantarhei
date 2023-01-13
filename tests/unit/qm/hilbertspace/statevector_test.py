# -*- coding: utf-8 -*-


import unittest



"""
*******************************************************************************


    Tests of the quantarhei.qm.hilbertspace.statevector package


*******************************************************************************
"""

from quantarhei import StateVector

class TestStateVector(unittest.TestCase):
    """Tests for the statevector package
    
    
    """
    
    def test_of_state_vector_creation(self):
        """Testing StateVector creation """
        
        psi = StateVector(3)
        
        self.assertEqual(psi.dim,3)
        
        psi = StateVector()
        
        self.assertFalse(psi._initialized)
        


    def test_of_sv_creation_from_data(self):
        """Testing StateVector creation with data"""
        
        import numpy
        
        vec = numpy.zeros((1,3), dtype=numpy.float)
        
        with self.assertRaises(Exception):
            psi = StateVector(data=vec)
            
        vec = numpy.zeros(3, dtype=numpy.float)
        psi = StateVector(data=vec)
        
        self.assertEqual(psi.dim,3)
        
    def test_creation_from_list(self):
        """Tests creation from non numpy array """
        import numpy
        
        vec = [0.1, 2.0, 0.0]
        
        psi = StateVector(data = vec)
        
        self.assertAlmostEqual(psi.norm(), numpy.sqrt(0.1**2 + 4.0))
        self.assertAlmostEqual(psi.dot(psi), psi.norm()**2)
        
        
    def test_of_scalar_product(self):
        """Test StateVector scalar product """
        
        import numpy
        vec1 = numpy.zeros(3, dtype=numpy.float)
        vec2 = numpy.zeros(3, dtype=numpy.float)
        
        vec1[1] = 1.0
        vec2[0] = 1.0
        
        psi1 = StateVector(data=vec1)
        psi2 = StateVector(data=vec2)
        psi3 = StateVector(3)
        
        psi3.data[1] = 0.5
        
        scl12 = psi1.dot(psi2)
        
        self.assertEqual(scl12, 0.0)
        
        scl13 = psi3.dot(psi1)
        
        self.assertEqual(scl13, 0.5)
        
        self.assertEqual(psi1.dot(psi1), 1.0)
        
        
    def test_norm(self):
        """Test StateVector norm 
        """
        
        import numpy
        vec = numpy.zeros(5, dtype=numpy.float)
        vec[1] = 2.0
        vec[4] = 3.0
        
        psi = StateVector(data=vec)
        
        self.assertAlmostEqual(psi.norm(), numpy.sqrt(13.0))
        
        

        