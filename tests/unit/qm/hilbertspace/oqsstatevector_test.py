# -*- coding: utf-8 -*-
import unittest


"""
*******************************************************************************


    Tests of the quantarhei.qm.hilbertspace.oqsstatevector package


*******************************************************************************
"""

from quantarhei import OQSStateVector
from quantarhei import REAL
import quantarhei as qr

class TestOQSStateVector(unittest.TestCase):
    """Tests for the statevector package
    
    
    """
    
    def test_of_state_vector_creation(self):
        """Testing StateVector creation """
        
        psi = OQSStateVector(3)
        
        self.assertEqual(psi.dim,3)
        
        psi = OQSStateVector()
        
        self.assertFalse(psi._initialized)
        


    def test_of_sv_creation_from_data(self):
        """Testing StateVector creation with data"""
        
        import numpy
        
        vec = numpy.zeros((1,3), dtype=REAL)
        
        with self.assertRaises(Exception):
            psi = OQSStateVector(data=vec)
            
        vec = numpy.zeros(3, dtype=REAL)
        psi = OQSStateVector(data=vec)
        
        self.assertEqual(psi.dim,3)
        

    def test_creation_from_list(self):
        """Tests creation from non numpy array """
        import numpy
        
        vec = [0.1, 2.0, 0.0]
        
        psi = OQSStateVector(data = vec)
        
        self.assertAlmostEqual(psi.norm(), numpy.sqrt(0.1**2 + 4.0))
        self.assertAlmostEqual(psi.puredot(psi), psi.norm()**2)
        
        
        
    def test_evolution(self):
        
        import numpy
        
        KK = numpy.array([[-1.0/200.0, 1.0/200.], 
                          [1.0/200.0, -1.0/200.0]],
                         dtype=REAL)
        
        timea = qr.TimeAxis(0.0, 1000, 1.0)
        
        #prop = qr.OQSStateVectorPropagator(timea, KK)
        
        #psii = qr.OQSStateVector(data=[numpy.sqrt(0.1),
        #                               numpy.sqrt(0.9)])
        
        
        
        #psit = prop.propagate(psii)
        
        
      
        

        