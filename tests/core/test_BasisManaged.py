# -*- coding: utf-8 -*-

import unittest

"""
*******************************************************************************


    Tests of the quantarhei.core.BasisManaged class


*******************************************************************************
"""        
        
import numpy

from quantarhei.core import BasisManaged
from quantarhei.core import BasisManagedReal
    
    
class BasisManagedObject(BasisManaged):
    """Test object

    This class should have methods for conversion of units

    """
    
    data = BasisManagedReal("data",shape=(2,2))

    def __init__(self,dat,name):
        self.data = dat
        self.name = name
        
    def __str__(self):
        return super().__str__()+" : "+self.name
        
        
    def diagonalize(self):
        dd, SS = numpy.linalg.eigh(self._data)
        
        self._data = numpy.diag(dd)
        
        return SS
        
    def transform(self,SS,inv=None):
        if inv is None:
            S1 = numpy.linalg.inv(SS)
        else:
            S1 = inv
        
        self._data = numpy.dot(S1,numpy.dot(self._data,SS))    

class TestBasisManaged(unittest.TestCase):
    
    def setUp(self):
        dat = numpy.zeros((2,2),dtype=numpy.float)
        self.u = BasisManagedObject(dat,"test")
        
    def test_inherited_method(self):
        """Test of the inheritance from the BasisManaged class
        
        """       
        cb = self.u.get_current_basis()
        self.assertEqual(cb,0)
        
        
        

