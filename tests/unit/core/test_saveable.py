# -*- coding: utf-8 -*-

import unittest
import h5py
import numpy

"""
*******************************************************************************


    Tests of the quantarhei.Aggregate class


*******************************************************************************
"""
from quantarhei import Saveable


class TSaveable(Saveable):
    
    def __init__(self):
        
        self.a = 1.0
        self.text = None
        self._txt = None
        self.b1 = False
        self.b2 = False
        self.b3 = False
        
        self.dat = None
        
    def set_a(self, a):
        self.a = a
        
    def set_str(self, text1, text2):
        self.text = text1
        self._txt = text2
        
    def set_bool(self, b1, b2, b3):
        self.b1 = b1
        self.b2 = b2
        self.b3 = b3
        
    def set_data(self, dat):
        self.dat = numpy.array(dat)
        
    def set_saveable(self, obj):
        self.obj = obj
        

class TestSaveable(unittest.TestCase):
    """Tests for the Saveable class
    
    
    """
    
    def setUp(self):
        
        self.obj1 = TSaveable()
        
        self.obj1.set_a(10.0)
        self.obj1.set_str("Hi, there", "Hello World")
        self.obj1.set_bool(True, False, True)
        
        dat = numpy.zeros((5,5), dtype=numpy.complex128)
        
        dat[3,4] = 1.3455 + 1j*0.3456
        
        self.obj1.set_data(dat)
        
        obj2 = TSaveable()
        obj2.set_str("I am from second tier", "and me, too")
        
        self.obj1.set_saveable(obj2)
        
    def test_saving_and_loading_1(self):
        """Testing saving an loading Saveable objects
        
        
        """
        
        obj1 = self.obj1
        
        with h5py.File("test_file_1",driver="core", 
                           backing_store=False) as f:
            
            obj1.save(f)
            
            
            obj2 = TSaveable()
            obj2.load(f)
            
            
        self.assertEqual(obj1.a,obj2.a)
        self.assertEqual(obj1.text,obj2.text)
        self.assertEqual(obj1._txt,obj2._txt)
        self.assertEqual(obj1.b1,obj2.b1)
        self.assertEqual(obj1.b2,obj2.b2)
        self.assertEqual(obj1.b3,obj2.b3)
        
        numpy.testing.assert_array_equal(obj1.dat,obj2.dat)
        
        #self.assertEqual(obj1.obj.text,obj2.obj.text)
            
            
            
        
