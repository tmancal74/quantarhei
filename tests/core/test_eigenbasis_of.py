# -*- coding: utf-8 -*-

import unittest

"""
*******************************************************************************


    Tests of the quantarhei.eigenbasis_of class


*******************************************************************************
"""        
        
import numpy

from quantarhei import eigenbasis_of
from .test_BasisManaged import BasisManagedObject
    
        

class TestEigenbasisOf(unittest.TestCase):
    
    def setUp(self):

        hval = numpy.array([[0.1, 1.0],
                            [1.0, 0.0]])

        self.H = BasisManagedObject(hval,"H")
        self.B = BasisManagedObject(numpy.diag(numpy.ones(2)),"B")
        self.B.data[1,1] = 0.1
        self.P = BasisManagedObject(hval.copy(),"P")
        self.P.data[0,0] = -0.2

        self.h_copy = self.H.data.copy()
        self.b_copy = self.B.data.copy()
        self.p_copy = self.P.data.copy()  
        
        
    def test_of_single_context_single_var(self):
        """Testing single operator in a single basis management context
        
        
        """
        x,Sh = numpy.linalg.eigh(self.H.data)
        Sh1 = numpy.linalg.inv(Sh)
        h_copy_tr = numpy.dot(Sh1,numpy.dot(self.h_copy,Sh))
        
        with eigenbasis_of(self.H):
            self.assertTrue(numpy.allclose(self.H.data,h_copy_tr))
        
        self.assertTrue(numpy.allclose(self.H.data,self.h_copy))
        
        
        
    def test_of_many_contexts_single_var(self):
        """Testing single operator in nested basis management contexts
        
        
        """
        x,Sh = numpy.linalg.eigh(self.H.data)
        Sh1 = numpy.linalg.inv(Sh)
        h_copy_tr1 = numpy.dot(Sh1,numpy.dot(self.h_copy,Sh))
        b_copy_tr1 = numpy.dot(Sh1,numpy.dot(self.b_copy,Sh))
        
        x,Sb = numpy.linalg.eigh(b_copy_tr1)
        Sb1 = numpy.linalg.inv(Sb)
        h_copy_tr2 = numpy.dot(Sb1,numpy.dot(h_copy_tr1,Sb))   
        
        with eigenbasis_of(self.H):
            # in context 1
            #print("In 1",self.H.manager.get_current_basis())
            self.assertTrue(numpy.allclose(self.H.data,h_copy_tr1))
            with eigenbasis_of(self.B):   
                # in context 2
                #print("In 2",self.B.manager.get_current_basis())
                self.assertTrue(numpy.allclose(self.H.data,h_copy_tr2))
                
            # in context 1
            #print("Out 2 in 1",self.H.manager.get_current_basis())
            self.assertTrue(numpy.allclose(self.H.data,h_copy_tr1))
        
        # outside contexts
        #print("Out 1",self.H.manager.get_current_basis())
        self.assertTrue(numpy.allclose(self.H.data,self.h_copy))


    def test_of_many_contexts_many_vars_1(self):
        """Testing several operators in nested basis management contexts
        
        
        """
        x,Sh = numpy.linalg.eigh(self.H.data)
        Sh1 = numpy.linalg.inv(Sh)
        h_copy_tr1 = numpy.dot(Sh1,numpy.dot(self.h_copy,Sh))
        b_copy_tr1 = numpy.dot(Sh1,numpy.dot(self.b_copy,Sh))
        
        x,Sb = numpy.linalg.eigh(b_copy_tr1)
        Sb1 = numpy.linalg.inv(Sb)
        h_copy_tr2 = numpy.dot(Sb1,numpy.dot(h_copy_tr1,Sb))   
        b_copy_tr2 = numpy.dot(Sb1,numpy.dot(b_copy_tr1,Sb))
        
        with eigenbasis_of(self.H):
            # in context 1
            self.assertTrue(numpy.allclose(self.H.data,h_copy_tr1))
            self.assertTrue(numpy.allclose(self.B.data,b_copy_tr1))
            with eigenbasis_of(self.B):   
                # in context 2
                #print("In 2",self.B.manager.get_current_basis())
                self.assertTrue(numpy.allclose(self.H.data,h_copy_tr2))
                self.assertTrue(numpy.allclose(self.B.data,b_copy_tr2))
            # in context 1
            #print("Out 2 in 1",self.H.manager.get_current_basis())
            self.assertTrue(numpy.allclose(self.H.data,h_copy_tr1))
            self.assertTrue(numpy.allclose(self.B.data,b_copy_tr1))
        
        # outside contexts
        #print("Out 1",self.H.manager.get_current_basis())
        self.assertTrue(numpy.allclose(self.H.data,self.h_copy))
        self.assertTrue(numpy.allclose(self.B.data,self.b_copy))

        
#b_copy = numpy.dot(Sh1,numpy.dot(b_copy,Sh))
#
#with eigenbasis_of(H):
#    print("In eigenstates of H: should be 0.0")
#    print(1)
#    print(H.data-h_copy)
#    print(2)
#    print(B.data-b_copy)
#
#    x,Sb = numpy.linalg.eigh(b_copy)
#    Sb1 = numpy.linalg.inv(Sb)    
#    h_copy = numpy.dot(Sb1,numpy.dot(h_copy,Sb))
#    b_copy = numpy.dot(Sb1,numpy.dot(b_copy,Sb))    
#    with eigenbasis_of(B):
#        print("In eigenstates of B: should be 0.0")
#        print(1)
#        print(H.data-h_copy)
#        print(2)
#        print(B.data-b_copy)
#        print(3)
#        print(P.data)
#    print("Back in eigenstates of H")
#    print(1)
#    print(H.data)
#    print(2)
#    print(B.data)
#    print(3)
#    print(P.data)
#    
#print("After")
#print(1)
#print(H.data)
#print(2)
#print(B.get_current_basis())
#print(B.data)
#print(3)
#print(P.get_current_basis())  # here is a problem
#print(P.data)

