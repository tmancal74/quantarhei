# -*- coding: utf-8 -*-

from .operators import SelfAdjointOperator

import numpy
import scipy


class TransitionDipoleMoment(SelfAdjointOperator):
    
    def __init__(self,dim=None,data=None):
        self.data = data

        if not self.check_selfadjoint():
            raise Exception("The data of this operator have"
            +" to be represented by 3 selfadjoint matrices") 
         
 
    def check_selfadjoint(self):
        a = numpy.allclose(numpy.transpose(numpy.conj(self.data[:,:,0])), 
         self.data[:,:,0])        
        b = numpy.allclose(numpy.transpose(numpy.conj(self.data[:,:,1])), 
         self.data[:,:,1])
        c = numpy.allclose(numpy.transpose(numpy.conj(self.data[:,:,2])), 
         self.data[:,:,2])
        return (a and b) and c   
        
    def transform(self,SS):
        """
        This function transforms the Operator into a different basis, using
        a given transformation matrix.
        """
        S1 = scipy.linalg.inv(SS)
        for i in range(3):
            self.data[:,:,i] = numpy.dot(S1,numpy.dot(self.data[:,:,i],SS))
        
          
    def dipole_strength(self,from_state=0, to_state=1):
        d = numpy.zeros(3,dtype=numpy.float64)        
        for i in range(3):        
            d[i] = self.data[from_state,to_state,i]
        return numpy.dot(d,d)
