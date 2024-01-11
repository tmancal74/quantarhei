# -*- coding: utf-8 -*-
import numpy

from ... import REAL, COMPLEX
from ..propagators.dmevolution import ReducedDensityMatrixEvolution

class OQSStateVectorEvolution:
     
    
    def __init__(self, timeaxis=None, psii=None):
        
        
        if timeaxis is not None:
                
            self.TimeAxis = timeaxis
            
            if psii is not None:
                self.dim = psii.data.shape[0]        
                
                self.set_initial_condition(psii)
            else: 
                self.dim = 0
            
                                    
            
    def set_initial_condition(self, psii):
        """
        
        
        """
        self.dim = psii.data.shape[0] 
        self.data = numpy.zeros((self.TimeAxis.length, \
                    self.dim),dtype=numpy.float64)
        self.data[0,:] = psii.data 
        
        
    def get_ReducedDensityMatrixEvolution(self, decoherence=False):
        """ Returns the corresponding reduced density matrix
        
        
        The present implementation cannot account for dephasing
        
        """
        if decoherence:
            raise Exception("Decoherence not yet implemented.")
            
        rhoi = numpy.zeros((self.dim, self.dim),
                           dtype=REAL)
        rhot = ReducedDensityMatrixEvolution(timeaxis=self.TimeAxis, rhoi=rhoi)
        for ii in range(self.dim):
            for jj in range(self.dim):
                rhot.data[:,ii,jj] = self.data[:,ii]*self.data[:,jj]
                
        return rhot
        