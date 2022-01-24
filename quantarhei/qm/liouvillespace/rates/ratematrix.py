# -*- coding: utf-8 -*-
import numpy

from ....core.matrixdata import MatrixData

class RateMatrix(MatrixData):
    """Represents a population transfer rate matrix
    
    
    """

    def __init__(self, dim=None, data=None):
        
        self.N = 0
        
        if dim is not None:
            self.N = dim
            
        if data is not None:
            # check if data are rectangular
            if data.shape[0] != data.shape[1]:
                raise Exception("Expecting rectangular matrix")
                
            if self.N == 0:
                self.N = data.shape[0]
                self.data = data
            elif self.N != dim:
                raise Exception("Inconsistent data dimensions")
                   
        else:
            if self.N == 0:
                raise Exception("One of the arguments has to be specified")
                
            else:
                self.data = numpy.zeros((self.N,self.N), dtype=numpy.float64)
            
           
        
    def set_rate(self, pos, value):
        """ Sets a value of a rate between two states
        
        Diagonal depopulation rates are automatically updated
        
        
        """
        N = pos[0]
        M = pos[1]
        orig_val = self.data[N,M]
        self.data[N,M] = value
        
        self.data[M,M] += orig_val
        self.data[M,M] -= value
        
        