# -*- coding: utf-8 -*-

import numpy

from .redfieldtensor import RedfieldRelaxationTensor

class LindbladForm(RedfieldRelaxationTensor):
    """Lindblad form of relaxation tensor
    
    We use the Redfield tensor class and reimplement its implementation
    
    """
    
    def __init__(self, ham, sbi, initialize=True,
                 as_operators=True, name=""):
                     
        super().__init__(ham, sbi, initialize=initialize,
                         as_operators=as_operators, name=name)    
        
        
    def _implementation(self, ham, sbi):
        """Very simple implementation of Lindblad form
        
        """
        # dimension of the operators 
        Na = ham.dim
        
        # number of operators            
        Nb = sbi.N   
        
        Km = numpy.zeros((Nb, Na, Na), dtype=numpy.float64)
        Ld = numpy.zeros((Nb, Na, Na), dtype=numpy.float64)
        for i in range(Nb):
            Km[i,:,:] = numpy.sqrt(sbi.rates[i])*sbi.KK[i,:,:]
            Ld[i,:,:] = numpy.transpose(Km[i,:,:])
        
        print("%%%%%%%%%%%%%%\n Km=", Km)
        self._post_implementation(Km, Km, Ld)

        print(self._is_initialized, self._data_initialized)