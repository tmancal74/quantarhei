# -*- coding: utf-8 -*-

import numpy

from ...core.managers import BasisManaged
from ...utils.types import BasisManagedComplexArray

class RelaxationTensor(BasisManaged):

    data = BasisManagedComplexArray("data")
    
    
    def secularize(self):
        """Secularizes the relaxation tensor


        """
        if self.as_operators:
            raise Exception("Cannot be secularized in an opeator form")
            
        else:
            N = self.data.shape[0]
            for ii in range(N):
                for jj in range(N):
                    for kk in range(N):
                        for ll in range(N):
                            if not (((ii == jj) and (kk == ll)) 
                                or ((ii == kk) and (jj == ll))) :
                                    self.data[ii,jj,kk,ll] = 0
                                        
                                        
    def transform(self,SS,inv=None):
        """Transformation of the tensor by a given matrix
        
        
        This function transforms the Operator into a different basis, using
        a given transformation matrix.
        
        Parameters
        ----------
         
        SS : matrix, numpy.ndarray
            transformation matrix
            
        inv : matrix, numpy.ndarray
            inverse of the transformation matrix
            
        """        
        if inv is None:
            S1 = numpy.linalg.inv(SS)
        else:
            S1 = inv
        dim = SS.shape[0]
        
        
        for c in range(dim):
            for d in range(dim):
                self._data[:,:,c,d] = \
                numpy.dot(S1,numpy.dot(self._data[:,:,c,d],SS))
                
        for a in range(dim):
            for b in range(dim):
                self._data[a,b,:,:] = \
                numpy.dot(S1,numpy.dot(self._data[a,b,:,:],SS))
        
#        RR = numpy.zeros((dim,dim,dim,dim), dtype=numpy.complex128)
#        for ag in range(dim):
#            for bg in range(dim):
#                for cg in range(dim):
#                    for dg in range(dim):
#                        rr = 0.0
#                        for a in range(dim):
#                            for b in range(dim):
#                                for c in range(dim):
#                                    for d in range(dim):
#                                        rr += S1[ag,a]*SS[b,bg]*self._data[a,b,c,d]*S1[c,cg]*SS[dg,d]
#                        RR[ag,bg,cg,dg] = rr
#        self._data = RR
#        