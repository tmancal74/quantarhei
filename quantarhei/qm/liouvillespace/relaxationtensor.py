# -*- coding: utf-8 -*-

import numpy

from ...core.managers import BasisManaged
from ...utils.types import BasisManagedComplexArray

class RelaxationTensor(BasisManaged):

    data = BasisManagedComplexArray("data")
    
    def __init__(self):
        
        # Set the currently used basis
        cb = self.manager.get_current_basis()
        self.set_current_basis(cb)
        # unless it is the basis outside any context
        if cb != 0:
            self.manager.register_with_basis(cb,self)
            
        self._data_initialized = False
        self.name = ""


    def secularize(self):
        """Secularizes the relaxation tensor


        """
        if self.as_operators:
            raise Exception("Cannot be secularized in an operator form")
            
        else:
            if self.data.ndim == 4:
                N = self.data.shape[0]
                for ii in range(N):
                    for jj in range(N):
                        for kk in range(N):
                            for ll in range(N):
                                if not (((ii == jj) and (kk == ll)) 
                                    or ((ii == kk) and (jj == ll))) :
                                        self.data[ii,jj,kk,ll] = 0
            else:  
                N = self.data.shape[1]
                for ii in range(N):
                    for jj in range(N):
                        for kk in range(N):
                            for ll in range(N):
                                if not (((ii == jj) and (kk == ll)) 
                                    or ((ii == kk) and (jj == ll))) :
                                        self.data[:,ii,jj,kk,ll] = 0                                        

                               
    def transform(self, SS, inv=None):
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
        

        if not self._data_initialized:
            print("NOT INITIALIZED")
            return
        
        if (self.manager.warn_about_basis_change):
                print("\nQr >>> Relaxation tensor '%s' changes basis"
                      %self.name)
           
        if inv is None:
            S1 = numpy.linalg.inv(SS)
        else:
            S1 = inv
        dim = SS.shape[0]
        
        if self._data.ndim == 4:
            for c in range(dim):
                for d in range(dim):
                    self._data[:,:,c,d] = \
                    numpy.dot(S1,numpy.dot(self._data[:,:,c,d],SS))
                    
            for a in range(dim):
                for b in range(dim):
                    self._data[a,b,:,:] = \
                    numpy.dot(S1,numpy.dot(self._data[a,b,:,:],SS))
        else:

            for tt in range(self._data.shape[0]):
                for c in range(dim):
                    for d in range(dim):
                        self._data[tt,:,:,c,d] = \
                            numpy.dot(S1,numpy.dot(self._data[tt,:,:,c,d],SS))
                    
                for a in range(dim):
                    for b in range(dim):
                        self._data[tt,a,b,:,:] = \
                            numpy.dot(S1,numpy.dot(self._data[tt,a,b,:,:],SS))            


    def updateStructure(self):
        """ Recalculates dephasing and depopulation rates
        
        """
        
        if self._data.ndim == 4:
            # depopulation rates 
            for nn in range(self.dim):
                #for ii in range(0,self.data.shape[1]):
                #    if ii != nn:
                #        self.data[nn,nn,nn,nn] -= self.data[ii,ii,nn,nn]
                self._data[nn,nn,nn,nn] -= (numpy.trace(self._data[:,:,nn,nn])
                                            - self._data[nn,nn,nn,nn])
                
            # dephasing rates 
            for nn in range(self.dim):    
                for mm in range(nn+1,self.dim):
                    self._data[nn,mm,nn,mm] = -(self._data[nn,nn,nn,nn]
                                              +self._data[mm,mm,mm,mm])/2.0
                    self._data[mm,nn,mm,nn] = self._data[nn,mm,nn,mm] 

        else:
            # depopulation rates 
            for nn in range(self.dim):
                #for ii in range(0,self.data.shape[1]):
                #    if ii != nn:
                #        self.data[nn,nn,nn,nn] -= self.data[ii,ii,nn,nn]
                self._data[:,nn,nn,nn,nn] -= (numpy.trace(self._data[:,:,:,nn,nn],
                                                          axis1=1,axis2=2)
                                            - self._data[:,nn,nn,nn,nn])
                
            # dephasing rates 
            for nn in range(self.dim):    
                for mm in range(nn+1,self.dim):
                    self._data[:,nn,mm,nn,mm] = -(self._data[:,nn,nn,nn,nn]
                                              +self._data[:,mm,mm,mm,mm])/2.0
                    self._data[:,mm,nn,mm,nn] = self._data[:,nn,mm,nn,mm] 
            