# -*- coding: utf-8 -*-

import numpy

from .superoperator import SuperOperator
from .secular import Secular
from ...core.saveable import Saveable

class RelaxationTensor(SuperOperator, Secular, Saveable):
    """Basic class representing a relaxation tensor
    
    
    """
    
    def __init__(self):
        
        self._initialize_basis()
            
        self._data_initialized = False
        self.name = ""
        self.as_operators = False

 
    def _initialize_basis(self):

        # Set the currently used basis
        cb = self.manager.get_current_basis()
        self.set_current_basis(cb)
        # unless it is the basis outside any context
        if cb != 0:
            self.manager.register_with_basis(cb, self)        



    def secularize(self, legacy=True):
        """Secularizes the relaxation tensor


        """
        
        if legacy:
            
            if self.as_operators:
                self.convert_2_tensor()
                #raise Exception("Cannot be secularized in the operator form")
                
            if True:
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
                                            
        else:
            
            super().secularize()


    def _set_population_rates_from_operators(self):
        raise Exception("Method _set_population_rates_from_operators()"+
                        " needs to be implemented in a class inheriting"+
                        " from Secular")


    def _set_population_rates_from_tensor(self):
        """ 
        
        """
        if self.data.ndim == 4:            
            self.secular_KK = numpy.einsum("iijj->ij", self.data)
        else:      
            self.secular_KK = numpy.einsum("hiijj->hij", self.data)            


    def _set_dephasing_rates_from_operators(self):
        raise Exception("Method _set_dephasing_rates_from_operators()"+
                        " needs to be implemented in a class inheriting"+
                        " from Secular")


    def _set_dephasing_rates_from_tensor(self):
        """
        
        """
        if self.data.ndim == 4:
            N = self.data.shape[0]
            
            self.secular_GG = numpy.einsum("ijij->ij", self.data)
            for ii in range(N):
                self.secular_GG[ii,ii] = 0.0
        else:
            N = self.data.shape[1]
            self.secular_GG = numpy.einsum("hijij->hij", self.data)            
            for ii in range(N):
                self.secular_GG[:,ii,ii] = 0.0

                               
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


    def convert_2_tensor(self):
        """Converts from operator to tensor form
        
        """
        pass


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
            
    
    def __mult__(self, scalar):
        """Multiplication of the Tensor by a scalar
        
        """
        import numbers
        
        if not isinstance(scalar, numbers.Number):
            raise Exception("Only multiplication by numbers is implemented")
            
        if self.as_operators:
            raise Exception("Multiplication in operator form not implemented")
            
        self._data = self._data*scalar
        return self

    
    def __rmult__(self, scalar):
        return self.__mult__(scalar)
    
    
    def __add__(self, other):
        self._data += other._data
        return self


    def __iadd__(self, other):
        return self.__add__(other)
    
    
    def _rhs(self, rho):
        """Applies the tensor to a given matrix
        
        """
        
        if self.as_operators:
            
            return self._rhs_as_operators(rho)
        
        else:
            
            return self._rhs_as_tensor(rho)
        
        
    def _rhs_as_operators(self, rho):
        """Applies the tensor in form of a set of operators to a given matrix
        
        Parameters
        ----------
        
        rho : complex or real array
            Array representing density matrix
        
        
        Returns
        -------
        
        Complex array
        
        
        """
        pass
    

    def _rhs_as_tensor(self, rho):
        """Applies the tensor to a given matrix
        
        Parameters
        ----------
        
        rho : complex or real array
            Array representing density matrix
        
        
        Returns
        -------
        
        Complex array
        
        
        """
        return numpy.tensordot(self.data, rho)
        
        