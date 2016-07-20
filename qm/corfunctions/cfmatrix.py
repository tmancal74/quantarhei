# -*- coding: utf-8 -*-
import numpy


from .correlationfunctions import c2h
from .correlationfunctions import c2g


class CorrelationFunctionMatrix: #(MatrixData,TimeDependent):
    """Matrix of correlation functions specifying cross-correlations 
    
    Parameters
    ----------
    timeaxis : cu.oqs.time.TimeAxis
        TimeAxis on which the correlation function is specified
    nob : int
        Number of individual diffent bath, i.e. the dimension
        of the correlation matrix
    nof : int
        Number of specific functions by which the matrix is defined
        
    Notes
    -----
    The correlation matrix specifies energy gap correlation functions for
    different molecules and their transitions. Some of the functions might
    have the same form and so they are stored only once. However, two
    molecules described by the same function can be completely uncorrelated,
    or have a cross-correlation specified by a completely different function
    (which then should be on the list of functions).
    
    """
    
    def __init__(self,timeaxis,nob,nof):        
        # Number of baths 
        self.nob = nob

        # Number of functions in the set 
        self.nof = nof 
        
        # Pointer pointing from position in the matrix
        # to the list of functions
        self.cpointer = numpy.zeros((nob,nob),dtype=numpy.int32)
        self.cpointer[:,:] = 0 
        
        self._A2 = numpy.zeros((nob,nob,nof+1),dtype=numpy.float64)
        self._A4 = None 
        
        # empty list for functions
        self.cfuncs = [None]*(nof+1)
        
        #FIXME: reorganization energies should be defined through _A2 and _A4 
        # reorganization energies
        self.lambdas = numpy.zeros(nof+1,dtype=numpy.float64)
        self.where = [None]*(nof+1)
        
        
        # Actual storage of functions
        # here we store correlation functions
        self._cofts = numpy.zeros((nof+1,timeaxis.length),
                                  dtype=numpy.complex64)
        # here we store their first integrals
        self._hofts = None
        # here we store their second integrals
        self._gofts = None
        
        self.data = self._cofts
        
        # TimeAxis on which the functions are defined
        self.timeAxis = timeaxis
        
        #FIXME: introduce cut-off time to reduce memory usage
        self._has_cutoff_time = False
        self.max_cutoff_index = 0
        
        #FIXME: implement transformation of the matrix
        self._is_transformed = False
        
    def get_correlation_function(self,n,m):
        """Returns correlation function from the matrix 
        
        Parameters
        ----------
        n, m : int
            indices of the matrix        
        """
        #return self.data[self.cpointer[n,m],:]
        return self.cfuncs[self.cpointer[n,m]]        
                 
    def get_reorganization_energy(self,n,m):
        """Returns reorganization energie in the site basis """
        #FIXME: should this should use _A2  ????
        return self.lambdas[self.cpointer[n,m]]

    def get_reorganization_energy4(self,a,b,c,d):
        """Returns reorganization energy in the transformed basis
        
        Parameters
        ----------
        a,b,c,d : int
            indices of states        
        
        """
        if self._is_transformed:
            lm = 0.0
            for k in range(self.nof):
                lm += self._A4[a,b,c,d,k]
            return lm
        else:
            raise Exception()

    def get_coft(self,n,m):
        """Returns correlation function from the matrix 
        
        Parameters
        ----------
        n, m : int
            indices of the matrix        
        """
        return self._cofts[self.cpointer[n,m],:]

    def get_hoft(self,n,m):
        """Returns first integral of the correlation function from the matrix 
        
        Parameters
        ----------
        n, m : int
            indices of the matrix        
        """
        return self._hofts[self.cpointer[n,m],:]

    def get_goft(self,n,m):
        """Returns the lineshape function from the matrix 
        
        Parameters
        ----------
        n, m : int
            indices of the matrix        
        """
        return self._gofts[self.cpointer[n,m],:]

    def get_coft4(self,a,b,c,d):
        if self._is_transformed:
            ret = numpy.zeros(self.timeAxis.length)
            for k in range(self.nof):
                ret += self._A4[a,b,c,d,k]*self._cofts[k+1,:]
            return ret
        else:
            raise Exception()
        
    def get_coft_matrix(self,t=None):   
        """Returns full matrix of correlation functions 

        It is expected that the matrix is transformed into new basis after
        it was specified in the site basis.

        Parameters
        ----------
        t : float, optional
            If t is specified, the routine returns the matrix
            at the particular time, otherwise a full matrix is returned.
                
        
        """         
        nob = self.nob
        if self._is_transformed:
            if t is None:
                Nt = self.max_cutoff_index
                ret = numpy.zeros((nob,nob,nob,nob,Nt),dtype=numpy.complex128)
                for k in range(self.nof):
                    for it in range(Nt):
                        ret[:,:,:,:,it] += self._A4[:,:,:,:,k]\
                        *self._cofts[k+1,it]
                return ret
            else:
                ret = numpy.zeros((nob,nob,nob,nob),dtype=numpy.complex128)
                it = self.timeAxis.nearest(t)
                for k in range(self.nof):
                    ret += self._A4[:,:,:,:,k]*self._cofts[k+1,it]
                return ret
        else:
            raise Exception()
        
    def get_hoft4(self,a,b,c,d):
        if self._is_transformed:
            ret = numpy.zeros(self.timeAxis.length)
            for k in range(self.nof):
                ret += self._A4[a,b,c,d,k]*self._hofts[k+1,:]
            return ret
        else:
            raise Exception()
 
    def get_hoft_matrix(self,t=None):
        """Returns full matrix of once integrated correlation functions 

        It is expected that the matrix is transformed into new basis after
        it was specified in the site basis.

        Parameters
        ----------
        t : float, optional
            If t is specified, the routine returns the matrix
            at the particular time, otherwise a full matrix is returned.
                
        
        """         
        nob = self.nob
        if self._is_transformed:
            if t is None:
                Nt = self.max_cutoff_index
                ret = numpy.zeros((nob,nob,nob,nob,Nt),dtype=numpy.complex128)
                for k in range(self.nof):
                    for it in range(Nt):
                        ret[:,:,:,:,it] += self._A4[:,:,:,:,k]\
                        *self._hofts[k+1,it]
                return ret
            else:
                ret = numpy.zeros((nob,nob,nob,nob),dtype=numpy.complex128)
                it = self.timeAxis.nearest(t)
                for k in range(self.nof):
                    ret += self._A4[:,:,:,:,k]*self._hofts[k+1,it]
                return ret
        else:
            raise Exception()
       
    def get_goft4(self,a,b,c,d):
        if self._is_transformed:
            ret = numpy.zeros(self.timeAxis.length)
            for k in range(self.nof):
                ret += self._A4[a,b,c,d,k]*self._gofts[k+1,:]
            return ret
        else:
            raise Exception()
        
    def get_goft_matrix(self,t=None):
        """Returns full matrix of lineshape functions 

        It is expected that the matrix is transformed into new basis after
        it was specified in the site basis.

        Parameters
        ----------
        t : float, optional
            If t is specified, the routine returns the matrix
            at the particular time, otherwise a full matrix is returned.
                
        
        """         
        nob = self.nob
        if self._is_transformed:
            if t is None:
                Nt = self.max_cutoff_index
                ret = numpy.zeros((nob,nob,nob,nob,Nt),dtype=numpy.complex128)
                for k in range(self.nof):
                    for it in range(Nt):
                        ret[:,:,:,:,it] += self._A4[:,:,:,:,k]\
                        *self._gofts[k+1,it]
                return ret
            else:
                ret = numpy.zeros((nob,nob,nob,nob),dtype=numpy.complex128)
                it = self.timeAxis.nearest(t)
                for k in range(self.nof):
                    ret += self._A4[:,:,:,:,k]*self._gofts[k+1,it]
                return ret
        else:
            raise Exception()
        

    def set_correlation_function(self,iof,fce,where):
        """Sets a correlation function to the matrix 
        
        Parameters
        ----------
        iof : int
            Index of the function in the list. Must be larger than 0
        fce : array of real values
            Correlation function represented by an array of its values.
            The length of the array must be the same as the lenth
            as the TimeAxis.
        where : list of tuples
            List of positions in the matrix to which the function should
            be essigned.
            
        """    
        if iof == 0:
            raise Exception("Zeros correlation function is \
            reserved for an emty function.")
        if iof >= self.nof+1:
            raise Exception("Index out of bounds")

        try:
            ic = fce.timeAxis.nearest(fce.cutoff_time)
        except:
            ic = fce.timeAxis.time[fce.timeAxis.length-1]
            
        if (ic > self.max_cutoff_index):
            self.max_cutoff_index = ic
            
        self.data[iof,:] = fce.data
        
        self.lambdas[iof] = fce.lamb
        self.cfuncs[iof]  = fce
        
        for wr in where:
            self._A2[wr[0],wr[1],iof] = fce.lamb
        self.where[iof] = where
        
        for loc in where: 
            ii = loc[0]
            jj = loc[1]
            self.cpointer[ii,jj] = iof

    def create_one_integral(self):
        self._hofts = numpy.zeros((self.nof+1,self.timeAxis.length),
                                  dtype=numpy.complex128)
        for ii in range(self.nof+1):
            self._hofts[ii,:] = c2h(self.timeAxis,self._cofts[ii,:])
    
    def create_double_integral(self):
        self._gofts = numpy.zeros((self.nof+1,self.timeAxis.length),
                                  dtype=numpy.complex128)
        for ii in range(self.nof+1):
            self._gofts[ii,:] = c2g(self.timeAxis,self._cofts[ii,:])
            
            
    def transform(self,SS):
        nob = self.nob
        nof = self.nof
        self._A4 = numpy.zeros((nob,nob,nob,nob,nof))

        for a in range(nob):
            for b in range(nob):
                for c in range(nob):
                    for d in range(nob):
                        for k in range(nof):
         
                           for n in range(nob):
                                for m in range(nob):
                                    self._A4[a,b,c,d,k] += SS[n,a]*SS[n,b]*\
                                    self._A2[n,m,k]*SS[m,c]*SS[m,d]
                                    
        self._is_transformed = True                            
            
