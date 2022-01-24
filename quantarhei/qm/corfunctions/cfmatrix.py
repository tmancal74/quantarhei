# -*- coding: utf-8 -*-
"""
    Quantarhei package (http://www.github.com/quantarhei)

    cfmatrix module

"""
import numpy

from ...core.saveable import Saveable
from ...core.time import TimeAxis

from .correlationfunctions import c2h
from .correlationfunctions import c2g
from ... import REAL, COMPLEX


class CorrelationFunctionMatrix(Saveable):
    """Matrix of correlation functions specifying cross-correlations

    Parameters
    ----------
    timeaxis : TimeAxis
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

    def __init__(self, timeaxis=TimeAxis(0.0,1,1.0), nob=0, nof=0):
        # Number of baths
        self.nob = nob

        # Number of functions in the set
        self.nof = nof

        # Pointer pointing from position in the matrix
        # to the list of functions
        self.cpointer = numpy.zeros((nob,nob),dtype=numpy.int32)
        self.cpointer[:,:] = 0

        self._A2 = None
        self._A4 = None

        # empty list for functions
        self.cfuncs = None

        #FIXME: reorganization energies should be defined through _A2 and _A4
        # reorganization energies
        self.lambdas = None
        self.where = None


        # Actual storage of functions
        # here we store correlation functions
        self._cofts = None

        # here we store their first integrals
        self._hofts = None
        # here we store their second integrals
        self._gofts = None

        self.data = None

        # TimeAxis on which the functions are defined
        self.timeAxis = timeaxis

        #FIXME: introduce cut-off time to reduce memory usage
        self._has_cutoff_time = False
        self.max_cutoff_index = 0

        #FIXME: implement transformation of the matrix
        self._is_transformed = False

        self._initiate_storage()


    def _initiate_storage(self):
        """Initializes storage of functions when object is created or updated

        """
        nob = self.nob
        nof = self.nof
        
        self._A2 = numpy.zeros((nob, nob, nof+1),
                               dtype=REAL)
        if self._is_transformed:
            self._A4 = numpy.zeros((nob, nob, nob, nob, nof),
                                   dtype=REAL)
        else:
            self._A4 = None

        # empty list for functions
        self.cfuncs = [None]*(nof+1)

        #FIXME: reorganization energies should be defined through _A2 and _A4
        # reorganization energies
        self.lambdas = numpy.zeros(nof+1, dtype=REAL)
        self.where = [[]]*(nof+1)

        # Actual storage of functions
        # here we store correlation functions
        self._cofts = numpy.zeros((nof+1, self.timeAxis.length),
                                  dtype=COMPLEX)

        self.data = self._cofts


    def _update_nof_storage(self):
        """Updates the storage of functions in the matrix
        
        When a new function is added, number of functions (self.nof) is
        updated and storage is realocated.
        
        """
        nof = self.nof

        # save all values
        save_A2 = self._A2
        if self._is_transformed:
            save_A4 = self._A4
        save_cfunc = self.cfuncs
        save_lambdas = self.lambdas
        save_where = self.where
        save_cofts = self._cofts

        # Add one place in the matrix
        self.nof += 1

        # reinitiate
        self._initiate_storage()

        # refill
        self._A2[:,:,0:nof+1] = save_A2
        if self._is_transformed:
            self._A4[:,:,:,:,0:nof] = save_A4
        for i in range(nof+1):
            self.cfuncs[i] = save_cfunc[i]
            self.lambdas[i] = save_lambdas[i]
            self.where[i] = save_where[i]
        self._cofts[0:nof+1,:] = save_cofts

    
    def _check_temperature_consistency(self):
        """Checks if the temperature of all correlation functions is the same
        
        Checks that all correlation functions have temperature specified and 
        that they are the same.
        
        
        """
        temp = -1.0
        none_count = 0
        for cf in self.cfuncs:
            if cf is not None:
                T = cf.get_temperature()
                if temp < 0.0:
                    temp = T
                if T != temp:
                    raise Exception("Temperature of "+
                                "CorrelationFunctionMatrix is not consistent")
            else:
                none_count += 1

        if none_count > 1:
            raise Exception()
            
        return temp
    

    def get_temperature(self):
        """Returns temperature of the system
        
        
        Temperature is speficied in the correlation functions. It is returned
        if it is the same for all correlation functions in the matrix. 
        Otherwise, exception is thrown.
        
        
        """
        temp = self._check_temperature_consistency()
        return temp
    

    def get_correlation_function(self,n,m):
        """Returns correlation function associated with sites n and m

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

    def get_correlation_time(self,n,m):
        return self.get_correlation_function(n,m).params[0]["cortime"]

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
        """Returns correlation function corresponding to two states n and m

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
            ret = numpy.zeros(self.timeAxis.length, dtype=COMPLEX)
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
                print(Nt)
                ret = numpy.zeros((nob,nob,nob,nob,Nt),dtype=COMPLEX)
                for k in range(self.nof):
                    for it in range(Nt):
                        ret[:,:,:,:,it] += self._A4[:,:,:,:,k]\
                        *self._cofts[k+1,it]
                return ret
            else:
                ret = numpy.zeros((nob,nob,nob,nob),dtype=COMPLEX)
                it = self.timeAxis.nearest(t)
                for k in range(self.nof):
                    ret += self._A4[:,:,:,:,k]*self._cofts[k+1,it]
                return ret
        else:
            raise Exception()


    def get_hoft4(self,a,b,c,d):
        if self._is_transformed:
            ret = numpy.zeros(self.timeAxis.length, dtype=COMPLEX)
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
                ret = numpy.zeros((nob,nob,nob,nob,Nt),dtype=COMPLEX)
                for k in range(self.nof):
                    for it in range(Nt):
                        ret[:,:,:,:,it] += self._A4[:,:,:,:,k]\
                        *self._hofts[k+1,it]
                return ret
            else:
                ret = numpy.zeros((nob,nob,nob,nob),dtype=COMPLEX)
                it = self.timeAxis.nearest(t)
                for k in range(self.nof):
                    ret += self._A4[:,:,:,:,k]*self._hofts[k+1,it]
                return ret
        else:
            raise Exception()


    def get_goft4(self,a,b,c,d):
        if self._is_transformed:
            ret = numpy.zeros(self.timeAxis.length, dtype=COMPLEX)
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
                ret = numpy.zeros((nob,nob,nob,nob,Nt),dtype=COMPLEX)
                for k in range(self.nof):
                    for it in range(Nt):
                        ret[:,:,:,:,it] += self._A4[:,:,:,:,k]\
                        *self._gofts[k+1,it]
                return ret
            else:
                ret = numpy.zeros((nob,nob,nob,nob),dtype=COMPLEX)
                it = self.timeAxis.nearest(t)
                for k in range(self.nof):
                    ret += self._A4[:,:,:,:,k]*self._gofts[k+1,it]
                return ret
        else:
            raise Exception()


    def set_correlation_function(self, fce, where, iof=None):
        """Sets a correlation function to the matrix

        Parameters
        ----------
        iof : int
            Index of the function in the list. Must be larger than 0
        fce : CorrelationFunction
            Correlation function represented by the CorrelationFunction object
        where : list of tuples
            List of positions in the matrix to which the function should
            be essigned.

        """

        if fce in self.cfuncs:

            # if the function is already in, iof parameter is ignored
            i = self.cfuncs.index(fce)
            iof = i
            self._update_where(iof, fce, where)

        else:

            # if iof is not specified explicitely, it is chosen as the next
            # available index
            if iof is None:
                iof = self.nof + 1
            while iof > self.nof:
                self._update_nof_storage()

            if iof == 0:
                raise Exception("Zeros correlation function is \
                                reserved for an empty function.")
            if iof >= self.nof+1:
                raise Exception("Index out of bounds")

            try:
                ic = fce.axis.nearest(fce.cutoff_time)
            except:
                ic = fce.axis.data[fce.axis.length-1]

            ic = int(ic)
            if ic > self.max_cutoff_index:
                self.max_cutoff_index = ic

            # this sets also self._confs; this can undergo transformations
            self.data[iof,:] = fce.data

            self.lambdas[iof] = fce.lamb
            self.cfuncs[iof]  = fce

            self._update_where(iof, fce, where)

        return iof


    def get_index_by_where(self, where):

        ret = -1
        for i in range(self.nof+1):
            if where in self.where[i]:
                ret = i

        return ret


    def _update_where(self, iof, fce, where):

        for loc in where:
            if loc in self.where[iof]:
                raise Exception("Location in correlation matrix already taken")
            else:
                self.where[iof].append(loc)

        for wr in where:
            self._A2[wr[0],wr[1],iof] = fce.lamb
            #print(wr[0], wr[1], iof, fce.lamb)

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
        self._A4 = numpy.zeros((nob,nob,nob,nob,nof), dtype=REAL)

        for a in range(nob):
            for b in range(nob):
                for c in range(nob):
                    for d in range(nob):
                        for k in range(nof):

                           for n in range(nob):
                                for m in range(nob):
                                    self._A4[a,b,c,d,k] += SS[n,a]*SS[n,b]*\
                                    self._A2[n,m,k+1]*SS[m,c]*SS[m,d]

        self._is_transformed = True
        #print(self._A4[1,1,1,1,0])
        #print(self._A2[1,1,1])
        #print(SS[:,1])
