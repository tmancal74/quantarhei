# -*- coding: utf-8 -*-
"""
    Quantarhei package (http://www.github.com/quantarhei)

    cfmatrix module

"""
import numpy

from ...core.saveable import Saveable
from ...core.time import TimeAxis
from ...core.managers import Manager

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

        # Number of system states (not all states are coupled
        #                          to the bath explicitely)
        
        # Pointer pointing from position in the matrix
        # to the list of functions
        self.cpointer = numpy.zeros((nob,nob),dtype=numpy.int32)
        self.cpointer[:,:] = 0

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
        self._has_hofts = False
        
        # here we store their second integrals
        self._gofts = None
        self._has_gofts = False

        self.data = None

        # TimeAxis on which the functions are defined
        self.timeAxis = timeaxis

        #FIXME: introduce cut-off time to reduce memory usage
        self._has_cutoff_time = False
        self.max_cutoff_index = 0


        # 
        #  THIS WILL BE TRANSPLATED TO OPEN SYSTEMS
        #

        #FIXME: implement transformation of the matrix
        self._is_transformed = False
        self._A2 = None
        self._A4 = None
        
        
        
        self._initiate_storage()


    def _initiate_storage(self):
        """Initializes storage of functions when object is created or updated

        """
        nob = self.nob
        nof = self.nof

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

        #
        # THIS WILL BE TRANSPLANTED TO OPEN SYSTEMS
        #
        
        nos = nob + 1 
        
        self._A2 = numpy.zeros((nos, nos, nof+1),
                               dtype=REAL)
        if self._is_transformed:
            self._A4 = numpy.zeros((nos, nos, nos, nos, nof),
                                   dtype=REAL)
        else:
            self._A4 = None


    def _update_nof_storage(self):
        """Updates the storage of functions in the matrix
        
        When a new function is added, number of functions (self.nof) is
        updated and storage is realocated.
        
        """
        nof = self.nof

        # save all values
        
        #
        #  THIS WILL bE TRANSPLANTED TO OPEN SYSTEMS
        #
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
        
        #
        #  THIS WILL BE TRANSPLANTED TO OPEN SYSTEM
        #
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
        """Returns reorganization energie in the site basis 
        """
        #FIXME: should this use _A2  ????
        man = Manager()
        return man.convert_energy_2_current_u(
                            self.lambdas[self.cpointer[n,m]])


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
            raise Exception("Correlation function matrix is not initialized.")


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
        nos = self.nob + 1
        if self._is_transformed:
            if t is None:
                
                Nt = self.max_cutoff_index + 1
                
                # FIXME: some consistent treatment of the cutoff time is needed  
                Nt = self._cofts.shape[1]
                
                ret = numpy.zeros((nos,nos,nos,nos,Nt),dtype=COMPLEX)
                #for k in range(self.nof):
                for k in range(self.nob):
                    for it in range(Nt):
                        #ret[:,:,:,:,it] += self._A4[:,:,:,:,k]\
                        #*self._cofts[k+1,it]
                        ret[:,:,:,:,it] += self._C4[:,:,:,:,k]\
                            *self._cofts[self.cpointer[k,k],it]

                return ret
            else:
                ret = numpy.zeros((nos,nos,nos,nos),dtype=COMPLEX)
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
        nos = self.nob + 1
        if self._is_transformed:
            if t is None:
                Nt = self.max_cutoff_index + 1
                
                # FIXME: some consistent treatment of the cutoff time is needed  
                Nt = self._cofts.shape[1]
                              
                ret = numpy.zeros((nos,nos,nos,nos,Nt),dtype=COMPLEX)
                #for k in range(self.nof):
                for k in range(self.nob):
                    for it in range(Nt):
                        #ret[:,:,:,:,it] += self._A4[:,:,:,:,k]\
                        #*self._hofts[k+1,it]
                        ret[:,:,:,:,it] += self._C4[:,:,:,:,k]\
                            *self._hofts[self.cpointer[k,k],it]
                return ret
            else:
                ret = numpy.zeros((nos,nos,nos,nos),dtype=COMPLEX)
                it = self.timeAxis.nearest(t)
                for k in range(self.nof):
                    ret += self._A4[:,:,:,:,k]*self._hofts[k+1,it]
                return ret
        else:
            raise Exception()


    def get_goft4(self,a,b,c,d):
        """Returns the matrix of the goft functions
        
        """
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
        nos = self.nob + 1
        if self._is_transformed:
            if t is None:
                Nt = self.max_cutoff_index + 1
                
                # FIXME: some consistent treatment of the cutoff time is needed
                Nt = self._cofts.shape[1]
                
                ret = numpy.zeros((nos,nos,nos,nos,Nt),dtype=COMPLEX)
                #for k in range(self.nof):
                for k in range(self.nob):
                    for it in range(Nt):
                        ret[:,:,:,:,it] += self._C4[:,:,:,:,k]\
                            *self._gofts[self.cpointer[k,k],it]
                        #ret[:,:,:,:,it] += self._A4[:,:,:,:,k]\
                        #*self._gofts[k+1,it]
                return ret
            else:
                ret = numpy.zeros((nos,nos,nos,nos),dtype=COMPLEX)
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
        """Get the index of the correlation function corresponding a tuple where
        
        """
        ret = -1
        for i in range(self.nof+1):
            if where in self.where[i]:
                ret = i

        return ret


    def _update_where(self, iof, fce, where):
        """Updates the location information for a correlation function
        
        """

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
        """Calculates a time integral of all correlation functions
        
        """
        if not self._has_hofts:
            self._hofts = numpy.zeros((self.nof+1,self.timeAxis.length),
                                  dtype=numpy.complex128)
            for ii in range(self.nof+1):
                self._hofts[ii,:] = c2h(self.timeAxis,self._cofts[ii,:])
            self._has_hofts = True


    def create_double_integral(self):
        """Calculates a double time integral of all correlation functions
        
        """
        if not self._has_gofts:
            self._gofts = numpy.zeros((self.nof+1,self.timeAxis.length),
                                      dtype=numpy.complex128)
            for ii in range(self.nof+1):
                self._gofts[ii,:] = c2g(self.timeAxis,self._cofts[ii,:])
            self._has_gofts = True


    #
    #
    #   TO BE TRANSPLANTED TO OPEN SYSTEMS
    #  
    #
    def init_site_mapping(self):
        """Initializes the internals to allow four-index correlation functions
        
        """
        
        nos = self.nob + 1
        one = numpy.eye(nos, dtype=REAL)
        
        self.transform(one)


    #
    #
    #   TO BE TRANSPLANTED TO OPEN SYSTEMS
    #  
    #
    def transform(self, SS):
        """Transform the system-bath interaction characteristics
        
        
        """
        
        # FIXME: This must go somewhere alse. It is system dependent
        
        # system size is by one larger than the number of sites (excitonic system)
        nos = self.nob + 1
        nof = self.nof
        
        self._A4 = numpy.zeros((nos,nos,nos,nos,nof), dtype=REAL)
        self._C4 = numpy.zeros((nos,nos,nos,nos,nos), dtype=REAL)

        for a in range(nos):
            for b in range(nos):
                for c in range(nos):
                    for d in range(nos):
                        for n in range(nos):
                            self._C4[a,b,c,d,n] = \
                                SS[n,a]*SS[n,b]*SS[n,c]*SS[n,d]
                            
                        for k in range(nof):

                           for n in range(1,nos):
                                for m in range(1,nos):
                                    self._A4[a,b,c,d,k] += SS[n,a]*SS[n,b]*\
                                    self._A2[n-1,m-1,k+1]*SS[m,c]*SS[m,d]
                                

        self._is_transformed = True

