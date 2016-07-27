# -*- coding: utf-8 -*-

import numpy
import scipy

from .systembathinteraction import SystemBathInteraction
from ..hilbertspace.hamiltonian import Hamiltonian

"""
*******************************************************************************

    REDFIELD RELAXATION TENSOR

*******************************************************************************
"""
#FIXME: This should go somewhere else
class RelaxationTensor:
    pass

class RedfieldRelaxationTensor(RelaxationTensor):
    
    def __init__(self,ham,sbi,initialize=True,cutoff_time=None):
    
        """Constructs the Redfield Relaxation Tensor
        
        
        Parameters
        ----------
    
        ham : cu.oqs.hilbertspace.Hamiltonian
            Hamiltonian of the system.
        
        sbi : cu.oqs.liouvillespace.SystemBathInteraction
            Object specifying the system-bath interaction
        
        initialize : bool
            If True, the tensor will be imediately calculated
        
        cutoff_time : float
            Time in femtoseconds after which the integral kernel in the definition
            of the relaxation tensor is assummed to be zero. 
        """
    
        if not isinstance(ham,Hamiltonian):
            raise Exception("First argument must be a Hamiltonian")
            
        if not isinstance(sbi,SystemBathInteraction):
            raise Exception

        self._is_initialized = False            
        self._has_cutoff_time = False
        if cutoff_time is not None:
            self.cutoff_time = cutoff_time
            self._has_cutoff_time = True            
            
        self.Hamiltonian = ham
        self.BilinearSystemBathInteraction = sbi
        
        if initialize:        
            self.__reference_implementation(ham,sbi)
            #dtn = Rrt.cython_implementation(self,ham,sbi)
            #self.data = dtn            
            self._is_initialized = True

        
    def __reference_implementation(self,ham,sbi):
        """ Reference implementation, completely in Python
        
        """
        
        # dimension of the Hamiltonian (includes excitons
        # with all multiplicities specified at its creation)
        Na = ham.data.shape[0]
        ta = sbi.TimeAxis
        multi_ex = False
        
        if sbi.aggregate is not None:
            agg = sbi.aggregate
            if agg.mult > 1:
                multi_ex = True
            
        
        if self._has_cutoff_time:
            tcut = ta.nearest(self.cutoff_time)
            tm = ta.time[0:tcut]
            length = tcut
        else:
            tm = ta.time
            length = ta.length

        # Eigen problem
        hD,SS = numpy.linalg.eigh(ham.data)   
            
        #  All eigenfrequencies 
        Om = numpy.zeros((Na,Na))
        for a in range(0,Na):
            for b in range(0,Na):
                Om[a,b] = hD[a] - hD[b]
                
        # Integrals of correlation functions from the set 
        #
        # length of the set is the number of bath == number of sites
        #Nb = sbi.CC.nob 
        #print(multi_ex)
        
        if not multi_ex:
            #Nb = Na
            Nb = sbi.N
            hset = numpy.zeros((Nb,Na,Na),dtype=numpy.complex128)
            
            for ns in range(0,Nb):
                
                rc1 = sbi.get_coft(ns,ns) 
                
                for a in range(0,Na):
                    for b in range(0,Na):
                        eexp = numpy.exp(1.0j*Om[a,b]*tm) 
                        rc = rc1[0:length]*eexp
                        rr = numpy.real(rc)
                        ri = numpy.imag(rc)
                        sr = scipy.interpolate.UnivariateSpline(tm,
                                    rr,s=0).antiderivative()(tm)
                        si = scipy.interpolate.UnivariateSpline(tm,
                                    ri,s=0).antiderivative()(tm)
                    
                        hset[ns,a,b] = sr[length-1] + 1.0j*si[length-1]
        
        else:
            
            #Nb = agg.Nel
            Nb = sbi.N
            hset = numpy.zeros((Nb,Nb,Na,Na),dtype=numpy.complex128)
            ns = 0
            for ens in agg.elsignatures(mult=agg.mult):
                print(ns)
                ms = 0
                for ems in agg.elsignatures(mult=agg.mult):
                
                    rc1 = sbi.get_coft_elsig(ens,ems)

                    for a in range(0,Na):
                        for b in range(0,Na):

                            eexp = numpy.exp(1.0j*Om[a,b]*tm) 
                            rc = rc1[0:length]*eexp
                            rr = numpy.real(rc)
                            ri = numpy.imag(rc)
                            sr = scipy.interpolate.UnivariateSpline(tm,
                                    rr,s=0).antiderivative()(tm)
                            si = scipy.interpolate.UnivariateSpline(tm,
                                    ri,s=0).antiderivative()(tm)
                    
                            hset[ns,ms,a,b] = (sr[length-1] +
                                   1.0j*si[length-1])
                            
                    ms += 1
                ns += 1
        
         
        # Site operators 
        K = numpy.zeros((Nb,Na,Na),dtype=numpy.complex128)
        S1 = scipy.linalg.inv(SS)
               
        
        for x in range(0,Nb):
            K[x,:,:] = sbi.KK[x,:,:]
            K[x,:,:] = numpy.dot(S1,numpy.dot(K[x,:,:],SS))
            
        
        if multi_ex:         
         
            """ Coefficients of the tensor """
            Hh = numpy.zeros((Na,Na,Na,Na),dtype=numpy.complex128)
            Hc = numpy.zeros((Na,Na,Na,Na),dtype=numpy.complex128)
            for a in range(0,Na):
                print(a)
                for b in range(a,Na):
                    for c in range(0,Na):
                        for d in range(0,Na):
                            # sum over sites
                            for x in range(0,Nb):
                                for y in range(0,Nb):
                                    xs = x
                                    ys = y
                                    Hh[a,b,c,d] += (K[x,a,b]*K[y,c,d]
                                        *hset[xs,ys,c,d]) 
                                    Hc[a,b,c,d] += (
                                        numpy.conj(K[x,a,b]*K[y,c,d]
                                        *hset[xs,ys,d,c]))
                            if a != b:
                                Hh[b,a,c,d] = Hh[a,b,c,d]
                                Hc[b,a,c,d] = Hc[a,b,c,d]

        else:
            
            """ Coefficients of the tensor """
            Hh = numpy.zeros((Na,Na,Na,Na),dtype=numpy.complex128)
            Hc = numpy.zeros((Na,Na,Na,Na),dtype=numpy.complex128)
            for a in range(0,Na):
                for b in range(a,Na):
                    for c in range(0,Na):
                        for d in range(0,Na):
                            # sum over sites
                            for x in range(0,Nb):
                                xs = x
                                Hh[a,b,c,d] += (K[x,a,b]*K[x,c,d]
                                        *hset[xs,c,d]) 
                                Hc[a,b,c,d] += (numpy.conj(K[x,a,b]*K[x,c,d]
                                        *hset[xs,d,c]))
                            if a != b:
                                Hh[b,a,c,d] = Hh[a,b,c,d]
                                Hc[b,a,c,d] = Hc[a,b,c,d]
            
            
        """ Preliminary sums """
        Hrh = numpy.zeros((Na,Na),dtype=numpy.complex128)
        Hrc = numpy.zeros((Na,Na),dtype=numpy.complex128)
        Hrh = numpy.einsum('ijjk',Hh)
        Hrc = numpy.einsum('ijki',Hc)

       
        """ Redfield tensor """
        RR = numpy.zeros((Na,Na,Na,Na),dtype=numpy.complex128)
        
        for a in range(0,Na):
            for b in range(0,Na):
                for c in range(0,Na):
                    for d in range(0,Na):
                        RR[a,b,c,d] = -Hc[a,c,d,b]-Hh[d,b,a,c]


        for a in range(0,Na):
            for b in range(0,Na):
                for c in range(0,Na):
                    #for d in range(0,Na):   
                    RR[a,b,c,b] += Hrh[a,c]
                    RR[b,a,b,c] += Hrc[a,c]

        
        self.data = -RR
        self._is_initialized = True


    def initialize(self):
        self.__reference_implementation(self.Hamiltonian,
                                        self.BilinearSystemBathInteraction)
        
    def secularize(self):
        
        N = self.data.shape[0]
        for ii in range(0,N):
            for jj in range(0,N):
                for kk in range(0,N):
                    for ll in range(0,N):
                        if not (((ii == jj) and (kk == ll)) 
                        or ((ii == kk) and (jj == ll))) :
                            self.data[ii,jj,kk,ll] = 0
                                        
