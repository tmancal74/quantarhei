# -*- coding: utf-8 -*-
"""
*******************************************************************************

      REDFIELD RATE MATRIX

*******************************************************************************
"""  

import numpy
import scipy

from ....core.implementations import implementation
from ....core.units import cm2int

from ...hilbertspace.hamiltonian import Hamiltonian
from ...liouvillespace.systembathinteraction import SystemBathInteraction

from ....core.time import TimeDependent

   
class TDRedfieldRateMatrix(TimeDependent):
    """Redfield relaxation rate matrix
    
    Redfield population relaxation rate matrix is calculated from the 
    Hamiltonian and system-system bath interation. The bath
    correlation functions are Fourier transformed by FFT and the negative
    frequency part is calculated from the expected thermodynamic
    symmetry. 
    
    Parameters
    ----------
    
    ham : Hamiltonian
        Hamiltonian object 
        
    sbi : SystemBathInteraction
        SystemBathInteraction object
        
    initialize : bool (default True)
        If true, the rates will be calculated when the object is created
        
    cutoff_time : float
        If cutoff time is specified, the tensor is integrated only up to the
        cutoff time
    
    
    """
    
    def __init__(self, ham, sbi, initialize=True, cutoff_time=None):
        
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
        self.sbi = sbi
        
        if initialize: 
            self._set_rates(ham,sbi)          
            self._is_initialized = True
                
                
    def _set_rates(self,ham,sbi):
        """ Reference implementation, completely in Python
        
        """
        
        # dimension of the Hamiltonian (includes excitons
        # with all multiplicities specified at its creation)
        Na = ham.data.shape[0]
     
        # number of components
        Nk = self.sbi.N  
        
        # number of steps in time
        time = self.sbi.TimeAxis
        tm = time.data
        Nt = time.length
        
        if Nk <= 0:
            raise Exception("No system bath intraction components present")
        
        # Eigen problem
        hD,SS = numpy.linalg.eigh(ham.data) 
        S1 = numpy.linalg.inv(SS)
        
        # component operators
        KI = self.sbi.KK.copy()
        
        #FIXME: This has to be made configurable
        # frequency cut-off
        freq_cutoff = 3000.0*cm2int
        #print("freq. cut-off",freq_cutoff)
     
        
        # transform interaction operators
        for i in range(Nk):
            KI[i,:,:] = numpy.dot(S1,numpy.dot(KI[i,:,:],SS))
            
        #  Find all eigenfrequencies 
        Om = numpy.zeros((Na,Na))
        for a in range(0,Na):
            for b in range(0,Na):
                Om[a,b] = hD[a] - hD[b]
                
        # calculate values of the spectral density at frequencies
        cc = numpy.zeros((Nt,Nk,Na,Na),dtype=numpy.complex)
        
        # loop over components
        for k in range(Nk):

            # correlation function
            cf = self.sbi.CC.get_correlation_function(k,k)
            
            # Ft correlation function
            # cw = cf.get_Fourier_transform()
            

            # Spectral density at all frequencies
            for i in range(Na):
                for j in range(Na):
                    if i != j:
                        if numpy.abs(Om[j,i]) > freq_cutoff:
                            cc[:,k,i,j] = 0.0
                        else:
                            ff = cf.data*numpy.exp(1.0j*Om[j,i]*tm)
                            rr = numpy.real(ff)
                            ri = numpy.imag(ff)
                            sr = scipy.interpolate.UnivariateSpline(tm,
                                    rr, s=0).antiderivative()(tm)
                            si = scipy.interpolate.UnivariateSpline(tm,
                                    ri, s=0).antiderivative()(tm)
                            cc[:,k,i,j] =  sr + 1.0j*si
                                

        """ To submit: 
                        Na
                        Nk
                        KI[Nk,Na,Na]
                        cc[Nk,Na,Na]
                        
             
            To return:
                        RR
                        
        """
        
        warning = []
        error = []
        rtol = 1.0e-6
        self.data = ssTDRedfieldRateMatrix(Na, Nk, Nt, KI, cc, rtol,
                                           warning, error)
        for w in error:
            print(w)     
                               
        self._is_initialized = True
        
        

@implementation("tdredfield","ssTDRedfieldRateMatrix",
                at_runtime=True,
                fallback_local=True,
                always_local=True)
def ssTDRedfieldRateMatrix(Na, Nk, Nt, KI, cc, rtol, warning, error):
    
    # output relaxatio rate matrix
    RR = numpy.zeros((Nt,Na,Na),dtype=numpy.float)
    
    # real part of FT correlation function = 2Re of the half FT of 
    # the correlation function - no factor of 2 here!
    cc = numpy.real(cc)
    
    # loop over components
    for k in range(Nk):
        
        # interaction operator
        KK = KI[k,:,:]

        for i in range(Na):
            for j in range(Na):
                
                # calculate rates, i.e. off diagonal elements
                if i != j:                                
                            
                    RR[:,i,j] += 2.0*numpy.real(cc[:,k,i,j]*KK[i,j]*KK[j,i])
                    
    # calculate the diagonal elements (the depopulation rates)            
    for i in range(Na):
        for j in range(Na):
            if i != j:
                RR[:,j,j] -= RR[:,i,j]
                    
    return RR
                    
                    

