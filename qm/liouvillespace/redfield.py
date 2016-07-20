# -*- coding: utf-8 -*-

import numpy

from ...core.implementations import implementation
from ...core.units import cm2int
from ...core.units import kB_intK

from ..hilbertspace.hamiltonian import Hamiltonian
from ..liouvillespace.systembathinteraction import SystemBathInteraction

        
        
class RedfieldRateMatrix:
    
    def __init__(self,ham,sbi,initialize=True,cutoff_time=None):
        
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
        
        if Nk <= 0:
            raise Exception("No system bath intraction components present")
        
        # Eigen problem
        hD,SS = numpy.linalg.eigh(ham.data) 
        S1 = numpy.linalg.inv(SS)
        
        # component operators
        KI = self.sbi.KK
        
        # frequency cut-off
        freq_cutoff = 3000.0*cm2int
        print("freq. cut-off",freq_cutoff)
        
        # temperature
        Temp = self.sbi.CC.get_correlation_function(0,0).T
        



        """ Might need some speed-up too"""
        
        # transform interaction operators
        for i in range(Nk):
            KI[i,:,:] = numpy.dot(S1,numpy.dot(KI[i,:,:],SS))
            
        #  Find all eigenfrequencies 
        Om = numpy.zeros((Na,Na))
        for a in range(0,Na):
            for b in range(0,Na):
                Om[a,b] = hD[a] - hD[b]
                
        # calculate values of the spectral density at frequencies
        cc = numpy.zeros((Nk,Na,Na),dtype=numpy.complex)
        
        # loop over components
        for k in range(Nk):

            # spectral density
            cf = self.sbi.CC.get_correlation_function(k,k)
            cf.convert_2_spectral_density()

            # Spectral density at all frequencies
            for i in range(Na):
                for j in range(Na):
                    if i != j:
                        if numpy.abs(Om[i,j]) > freq_cutoff:
                            cc[k,i,j] = 0.0
                        else:
                            if Om[j,i] < 0.0:
                                cc[k,i,j] = (cf.interp_data(Om[i,j])
                                *numpy.exp(Om[j,i]/(kB_intK*Temp)))
                            else:
                                cc[k,i,j] = cf.interp_data(Om[j,i])
                                

        """ To submit: 
                        Na
                        Nk
                        KI[Nk,Na,Na]
                        cc[Nk,Na,Na]
                        
             
            To return:
                        RR
                        
        """
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        #print(KI[0,:,:])
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        warning = []
        error = []
        rtol = 1.0e-6
        self.data = ssRedfieldRateMatrix(Na,Nk,KI,cc,rtol,warning,error)
        for w in error:
            print(w)                                    
        self._is_initialized = True
        
        

@implementation("redfield","ssRedfieldRateMatrix",
                at_runtime=True,
                fallback_local=True,
                always_local=True)
def ssRedfieldRateMatrix(Na,Nk,KI,cc,rtol,warning,error):
    
    # output relaxatio rate matrix
    RR = numpy.zeros((Na,Na),dtype=numpy.float)
    cc = 2.0*numpy.real(cc)
    
    # loop over components
    for k in range(Nk):
        
        # interaction operator
        KK = KI[k,:,:]

        for i in range(Na):
            for j in range(Na):
                if i != j:                                
                            
                    RR[i,j] += (cc[k,i,j]*KK[i,j]*KK[j,i])
                    
                    if RR[i,j] < 0.0:
                        warning.append("\n%i %i %r smaller than zero"
                                    % (i,j,RR[i,j]))
                        if numpy.abs(RR[i,j]) < rtol:
                            RR[i,j] = 0.0
                        else:
                            error.append(" ... significantly")
                            
                    RR[j,j] -= RR[i,j]
                    
    return RR
                    
                    