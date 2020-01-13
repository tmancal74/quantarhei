# -*- coding: utf-8 -*-

"""
*******************************************************************************

      MODIFIED REDFIELD RATE MATRIX

*******************************************************************************
"""  

import numpy
from scipy import integrate

from quantarhei.core.implementations import implementation
from quantarhei.core.units import cm2int
from quantarhei.core.units import kB_intK

from quantarhei.qm.hilbertspace.hamiltonian import Hamiltonian
from quantarhei.qm.liouvillespace.systembathinteraction import SystemBathInteraction

import quantarhei as qr
import itertools as it

class ModifiedRedfieldRateMatrix:
    """Modifield Redfield relaxation rate matrix
    
    Modified Redfield population relaxation rate matrix is calculated from the 
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
    
    def __init__(self, ham, sbi, time, initialize=True, cutoff_time=None):
        
        if not isinstance(ham,Hamiltonian):
            raise Exception("First argument must be a Hamiltonian")
            
        if not isinstance(sbi,SystemBathInteraction):
            raise Exception
            
        self._is_initialized = False            
        self._has_cutoff_time = False
        
        if cutoff_time is not None:
            self.cutoff_time = cutoff_time
            self._has_cutoff_time = True            
            
        self.ham = ham
        self.sbi = sbi
        self.tt = time.data
        
        if initialize: 
            self._set_rates()          
            self._is_initialized = True
            
    def _set_rates(self):
        """
        
        """

        Na = self.ham.dim
        Nc = self.sbi.N 
        tt = self.tt
    
        # Eigen problem
        hD,SS = numpy.linalg.eigh(self.ham._data) #hD=eigenvalues, SS=eigenvectors
        
        Nt = self.sbi.CC.timeAxis.length
        
        lam4 = numpy.zeros((Na-1,Na-1,Na-1,Na-1),dtype=qr.REAL)
        #lam4 = numpy.zeros((Na,Na,Na,Na),dtype=qr.REAL)
        
        for a in range(Na-1):
        #for a in range(1,Na):    
            for b in range(Na-1):
            #for b in range(1,Na):
                for c in range(Na-1):
                #for c in range(1,Na):   
                    for d in range(Na-1):
                    #for d in range(1,Na):   
                         lam4[a,b,c,d] = self.sbi.CC.get_reorganization_energy4(a,b,c,d)
            
        self.sbi.CC.create_double_integral() #g(t)
        self.sbi.CC.create_one_integral()  #g_dot(t)
        self.sbi.CC.transform(SS)
        #g4_1value = self.sbi.CC.get_goft4(1,2,3,4)
        g4 = self.sbi.CC.get_goft_matrix()   #g_{abcd}(t), dimensions (Na, Na, Na, Na, Nt-1)
        h4 = self.sbi.CC.get_hoft_matrix()   #g_dot_{abcd}(t), dimensions (Na, Na, Na, Na, Nt-1)
        c4 = self.sbi.CC.get_coft_matrix()   #g_dotdot_{abcd}(t) = C(t) in exciton basis, dimensions (Na, Na, Na, Na, Nt-1)
        
        rates = ssModifiedRedfieldRateMatrix(Na, Nc, Nt, hD, lam4, g4, h4, c4, tt)
        self.rates = rates

def ssModifiedRedfieldRateMatrix(Na, Nc, Nt, hD, lam4, g4, h4, c4, tt): #, Ee, SS, prt, gg, hh, cf, tt, ls,
                                 #rtol, werror, RR):
        """Standard redfield rates
    
    
        Parameters
        ----------
    
        Na : integer
        Rank of the rate matrix, number of excitons
        
        Nc : integer
        Number of components of the interaction Hamiltonian (number of sites
        or number of distinct correlation functions)
    
        Ee : float array
        Eigen energies of the Hamiltonian
        
        SS : float array
        Transformation matrix 
        (components of the interaction operator in the exciton basis)
        
        prt : integer array
        Pointer between site index and correlation function index
    
        gg : float array
        Line shape functions
        
        hh : float array
        derivative of the line shape function
    
        cf : float array
        second derivatives of the line shape functions (correlation functions)
        
        tt : float array
        values of time
        
        ls : float array
        reorganization energies of the sites
        
        RR : real array
        Relaxation rate matrix (to be calculated and returned)
    
        rtol : float array
        tolerances
        
        werror : integer array
        warnings and errors
        
        """
        
        print("***")
        
        E_0k = numpy.zeros(Na-1,dtype=numpy.float)
        #E_0k = numpy.zeros(Na,dtype=numpy.float)
        
        for ii in range(Na-1):
        #for ii in range(1,Na):
            E_0k[ii] = hD[ii] - lam4[ii,ii,ii,ii]  
            
        F_k_t = numpy.zeros((Na-1,Nt-1),dtype=numpy.complex)    
        A_k_t = numpy.zeros((Na-1,Nt-1),dtype=numpy.complex) 
        N_kl_t = numpy.zeros((Na-1,Na-1,Nt-1),dtype=numpy.complex) 
        
        #F_k_t = numpy.zeros((Na,Nt-1),dtype=numpy.complex)    
        #A_k_t = numpy.zeros((Na,Nt-1),dtype=numpy.complex) 
        #N_kl_t = numpy.zeros((Na,Na,Nt-1),dtype=numpy.complex) 
        
        for a in range(Na-1):
        #for a in range(1,Na):
            for ti in range(Nt-1):
                F_k_t[a,ti] = numpy.exp(-1j*(E_0k[a] - lam4[a,a,a,a])*tt[ti] - numpy.conjugate(g4[a,a,a,a,ti])) 
                A_k_t[a,ti] = numpy.exp(-1j*(E_0k[a] + lam4[a,a,a,a])*tt[ti] - g4[a,a,a,a,ti]) 
                                
        for a in range(Na-1):
        #for a in range(1,Na):
            for b in range(Na-1):
            #for b in range(1,Na):
                for ti in range(Nt-1): 
                    N_kl_t[a,b,ti] = (c4[b,a,a,b,ti] - (h4[b,a,a,a,ti] -  h4[b,a,b,b,ti] - 2j*lam4[b,a,b,b])*(h4[a,b,a,a,ti] - h4[a,b,b,b,ti] - 2j*lam4[a,b,b,b]))*numpy.exp(2*(g4[a,a,b,b,ti] + 1j*lam4[a,a,b,b]*tt[ti]))
                  
        f = numpy.zeros((Na-1,Na-1,Nt),dtype=numpy.complex) 
        RR = numpy.zeros((Na-1,Na-1),dtype=numpy.complex) 
        
        #f = numpy.zeros((Na,Na,Nt),dtype=numpy.complex) 
        #RR = numpy.zeros((Na,Na),dtype=numpy.complex) 
        
        for a in range(Na-1):
        #for a in range(1,Na):
            for b in range(Na-1):
            #for b in range(1,Na):
                for ti in range(Nt-1):
                    f[a,b,ti] = numpy.conjugate(F_k_t[b,ti])*A_k_t[a,ti]*N_kl_t[a,b,ti]
                    RR[a,b] = 2*numpy.real(integrate.simps(f[a,b,:],tt))
                    
        for a in range(Na-1):
        #for a in range(1,Na):
            RR[a,a] = 0
                   
        RR_bbaa = -numpy.sum(RR, axis = 0)
        print('RR_bbaa is:')
        print(RR_bbaa)
        #RR = numpy.diag(numpy.diag(RR_bbaa)) + RR  
        RR = numpy.diag(RR_bbaa) + RR  
        print('diag RR_bbaa is:')
        print(numpy.diag(RR_bbaa))
        print('RR is:')
        print(RR)
                        
        print("I am called from outside")
        return RR

        qr.stop()
                    
                    
