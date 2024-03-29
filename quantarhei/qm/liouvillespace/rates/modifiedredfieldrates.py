# -*- coding: utf-8 -*-

"""
*******************************************************************************

      MODIFIED REDFIELD RATE MATRIX

*******************************************************************************
"""  

import numpy
from scipy import integrate

#from quantarhei.core.implementations import implementation

from quantarhei.qm.hilbertspace.hamiltonian import Hamiltonian
from quantarhei.qm.liouvillespace.systembathinteraction \
    import SystemBathInteraction

from quantarhei import REAL
from quantarhei import COMPLEX

#import itertools as it

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
            
        self.ham = ham
        self.sbi = sbi
        
        if initialize: 
            self._set_rates()          
            self._is_initialized = True
            
            
    def _set_rates(self):
        """ Setting Modified Redfield rates for an electronic system
        
        We assume a single exciton band only!!!
        
        
        
        """

        Na = self.ham.dim
        Nc = self.sbi.N 
        
        # this is the system for which we calculate rates (OpenSystem) 
        sys = self.sbi.system
        
        # The size of the Hamiltonian must be by 1 larger than
        # the number of sites - this is the excitonic situation
        if Na != Nc + 1:
            raise Exception("Modified Redfield is implemented for "
                            +"single exciton systems only")        
        
        tt = self.sbi.TimeAxis.data
        Nt = self.sbi.CC.timeAxis.length
    
        # Eigen problem 
        hD,SS = numpy.linalg.eigh(self.ham._data) 
        
        # 
        # THIS WILL BE DONE ON THE OPEN SYSTEM
        #
        self.sbi.CC.transform(SS)
        
        #
        # HIDE THIS INTO: get_reorganization_energy_matrix()
        #
        lam4 = numpy.zeros((Na,Na,Na,Na),dtype=REAL)
        for a in range(Na):
            for b in range(Na):
                for c in range(Na):
                    for d in range(Na):
                        lam4[a,b,c,d] = \
                             self.sbi.CC.get_reorganization_energy4(a,b,
                                                                    c,d)
        
        #
        #  THESE CALLS HAVE TO BE DONE AUTOMATICALY BEFORE get_Xoft_matrix()
        #
        self.sbi.CC.create_double_integral() #g(t)
        self.sbi.CC.create_one_integral()  #g_dot(t)
        
        # g_{abcd}(t), dimensions (Na, Na, Na, Na, Nt-1)
        g4 = self.sbi.CC.get_goft_matrix()
        
        # g_dot_{abcd}(t), dimensions (Na, Na, Na, Na, Nt-1)
        h4 = self.sbi.CC.get_hoft_matrix()  
        
        # g_dotdot_{abcd}(t) = C(t) in exciton basis, 
        # dimensions (Na, Na, Na, Na, Nt-1)
        c4 = self.sbi.CC.get_coft_matrix()   


        rates = ssModifiedRedfieldRateMatrix(Na, Nc, Nt, hD,
                                             lam4, g4, h4, c4, tt)
        
        self.data = rates


def ssModifiedRedfieldRateMatrix(Na, Nc, Nt, hD, lam4, g4, h4, c4, tt): 
                                 #rtol, werror, RR):
        """Modified redfield rates
    
    
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
        
        E_0k = numpy.zeros(Na,dtype=REAL)
        
        for ii in range(Na):
            E_0k[ii] = hD[ii] - lam4[ii,ii,ii,ii]  
            
        F_k_t = numpy.zeros((Na,Nt),dtype=COMPLEX)    
        A_k_t = numpy.zeros((Na,Nt),dtype=COMPLEX) 
        N_kl_t = numpy.zeros((Na,Na,Nt),dtype=COMPLEX) 
        
        for a in range(Na):
            for ti in range(Nt):
                F_k_t[a,ti] = (numpy.exp(-1j*(E_0k[a] - lam4[a,a,a,a])*tt[ti]
                                         - numpy.conjugate(g4[a,a,a,a,ti])) )
                A_k_t[a,ti] = (numpy.exp(-1j*(E_0k[a] + lam4[a,a,a,a])*tt[ti]
                                         - g4[a,a,a,a,ti]) )
                                
        for a in range(Na):
            for b in range(Na):
                for ti in range(Nt): 
                    N_kl_t[a,b,ti] = ((c4[b,a,a,b,ti] - (h4[b,a,a,a,ti]
                                     -  h4[b,a,b,b,ti]
                                     - 2j*lam4[b,a,b,b])*(h4[a,b,a,a,ti] 
                                     - h4[a,b,b,b,ti] - 2j*lam4[a,b,b,b]))
                                     *numpy.exp(2*(g4[a,a,b,b,ti]
                                        + 1j*lam4[a,a,b,b]*tt[ti])))
                  
        f = numpy.zeros((Na,Na,Nt),dtype=COMPLEX) 
        RR = numpy.zeros((Na,Na),dtype=COMPLEX) 
        
        for a in range(Na):
            for b in range(Na):

                for ti in range(Nt):
                    f[a,b,ti] = (numpy.conjugate(F_k_t[b,ti])
                                 *A_k_t[a,ti]*N_kl_t[a,b,ti])
                    RR[a,b] = 2*numpy.real(integrate.simps(f[a,b,:],tt))
                    
        for a in range(Na):
            RR[a,a] = 0
                   
        RR_bbaa = -numpy.sum(RR, axis = 0)

        RR = numpy.diag(RR_bbaa) + RR  

        return RR

                    
                    
