# -*- coding: utf-8 -*-

"""
*******************************************************************************

      MODIFIED REDFIELD RATE MATRIX

*******************************************************************************
"""  

import numpy

from ....core.implementations import implementation
from ....core.units import cm2int
from ....core.units import kB_intK

from ...hilbertspace.hamiltonian import Hamiltonian
from ...liouvillespace.systembathinteraction import SystemBathInteraction

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
    p
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
        """
        
        """
        Na = self.ham.dim
        Nc = self.sbi.N 
    
        # Eigen problem
        hD,SS = numpy.linalg.eigh(self.ham._data) #hD=eigenvalues, SS=eigenvectors

        # loop over components
        for k in range(Nc):

            # correlation function
            cf = self.sbi.CC.get_correlation_function(k,k)
            
        Nt = cf.axis.length
        
        ssModifiedRedfieldRateMatrix(Na, Nc, Nt, hD, SS)


def ssModifiedRedfieldRateMatrix(Na, Nc, Nt, hD, SS): #, Ee, SS, prt, gg, hh, cc, tt, ls,
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
    
        cc : float array
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
   
        ls_n = numpy.zeros(Na)
        ls_n[0] = 21.2799
        ls_n[1] = 31.5194
        ls_n[2] = 22.8585
        ls_n[3] = 17.8808
        ls_n[4] = 15.3586
        ls_n[5] = 23.1838
        #ls_n(6) = 24.8873
        #ls_n(7) = 40.9997

        RR = numpy.zeros((Na,Na),dtype=numpy.float)
        m_ket = numpy.identity(Na)
    
        g_a = numpy.zeros((Nt,Na),dtype=numpy.float)
        ls_a = numpy.zeros((Na,1),dtype=numpy.float)

        for a in range(Na):
            for m in range(Na):
        #        g_a_temp = (m_ket[:,m].reshape(-1,1)*SS[:,a])**4*g_t[:,m]        
        #        g_a[:,a] = g_a[:,a] + g_a_temp             
        #        ls_a_temp = (m_ket[:,m]*SS[:,a].reshape(-1,1))**4*ls_n[m]
                 ls_a_temp = (m_ket[:,m]*SS[:,a].reshape(-1,1))
        #        ls_a_temp = 0
                 ls_a[a] = ls_a[a] + ls_a_temp
            
        E_0k = hD - ls_a  #a[alpha,beta] = numpy.conj(c)*c

        pairs = numpy.array(list(it.combinations(range(1,Na+1),2)))

        g_baab_dot2 = numpy.zeros((Nt,len(pairs[:,1])),dtype=numpy.float)
        g_abaa_dot = numpy.zeros((Nt,len(pairs[:,1])),dtype=numpy.float)
        g_babb_dot = numpy.zeros((Nt,len(pairs[:,1])),dtype=numpy.float)
        g_abbb_dot = numpy.zeros((Nt,len(pairs[:,1])),dtype=numpy.float)
        g_baaa_dot = numpy.zeros((Nt,len(pairs[:,1])),dtype=numpy.float)
        g_aabb = numpy.zeros((Nt,len(pairs[:,1])),dtype=numpy.float)
        ls_babb = numpy.zeros((1,len(pairs[:,1])),dtype=numpy.float)
        ls_abaa = numpy.zeros((1,len(pairs[:,1])),dtype=numpy.float)
        ls_abbb = numpy.zeros((1,len(pairs[:,1])),dtype=numpy.float)
        ls_baaa = numpy.zeros((1,len(pairs[:,1])),dtype=numpy.float)
        ls_aabb = numpy.zeros((1,len(pairs[:,1])),dtype=numpy.float)

        for p in range(len(pairs[:,1])):
            for m in range(Na):
                 a = pairs[p,0]  #smaller number
                 b = pairs[p,1]  #larger number
                 psi_a_mv = m_ket[:,m]*SS[:,a-1].reshape(-1,1)  #dimensionless, expansion (superposition) coeff 
                 psi_b_mv = m_ket[:,m]*SS[:,b-1].reshape(-1,1)  #dimensionless
              
        #         g_baab_dot2_temp = psi_b_mv*psi_a_mv*psi_a_mv*psi_b_mv*g_t_dot2[:,m]
        #         g_baab_dot2[:,p] = g_baab_dot2[:,p] + g_baab_dot2_temp
                        
        #         g_baaa_dot_temp = psi_b_mv*psi_a_mv*psi_a_mv*psi_a_mv*g_t_dot[:,m]
        #         g_baaa_dot[:,p] = g_baaa_dot[:,p] + g_baaa_dot_temp         
             
        #         g_babb_dot_temp = psi_b_mv*psi_a_mv*psi_b_mv*psi_b_mv*g_t_dot[:,m]
        #         g_babb_dot[:,p] = g_babb_dot[:,p] + g_babb_dot_temp    
             
        #         g_abaa_dot_temp = psi_a_mv*psi_b_mv*psi_a_mv*psi_a_mv*g_t_dot[:,m]
        #         g_abaa_dot[:,p] = g_abaa_dot[:,p] + g_abaa_dot_temp         
             
        #         g_abbb_dot_temp = psi_a_mv*psi_b_mv*psi_b_mv*psi_b_mv*g_t_dot[:,m]
        #         g_abbb_dot[:,p] = g_abbb_dot[:,p] + g_abbb_dot_temp   
                                      
        #         g_aabb_temp = psi_a_mv*psi_a_mv*psi_b_mv*psi_b_mv*g_t[:,m]
        #         g_aabb[:,p] = g_aabb[:,p] + g_aabb_temp
          
        #         ls_babb_temp = psi_b_mv*psi_a_mv*psi_b_mv*psi_b_mv*ls_n[m]
        #         ls_babb[p] = ls_babb[p] + ls_babb_temp
        #         ls_abaa_temp = psi_a_mv*psi_b_mv*psi_a_mv*psi_a_mv*ls_n[m]
        #         ls_abaa[p] = ls_abaa[p] + ls_abaa_temp
             
        #         ls_abbb_temp = psi_a_mv*psi_b_mv*psi_b_mv*psi_b_mv*ls_n[m]
        #         ls_abbb[p] = ls_abbb[p] + ls_abbb_temp
        #         ls_baaa_temp = psi_b_mv*psi_a_mv*psi_a_mv*psi_a_mv*ls_n[m]
        #         ls_baaa[p] = ls_baaa[p] + ls_baaa_temp
             
        #         ls_aabb_temp = psi_a_mv*psi_a_mv*psi_b_mv*psi_b_mv*ls_n[m]
        #         ls_aabb[p] = ls_aabb[p] + ls_aabb_temp
             
             
        g_abba_dot2 = g_baab_dot2
        g_bbaa = g_aabb
        ls_bbaa = ls_aabb

        #F_k_t = numpy.exp(-1j*(E_0k - ls_a)*tlist - numpy.conjugate(g_a).reshape(-1,1)) ##dimensionless
        #A_k_t = numpy.exp(-1j*(E_0k + ls_a)*tlist - g_a.reshape(-1,1)) #dimensionless
        ##For downhill exciton transfer (b -> a), b>a   N_ab
        #N_kl1_t = (g_baab_dot2 - (g_baaa_dot - g_babb_dot - 2j*ls_babb)*(g_abaa_dot - g_abbb_dot - 2j*ls_abbb)).reshape(-1,1)*numpy.exp(2*(g_aabb.reshape(-1,1) + 1j*ls_aabb.reshape(-1,1)*tlist))
        ##For uphill exciton transfer (a -> b), b>a    N_ba
        #N_kl2_t = (g_abba_dot2 - (g_abbb_dot - g_abaa_dot - 2j*ls_abaa)*(g_babb_dot - g_baaa_dot - 2j*ls_baaa)).reshape(-1,1)*numpy.exp(2*(g_bbaa.reshape(-1,1) + 1j*ls_bbaa.reshape(-1,1)*tlist))

        #for k' < k
        f1 = numpy.zeros((len(pairs[:,1]),Nt),dtype=numpy.float)
        #for k' > k
        f2 = numpy.zeros((len(pairs[:,1]),Nt),dtype=numpy.float)       
             
        RR = numpy.zeros((Na,Na)) #2nd subscript is initial state, 1st subscript is final state - eg K12 means rate for transfer fr exciton 2 to 1; %upper matrix = downhill transfer; lower matrix = uphill
        RR_temp = numpy.zeros((Na,Na))

        #for p in range(len(pairs[:,1])):
        #    f1[p,:] = numpy.conjugate(F_k_t[pairs[p,2],:])*A_k_t[pairs[p,1],:]*N_kl1_t[p,:] #ps-2, downhill transfer
        #    f2[p,:] = numpy.conjugate(F_k_t[pairs(p,1),:])*A_k_t[pairs[p,2],:]*N_kl2_t[p,:] #ps-2, uphill transfer
        #    RR[pairs[p,1],pairs[p,2]] = 2*numpy.real(spi.simps(f1[p,:],tlist[:-2])) #for downhill transfer, upper triangular elements, ps-1
        #    RR[pairs[p,2],pairs[p,1]] = 2*numpy.real(spi.simps(f2[p,:],tlist[:-2])) #for uphill transfer, lower triangular elements, ps-1      
        
        RR_bbaa = -numpy.sum(RR, axis = 0)
        RR = numpy.diag(numpy.diag(RR_bbaa)) + RR  #2nd subscript is initial state, 1st subscript is final state - eg K12 means rate for transfer fr exciton 2 to 1
                

        print("I am called from outside")
        print(m_ket[:,m].reshape(-1,1))
        print(SS[:,1])
       # print(g_baab_dot2)
       # print(g_abaa_dot)
       # print(g_babb_dot)
       # print(g_abbb_dot)
       # print(g_baaa_dot)
       # print(g_aabb)
       # print(ls_babb)
       # print(ls_abaa)
        
       # print(ls_abbb)
       # print(ls_baaa)
       # print(ls_aabb)
    
       # print("E_0k is:")
       # print(E_0k)
    
       # print("pairs are:")
       # print(pairs)
  
       # print(hD)
       # print(SS)
       # print(RR)
    
       # print(m_ket)
       # print(g_a)
       # print(ls_a)
       # print(Na)
       # print(Nc)
       # print(Nt)
    
        qr.stop()
                    
                    
