# -*- coding: utf-8 -*-

"""
*******************************************************************************

      MODIFIED REDFIELD RATE MATRIX

*******************************************************************************
"""  

import numpy
from scipy import integrate
import scipy

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
            self._set_rates2()          
            self._is_initialized = True

    def _set_rates2(self,force_detbalance=True):
        """
        
        """

        Na = self.ham.dim
        Nc = self.sbi.N 
        tt = self.tt
    
        # Eigen problem
        hD,SS = numpy.linalg.eigh(self.ham._data) #hD=eigenvalues, SS=eigenvectors
        
        Nt = self.sbi.CC.timeAxis.length
        dt = self.sbi.CC.timeAxis.step 
        time = self.sbi.CC.timeAxis

        time_full = numpy.zeros(Nt*2)
        time_full[time._length:] = tt.copy()
        time_full[1:time._length] = -tt[:0:-1]
        time_full[0] = time_full[1] - dt


        # fill the corr. matrix 
        Cmat = numpy.zeros((Na,Nt),dtype=numpy.complex128)
        lambdas = numpy.zeros(Na,dtype=numpy.float64)
        for ii in range(Na):
            Cmat[ii,:] = self.sbi.get_coft(ii,ii)
        for ii in range(Na-1):   # assume single excitation per monomer - not correct one should replace it with: aggreg.Nb[1]
            lambdas[ii+1] = self.sbi.get_reorganization_energy(ii)

        # Rate matrix init
        rates = numpy.zeros((Na,Na),dtype=numpy.float64)
        
        for M in range(Na):
            if force_detbalance:
                minL = M+1
            else:
                minL=0
            for L in range(minL,Na):
                if M==L:
                    continue
                # Here we assume noncorrelated sites (diagonal correlation function matrix)
                # prepare individual quantities
                cMMMM = numpy.dot(SS[:,M]**4,Cmat)
                cLLLL = numpy.dot(SS[:,L]**4,Cmat)
                cMLLL = numpy.dot((SS[:,L]**3)*SS[:,M],Cmat)
                cMMML = numpy.dot((SS[:,M]**3)*SS[:,L],Cmat)
                cMMLL = numpy.dot((SS[:,M]**2)*(SS[:,L]**2),Cmat)

                gMMMM = _c2g(time,cMMMM)
                gLLLL = _c2g(time,cLLLL)
                gMMLL = _c2g(time,cMMLL)

                hMLLL = _c2h(time,cMLLL)
                hMMML = _c2h(time,cMMML)

                lambdaLLLL = numpy.dot(SS[:,L]**4,lambdas)
                lambdaMMLL = numpy.dot((SS[:,M]**2)*(SS[:,L]**2),lambdas)
                lambdaMLLL = numpy.dot(SS[:,M]*(SS[:,L]**3),lambdas)

                omML = hD[L]-hD[M]

                # Compute rates (time dependent function)
                om = omML-2*lambdaLLLL+2*lambdaMMLL
                ft_full = numpy.zeros(Nt*2,dtype=numpy.complex128)
                ft = numpy.exp(- gLLLL - gMMMM + 2*gMMLL)
                ft *= (cMMLL - (hMLLL - hMMML + 2*1j*lambdaMLLL)**2 )

                # Replace the time integration over very oscillating function by the fourier transform
                # Prepare the negative part of the function - because we want the real part of the positive
                # part of the function the negative has to be the complex conjugate and we get exactly what
                # we need\
                ft_full[Nt:] = ft.copy()
                ft_full[1:Nt] = numpy.conj(ft[:0:-1])
                # Perform the fourier transform
                fw_full = 2*Nt*numpy.fft.fftshift(numpy.fft.ifft(numpy.fft.fftshift(ft_full)))*dt
                fw_full = numpy.real(fw_full)
                faxis = time.get_FrequencyAxis()
                # Define the result as DFunction to enable the interpolation
                Fw = qr.DFunction(x=faxis,y=fw_full)
                
                rates[M,L] = Fw.at(om)
        
        if force_detbalance:
            T = self.sbi.get_temperature()
            for M in range(Na):
                for L in range(M+1,Na):
                    if M==L:
                        continue
                    
                    om00 = hD[M]-hD[L] + numpy.dot(SS[:,L]**4,lambdas) - numpy.dot(SS[:,M]**4,lambdas)
                    
                    rates[L,M] = rates[M,L] * numpy.exp(om00/kB_intK/T)
            
        # Compute the diagonal elements
        for L in range(Na):
            rates[L,L] = -numpy.sum(rates[:,L])
        
        self.rates = rates
            
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
                    
def _c2g(timeaxis,coft):
    """ Converts correlation function to lineshape function

    Explicit numerical double integration of the correlation
    function to form a lineshape function.

    Parameters
    ----------

    timeaxis : cu.oqs.time.TimeAxis
        TimeAxis of the correlation function

    coft : complex numpy array
        Values of correlation function given at points specified
        in the TimeAxis object


    """

    ta = timeaxis
    rr = numpy.real(coft)
    ri = numpy.imag(coft)
    sr = scipy.interpolate.UnivariateSpline(ta.data,
                        rr,s=0).antiderivative()(ta.data)
    sr = scipy.interpolate.UnivariateSpline(ta.data,
                        sr,s=0).antiderivative()(ta.data)
    si = scipy.interpolate.UnivariateSpline(ta.data,
                        ri,s=0).antiderivative()(ta.data)
    si = scipy.interpolate.UnivariateSpline(ta.data,
                        si,s=0).antiderivative()(ta.data)
    gt = sr + 1j*si
    return gt

def _c2h(timeaxis,coft):
    """ Converts correlation function to derivative of lineshape function

    Explicit numerical single integration of the correlation
    function to form a time redivative of lineshape function.

    Parameters
    ----------

    timeaxis : cu.oqs.time.TimeAxis
        TimeAxis of the correlation function

    coft : complex numpy array
        Values of correlation function given at points specified
        in the TimeAxis object


    """

    ta = timeaxis
    rr = numpy.real(coft)
    ri = numpy.imag(coft)
    sr = scipy.interpolate.UnivariateSpline(ta.data,
                        rr,s=0).antiderivative()(ta.data)
    si = scipy.interpolate.UnivariateSpline(ta.data,
                        ri,s=0).antiderivative()(ta.data)
    ht = sr + 1j*si
    return ht
