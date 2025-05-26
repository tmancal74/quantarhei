# -*- coding: utf-8 -*-

"""
*******************************************************************************

      MODIFIED REDFIELD RATE MATRIX

*******************************************************************************
"""  
import time

import numpy
from scipy import integrate
import scipy

from scipy.integrate import simpson as simps 

#from quantarhei.core.implementations import implementation

from ...hilbertspace.hamiltonian import Hamiltonian
from ...liouvillespace.systembathinteraction import SystemBathInteraction

from .... import REAL
from .... import COMPLEX

from ....core.units import convert

import matplotlib.pyplot as plt


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
        
        #c1 = self.sbi.CC.get_coft(0,0)
        
        #print("Lambda:")
        #print(lam4[2,2,2,2])

        
        #plt.plot(tt,c1,"-k")
        #plt.plot(tt,c4[1,1,1,1,:],"--b")
        #plt.show()
        
        #quit()


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
    
        Nt : integer
        
        
        hD : real
        Eigenenergies of the Hamiltonian
        
        
    
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
        F_k_t = numpy.zeros((Na,Nt),dtype=COMPLEX)    
        A_k_t = numpy.zeros((Na,Nt),dtype=COMPLEX) 
        N_kl_t = numpy.zeros((Na,Na,Nt),dtype=COMPLEX) 
        f = numpy.zeros((Na,Na,Nt),dtype=COMPLEX) 
        RR = numpy.zeros((Na,Na),dtype=COMPLEX)  
        RR1 = numpy.zeros((Na,Na),dtype=COMPLEX)
        
        
        #lam1 = (convert(numpy.imag(-h4[2,1,1,2,Nt-1]),"int","1/cm"))
        #lam2 = convert(lam4[2,1,1,2],"int","1/cm")
        #fac = lam2/lam1
        #print(fac)
        #lam1 = (convert(numpy.imag(-h4[1,2,2,1,Nt-1]),"int","1/cm"))
        #lam2 = convert(lam4[1,2,2,1],"int","1/cm")
        #fac = lam2/lam1
        #print(fac)        
        #quit()
        
        #h4 = fac*h4
        #g4 = fac*g4
        #c4 = fac*c4
        
        
        for ii in range(Na):
            E_0k[ii] = hD[ii] - lam4[ii,ii,ii,ii]
            
        for a in range(Na):
            F_k_t[a,:] = (numpy.exp(-1j*(E_0k[a] - lam4[a,a,a,a])*tt[:]
                                    - g4[a,a,a,a,:]) )
                                    # - numpy.conjugate(g4[a,a,a,a,:])) )
            A_k_t[a,:] = (numpy.exp(-1j*(E_0k[a] + lam4[a,a,a,a])*tt[:]
                                     - g4[a,a,a,a,:]) )
              
        for a in range(Na):
            for b in range(Na):
                N_kl_t[a,b,:] = ((c4[b,a,a,b,:] 
                    - (h4[b,a,b,b,:] - h4[b,a,a,a,:] + 2j*lam4[b,a,b,b])
                     *(h4[b,b,a,b,:] - h4[a,a,a,b,:] + 2j*lam4[b,b,a,b]))
                     *numpy.exp(2*(g4[a,a,b,b,:] + 2j*lam4[a,a,b,b]*tt[:])))
 
        """
        for a in range(Na):
            for b in range(Na):
                
                #f[a,b,:] = numpy.exp(1j*(hD[b]-hD[a])*tt[:])*c4[b,a,a,b,:]
                
                
                f[b,a,:] = (numpy.exp(1j*(E_0k[a]-E_0k[b])*tt[:])
                *(c4[a,b,b,a,:]-
                (h4[b,b,b,a,:]-h4[a,a,b,a,:]+2.0*1j*lam4[b,a,a,a])
                *(h4[a,b,a,a,:]-h4[a,b,b,b,:]+2.0*1j*lam4[a,b,a,a]))
                *numpy.exp(-g4[a,a,a,a,:]-g4[b,b,b,b,:]
                          +g4[a,a,b,b,:]+g4[b,b,a,a,:]
                          +2.0*1j*(lam4[a,a,b,b]-lam4[a,a,a,a])*tt[:]))
                
        """
                
        f1 = numpy.zeros((Na,Na,Nt),dtype=COMPLEX)                                 
        for a in range(Na):
            for b in range(Na):
                f1[a,b,:] = (numpy.conjugate(F_k_t[b,:])
                            *A_k_t[a,:]*N_kl_t[a,b,:])

                RR[a,b] = 2.0*numpy.real(simps(f1[a,b,:],tt))
 
                  
        #print(RR-RR1)
                
        #for a in range(Na):
        #    RR[a,a] = 0
        #RR_bbaa = -numpy.sum(RR, axis = 0)
        #RR = numpy.diag(RR_bbaa) + RR  


        return RR

        #qr.stop()
                    
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
                 
def ssModifiedRRM(Na, Nt, en, ln, SS, tt, ct, ht, gt, method="ModifiedRedfield"):
    """ Modifield Redfield rates
    
    Na : integer
        Number of sites
        
    Nt : integer
        Number of time steps
        
    en : real
        Eigenenergies 
        
    SS : transformation matrix
    
    tt : time
    
    ct : complex vector
    
    ht : complex vector
    
    gt : conmplex vector
    
    
    """
        
    RR = numpy.zeros((Na, Na), dtype=REAL)   
    f = numpy.zeros((Na, Na, Nt), dtype=COMPLEX)
    Ab = numpy.zeros(Nt, dtype=COMPLEX)
    Fl = numpy.zeros(Nt, dtype=COMPLEX)
    Nn = numpy.zeros(Nt, dtype=COMPLEX)
    
    if method == "ModifiedRedfield":
        
       
        for a in range(Na):
            for b in range(Na):
                Ab[:] = 0.0
                Fl[:] = 0.0
                Nn[:] = 0.0
                f[a,b,:] = Ab[:]*Fl[:]*Nn[:]
                RR[a,b] = 2.0*numpy.real(simps(f[a,b,:],tt))
        
    
    elif method == "Redfield":
    
        RR[:,:] = 0.0
    
    return RR
                    
