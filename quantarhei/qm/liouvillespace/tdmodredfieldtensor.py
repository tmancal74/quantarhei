# -*- coding: utf-8 -*-
import time
import numpy
#import scipy.interpolate as interp

from ..hilbertspace.hamiltonian import Hamiltonian
from ..liouvillespace.systembathinteraction import SystemBathInteraction
from .relaxationtensor import RelaxationTensor
#from .rates.foersterrates import FoersterRateMatrix
from ...core.time import TimeDependent
from ...core.managers import energy_units, eigenbasis_of
from ... import COMPLEX, REAL

class TDModRedfieldRelaxationTensor(RelaxationTensor, TimeDependent):
    """Time-dependent Modifield Redfield Tensor
    
    
    """
    def __init__(self, ham, sbi, initialize=True, cutoff_time=None, 
                 theory="mod_eq"):
        
        self._initialize_basis()
        
        if not isinstance(ham, Hamiltonian):
            raise Exception("First argument must be a Hamiltonian")
            
        if not isinstance(sbi, SystemBathInteraction):
            raise Exception("Second argument must be of"
                           +" type SystemBathInteraction")
            
        self._is_initialized = False
        self._has_cutoff_time = False
        self.as_operators = False
        
        self.theory = theory
        
        if cutoff_time is not None:
            self.cutoff_time = cutoff_time
            self._has_cutoff_time = True            
            
        self.Hamiltonian = ham
        self.dim = ham.dim
        self.SystemBathInteraction = sbi
        self.TimeAxis = sbi.TimeAxis
    
        
        if initialize:
            self.initialize()
            self._data_initialized = True 
            self._is_initialized = True
            
        else:
            self._data_initialized = False
        
        
        self.has_Iterm = False
        
        self.is_time_dependent = True
        

    def initialize(self):
        
        tdep = True
        
        #
        # get elementary properties
        #
        HH = self.Hamiltonian
        sbi = self.SystemBathInteraction
        timea = self.TimeAxis
        if HH.Nblocks == 2:
            Na = HH._data.shape[0]-1
        else:
            raise Exception("Theory is not implemented beyond single exciton block")

        ll = numpy.zeros(Na+1,dtype=REAL)
        cfce = sbi.CC
        
        cfce.create_one_integral()
        cfce.create_double_integral() #g(t)

        #
        # Tensor data
        #
        if tdep:
            self.data = numpy.zeros((timea.length, Na+1,Na+1,Na+1,Na+1),
                                    dtype=COMPLEX)
        else:
            self.data = numpy.zeros((Na+1,Na+1,Na+1,Na+1),dtype=COMPLEX)
        
        
        with energy_units("int"):

            for ii in range(Na):
                ll[ii+1] = cfce.get_reorganization_energy(ii,ii)
            ee, ss = numpy.linalg.eigh(HH._data)   

            #t1 = time.time()
            RR, Iterm = ssmodr(Na, ee, ll, ss, timea.data, cfce._cofts,
                        cfce._hofts, cfce._gofts,
                        cfce.cpointer, tdep=tdep, theory=self.theory)
            self.RR = RR
            #t2 = time.time()
            #print("Redfield: '"+self.theory+"' calculated in", t2-t1, "sec")
            
    
            #
            # Transfer rates
            #                     
            with eigenbasis_of(HH):                                     
                for aa in range(self.dim):
                    for bb in range(self.dim):
                        self.Iterm[:,aa,bb] = Iterm[aa,bb,:]
                        if aa != bb:
                            self.data[:,aa,aa,bb,bb] = RR[aa,bb,:]
                
                #  
                # calculate dephasing rates and depopulation rates
                #
                #self.updateStructure()
                # depopulation rates 
                for nn in range(self.dim):
                    self.data[:,nn,nn,nn,nn] -= (numpy.trace(self.data[:,:,:,nn,nn],
                                                              axis1=1,axis2=2)
                                                - self.data[:,nn,nn,nn,nn])
                    
                # dephasing rates 
                for nn in range(self.dim):    
                    for mm in range(nn+1,self.dim):
                        self.data[:,nn,mm,nn,mm] = (self.data[:,nn,nn,nn,nn]
                                                  +self.data[:,mm,mm,mm,mm])/2.0
                        self.data[:,mm,nn,mm,nn] = self.data[:,nn,mm,nn,mm] 
                        



def intfull(F, dt, omega):
    """Discrete integration of an oscillating function
    
    """

    Nt = F.shape[0]
    x = (1.0 - numpy.exp(-1j*omega*dt))/(1j*omega*dt)
    y = 1j*(1.0 - x)/omega
    I = 2.0*numpy.real(y)*numpy.sum(F)-numpy.conj(y)*F[Nt-1] - y*F[0]

    return I


def intrun(F, dt, omega):
    """Running intergration of an oscillating function of one parameter"""

    Nt = F.shape[0]
    I = numpy.zeros(Nt, dtype=COMPLEX)
    x = (1.0 - numpy.exp(-1j*omega*dt))/(1j*omega*dt)
    y = (1.0 - x)/(1j*omega)   

    for ii in range(1,Nt):
        I[ii] = I[ii-1] + y*F[ii] + numpy.conj(y)*F[ii-1]

    print(dt, omega, numpy.conj(y), y)

    return I

def intruntrap(F, dt):
    Nt = F.shape[0]
    I = numpy.zeros(Nt, dtype=COMPLEX)
    for ii in range(1,Nt):
        I[ii] = I[ii-1]+ F[ii]*dt/2.0 + F[ii-1]*dt/2.0

    return I



def intfullsec(F, dt, omega, Nt):
    """Integrate only a section of a give function F
    
    This is used when the function depends explicitely on time (the upper limit)

    """
    x = (1.0 - numpy.exp(-1j*omega*dt))/(1j*omega*dt)
    y = 1j*(1.0 - x)/omega
    if Nt > 0:
        I = 2.0*numpy.real(y)*numpy.sum(F[:Nt])-numpy.conj(y)*F[Nt-1] - y*F[0]
    else:
        I = 0.0
    return I




def ssmodr(Na, ee, ll, ss, tt, cofts, hofts, gofts, pntr, tdep=False, theory="std"):
    """Various version of the Redfield relaxation rate calculation

        1) Standard Redfield theory in time-dependent and time-independent versions
        2) Modified Redfield theory in time-dependent and time-independent version
        3) Non-equilibrium Modified Redfield theory in time-dependent version


    Parameters
    ----------

    Na : int
        Number of molecular sites. Na + 1 gives the dimension of the tensor, because 
        the ground state is also included

    ee : float array
        Array of Hamiltonian eigenenergies

    ll : float array
        Array of site reorganization energies

    ss : float array 2D
        Transformation matrix composed of Hamiltonian eigenvectors

    tt : float array
        Time axis. Defines time step as dt = tt[1] - tt[0] and the number of 
        time steps as Nt = tt.shape[0]
        
    cofts : float array
        Array of correlation functions

    hofts : float array
        Array of first integrals of the correlation functions

    gofts : float array
        Array of second integrals of the correlation functions, line-shape functions

    pntrs : integer array 2D
        Array of mapping indices between site pairs and the correlation function matrices

    tdep : bool
        True of the output will be time-dependent rate matrix

    theory : string
        Which theory is to be used for the calculation
        "std" - standard Redfield theory
        "mod_eq" - equilibrium Modified Redfield theory
        "mod_noneq" - non-equilibrium Modified Redfield theory
        "std_test" - testing mode, one should get the same results 
                     as for "std" with the option
        

    
    """

    #ee, ss = numpy.linalg.eigh(hh)

    # number of time steps
    Nt = tt.shape[0]

    # Relaxation rate matrix (may be time-dependent)
    if tdep:
        RR = numpy.zeros((Na+1,Na+1,Nt), dtype=REAL)
        Iterm = numpy.zeros((Na+1,Na+1,Nt), dtype=COMPLEX)
    else:
        RR = numpy.zeros((Na+1,Na+1), dtype=REAL)
        Iterm = numpy.zeros((Na+1,Na+1), dtype=REAL)
    
    # Various storage variables
    cbaab = numpy.zeros(Nt, dtype=COMPLEX)
    Nab = numpy.zeros(Nt, dtype=COMPLEX)
    
    # time step
    dt = tt[1]-tt[0]

    if theory == "std":
        #
        # Standard Redfield theory 
        # [FINISHED]
        #
        for a in range(Na+1):

            aa = -1j*ee[a]*tt
            Aa = numpy.exp(aa)

            for b in range(Na+1):
                
                if a != b:

                    fb = -1j*ee[b]*tt
                    Fb = numpy.exp(fb)    

                    cbaab[:] = 0.0
                    for n in range(Na):
                        cbaab += (ss[b,n]**2)*(ss[a,n]**2)*cofts[pntr[n,n],:]

                    Nab = cbaab
                    
                    kab = numpy.conj(Fb)*Nab*Aa

                    om = ee[b]-ee[a]
                    if tdep:
                        RR[a,b,:] = 2.0*numpy.real(intrun(kab, dt, om))
                    else:
                        RR[a,b] = 2.0*numpy.real(intfull(kab, dt, om))


    elif theory == "std_test":
        #
        # Standard Redfield by integrating each time from scratch
        # This is only to test intfullsec() routine and make shure the 
        # [FINISHED]
        #
        for a in range(Na+1):
            
            aa = -1j*ee[a]*tt 
            Aa = numpy.exp(aa)

            for b in range(Na+1):

                if a != b:

                    cbaab[:] = 0.0
                    
                    for n in range(Na):
                        cbaab += (ss[b,n]**2)*(ss[a,n]**2)*cofts[pntr[n,n],:]

                    Nab[:] = 0.0
                    for it in range(Nt):

                        tm = it + 1
                        fb = -1j*ee[b]*tt[:tm]
                        Fb = numpy.exp(fb)

                        Nab[:tm] = cbaab[:tm] #
                    
                        kab = numpy.conj(Fb)*Nab[:(tm)]*Aa[:(tm)]  

                        om = ee[b]-ee[a]
                        if tdep:
                            RR[a,b,it] = 2.0*numpy.real(intfullsec(kab, dt, om, it+1))
                        else:
                            RR[a,b] = 2.0*numpy.real(intfullsec(kab, dt, om, it+1))                     


    elif theory == "mod_eq":
        #
        # Modified Redfield theory (equilibrium version)
        # 
        #
        mm = numpy.zeros(Nt, dtype=COMPLEX)

        for a in range(Na+1):

            aa = -1j*ee[a]*tt
            for n in range(Na):
                aa -= (ss[a,n]*ss[a,n]*ss[a,n]*ss[a,n])*gofts[pntr[n,n],:]
            Aa = numpy.exp(aa)

            for b in range(Na+1):

                if a != b:

                    fb = -1j*ee[b]*tt
                    for n in range(Na):
                        fb -= (ss[b,n]*ss[b,n]*ss[b,n]*ss[b,n]) \
                        *(numpy.conj(gofts[pntr[n,n],:])-2.0*1j*ll[n]*tt)
                    Fb = numpy.exp(fb)    

                    cbaab[:] = 0.0
                    for n in range(Na):
                        cbaab += (ss[b,n]**2)*(ss[a,n]**2)*cofts[pntr[n,n],:]

                    Nab[:] = 0.0  # we use Nab storage for nab quantity
                    for n in range(Na):
                        Nab += 2.0*(ss[b,n]**2)*(ss[a,n]**2)*(gofts[pntr[n,n],:]+ 1j*ll[n]*tt)
                    
                    mm[:] = 0.0
                    for n in range(Na):
                        mm += ss[a,n]*ss[b,n]*((ss[b,n]**2) - (ss[a,n]**2))*hofts[pntr[n,n],:] \
                               + 2.0*1j*ss[a,n]*ss[b,n]*(ss[b,n]**2)*ll[n]

                    Nab = numpy.exp(Nab)*(cbaab - mm**2)
                    
                    kab = numpy.conj(Fb)*Nab*Aa

                    om = ee[b]-ee[a]
                    if tdep:
                        RR[a,b,:] = 2.0*numpy.real(intrun(kab, dt, om))
                    else:
                        RR[a,b] = 2.0*numpy.real(intfull(kab, dt, om))


    elif theory == "mod_eq_test":
        #
        # Modified Redfield theory (equilibrium version)
        # Testing integration from scratch
        #
        mm = numpy.zeros(Nt, dtype=COMPLEX)
        fb = numpy.zeros(Nt, dtype=COMPLEX)

        for a in range(Na+1):

            aa = -1j*ee[a]*tt
            for n in range(Na):
                aa -= (ss[a,n]*ss[a,n]*ss[a,n]*ss[a,n])*gofts[pntr[n,n],:]
            Aa = numpy.exp(aa)

            for b in range(Na+1):

                if a != b:  

                    cbaab[:] = 0.0
                    for n in range(Na):
                        cbaab += (ss[b,n]**2)*(ss[a,n]**2)*cofts[pntr[n,n],:]

                    for it in range(Nt):

                        tm = it + 1

                        fb[:tm] = -1j*ee[b]*tt[:tm]
                        for n in range(Na):
                            fb[:tm] -= (ss[b,n]*ss[b,n]*ss[b,n]*ss[b,n]) \
                            *(numpy.conj(gofts[pntr[n,n],:tm])-2.0*1j*ll[n]*tt[:tm])
                        Fb = numpy.exp(fb) 


                        Nab[:] = 0.0  # we use Nab storage for nab quantity
                        for n in range(Na):
                            Nab[:tm] += 2.0*(ss[b,n]**2)*(ss[a,n]**2)*(gofts[pntr[n,n],:tm]+ 1j*ll[n]*tt[:tm])
                        
                        mm[:] = 0.0
                        for n in range(Na):
                            # the complex conjugation fixes the result of Seibt et al. 2017
                            mm[:tm] += ss[a,n]*ss[b,n]*((ss[b,n]**2)*hofts[pntr[n,n],:tm] 
                                                - (ss[a,n]**2)*(hofts[pntr[n,n],:tm])) \
                                    + 2.0*1j*ss[a,n]*ss[b,n]*(ss[b,n]**2)*ll[n]

                        Nab[:tm] = numpy.exp(Nab[:tm])*(cbaab[:tm] - (mm[:tm]**2))

                        kab = numpy.conj(Fb)*Nab*Aa

                        om = ee[b]-ee[a]
                        if tdep:
                            RR[a,b,it] = 2.0*numpy.real(intfullsec(kab, dt, om, tm))
                        else:
                            RR[a,b] = 2.0*numpy.real(intfullsec(kab, dt, om, tm))


    elif theory == "mod_noneq":
        #
        # Modified Redfield theory (non-equilibrium version)
        # 
        #
        mm = numpy.zeros(Nt, dtype=COMPLEX)
        mr = numpy.zeros(Nt, dtype=COMPLEX)
        fb = numpy.zeros(Nt, dtype=COMPLEX)

        for a in range(Na+1):

            aa = -1j*ee[a]*tt
            for n in range(Na):
                aa -= (ss[a,n]*ss[a,n]*ss[a,n]*ss[a,n])*gofts[pntr[n,n],:]
            Aa = numpy.exp(aa)

            for b in range(Na+1):

                if a != b:  

                    cbaab[:] = 0.0
                    for n in range(Na):
                        cbaab += (ss[b,n]**2)*(ss[a,n]**2)*cofts[pntr[n,n],:]

                    # here we integrate over time
                    for it in range(Nt):

                        tm = it + 1

                        Nab[:] = 0.0  # we use Nab storage for nab quantity
                        mm[:] = 0.0
                        mr[:] = 0.0

                        fb[:tm] = -1j*ee[b]*tt[:tm]
                        for n in range(Na):

                            hrev = numpy.flipud(hofts[pntr[n,n],:tm])
                            grev = numpy.flipud(gofts[pntr[n,n],:tm])
                            reorg_term = numpy.imag(gofts[pntr[n,n],it] - grev)

                            fb[:tm] -= (ss[b,n]*ss[b,n]*ss[b,n]*ss[b,n]) \
                                *(numpy.conj(gofts[pntr[n,n],:tm])
                                +2.0*1j*reorg_term)
                         
                            Nab[:tm] += 2.0*(ss[b,n]**2)*(ss[a,n]**2)*(gofts[pntr[n,n],:tm] 
                                                                           - 1j*reorg_term)
                        
                            mm[:tm] += ss[a,n]*ss[b,n]*((ss[b,n]**2) - (ss[a,n]**2))*hofts[pntr[n,n],:tm] \
                                    - 2.0*1j*ss[a,n]*ss[b,n]*(ss[b,n]**2)*numpy.imag(hofts[pntr[n,n],it])
                                                                                                                                             
                            mr[:tm] += ss[a,n]*ss[b,n]*((ss[b,n]**2) - (ss[a,n]**2))*hofts[pntr[n,n],:tm] \
                                    - 2.0*1j*ss[a,n]*ss[b,n]*(ss[b,n]**2)*numpy.imag(hrev)                                                                                     

                        Fb = numpy.exp(fb)
                        Nab[:tm] = numpy.exp(Nab[:tm])*(cbaab[:tm] - mm[:tm]*mr[:tm])
                        
                        kab = numpy.conj(Fb)*Nab*Aa

                        om = ee[b]-ee[a]
                        if tdep:
                            RR[a,b,it] = 2.0*numpy.real(intfullsec(kab, dt, om, tm))
                        else:
                            RR[a,b] = 2.0*numpy.real(intfullsec(kab, dt, om, tm))


    else:

        raise Exception("Theory type not implemented")

    #
    # diagonal elements
    #
    if tdep:
        for a in range(Na+1):
            for b in range(Na+1):
                if b != a:
                    RR[a,a,:] -= RR[b,a,:]
    else:
        for a in range(Na+1):
            for b in range(Na+1):
                if b != a:
                    RR[a,a] -= RR[b,a]

    # FIXME: where is the factor of 2.0 comming from?
    RR = 2.0*RR
    
    return RR, Iterm

        

