# -*- coding: utf-8 -*-
import numpy
import scipy.interpolate as interp

from ..hilbertspace.hamiltonian import Hamiltonian
from ..liouvillespace.systembathinteraction import SystemBathInteraction
from .relaxationtensor import RelaxationTensor
from ..corfunctions.correlationfunctions import c2g
from ...core.managers import energy_units

class FoersterRelaxationTensor(RelaxationTensor):
    """Weak resonance coupling relaxation tensor by Foerster theory
    
    
    
    """
    def __init__(self, ham, sbi, initialize=True, cutoff_time=None):
        
        super().__init__()
        
        if not isinstance(ham, Hamiltonian):
            raise Exception("First argument must be a Hamiltonian")
            
        if not isinstance(sbi, SystemBathInteraction):
            raise Exception("Second argument must be of"
                           +" type SystemBathInteraction")
            
        self._is_initialized = False
        self._has_cutoff_time = False
        
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

    def initialize(self):
        
        #
        # Tensor data
        #
        Na = self.dim
        self.data = numpy.zeros((Na,Na,Na,Na),dtype=numpy.complex128)
      
        
        with energy_units("int"):
            
            # Hamiltonian matrix
            HH = self.Hamiltonian.data
            tt = self.SystemBathInteraction.TimeAxis.data
            sbi = self.SystemBathInteraction 
            
            # line shape functions
            gt = numpy.zeros((Na, sbi.TimeAxis.length),
                             dtype=numpy.complex64)
    
            # SBI is defined with "sites"
            for ii in range(1, Na):
                gt[ii,:] = c2g(sbi.TimeAxis, sbi.CC.get_coft(ii-1,ii-1))
            
            # reorganization energies
            ll = numpy.zeros(Na)
            for ii in range(1, Na):
                ll[ii] = sbi.CC.get_reorganization_energy(ii-1,ii-1)
                        
            KK = _reference_implementation(Na, HH, tt, gt, ll)
   
            #
            # Transfer rates
            #                                                          
            for aa in range(self.dim):
                for bb in range(self.dim):
                    if aa != bb:
                        self.data[aa,aa,bb,bb] = KK[aa,bb]
                
            #  
            # calculate dephasing rates and depopulation rates
            #
            self.updateStructure()
        
        


def _reference_implementation(Na, HH, tt, gt, ll):
    """Reference implementation of Foerster rates 
    
    Calculate the rates between specified sites using standard Foerster
    theory. 
    
    Reference:  
    L. Valkunas, D. Abramavicius, and T. Mančal, Molecular Excitation
    Dynamics and Relaxation, Wiley-VCH, Berlin (2013), page: 

    Parameters
    ----------
    
    Na : integer
        Number of sites in the problem (rank of the rate matrix)
        
    HH : float array
        Hamiltonian matrix
        
    tt : float array
        Time points in which the line shape functions are given
        
    gt : complex array
        Line shape functions values at give time points.
        First index corresponds to the site, the second to the time point
        
    ll : array
        Reorganization energies on sites
        
    Returns
    -------
    
    KK : float array
        Rate matrix with zeros on the diagonal

    """
                                     
    #
    # Rates between states a and b
    # 
    KK = numpy.zeros((Na,Na), dtype=numpy.float64)
    for a in range(Na):
        for b in range(Na):
            if a != b:
                ed = HH[b,b] # donor
                ea = HH[a,a] # acceptor
                KK[a,b] = (HH[a,b]**2)*_fintegral(tt, gt[a,:], gt[b,:],
                                                  ed, ea, ll[b])
   
    return KK

    
    
def _fintegral(tt, gtd, gta, ed, ea, ld):
    """Foerster integral
    
    
    Parameters
    ----------
    tt : numpy array
        Time 
        
    gtd : numpy array
        lineshape function of the donor transition

    gta : numpy array
        lineshape function of the acceptor transition 
        
    ed : float
        Energy of the donor transition
        
    ea : float
        Energy of the acceptor transition

    ld : float
        Reorganization energy of the donor             

    Returns
    -------
    ret : float
        The value of the Foerster integral            
    
    """
    #fl = numpy.exp(-gtd +1j*(ed-2.0*ld)*tm.data)
    #ab = numpy.exp(-gta -1j*ea*tm.data)
    #prod = ab*fl

    prod = numpy.exp(-gtd-gta +1j*((ed-ea)-2.0*ld)*tt)
    
    preal = numpy.real(prod)
    pimag = numpy.imag(prod)
    splr = interp.UnivariateSpline(tt,
                               preal, s=0).antiderivative()(tt)
    spli = interp.UnivariateSpline(tt,
                               pimag, s=0).antiderivative()(tt)
    hoft = splr + 1j*spli


    ret = 2.0*numpy.real(hoft[len(tt)-1])
    
    return ret
    