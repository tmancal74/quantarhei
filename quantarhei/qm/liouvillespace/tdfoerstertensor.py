# -*- coding: utf-8 -*-
import numpy
#import matplotlib.pyplot as plt
import scipy.interpolate as interp

#from ..hilbertspace.hamiltonian import Hamiltonian
#from ..liouvillespace.systembathinteraction import SystemBathInteraction
from .foerstertensor import FoersterRelaxationTensor
from ..corfunctions.correlationfunctions import c2g
from ...core.managers import energy_units

from ...core.time import TimeDependent

class TDFoersterRelaxationTensor(FoersterRelaxationTensor, TimeDependent):
    """Weak resonance coupling relaxation tensor by Foerster theory
    
    Time-dependent version
    
    
    
    """
    def __init__(self, ham, sbi, initialize=True, cutoff_time=None):
        
        super().__init__(ham, sbi, initialize, cutoff_time)
        
        self.is_time_dependent = True
        
    def initialize(self):
        
        tt = self.SystemBathInteraction.TimeAxis.data
        Nt = len(tt)
        #
        # Tensor data
        #
        Na = self.dim
        self.data = numpy.zeros((Nt,Na,Na,Na,Na),dtype=numpy.complex128)
        
        
        with energy_units("int"):

            # Hamiltonian matrix
            HH = self.Hamiltonian.data

            sbi = self.SystemBathInteraction
            Na = self.dim
            
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
                        
            KK = self.td_reference_implementation(Na, Nt, HH, tt, gt, ll)
   
            #
            # Transfer rates
            #                                                          
            for aa in range(self.dim):
                for bb in range(self.dim):
                    if aa != bb:
                        self.data[:,aa,aa,bb,bb] = KK[:,aa,bb]
                
            #  
            # calculate dephasing rates and depopulation rates
            #
            self.updateStructure() 
            
            # additional pure dephasing 
            self.add_dephasing()
            
           
    def add_dephasing(self):
       

        # line shape function derivatives
        sbi = self.SystemBathInteraction
        Na = self.dim
        
        ht = numpy.zeros((Na, sbi.TimeAxis.length),
                         dtype=numpy.complex64)
        
        sbi.CC.create_one_integral()
    
        for ii in range(1, Na):
           ht[ii,:] = sbi.CC.get_hoft(ii-1,ii-1)
        
                        
        for aa in range(self.dim):
            for bb in range(self.dim):
                if aa != bb:
 
                    self.data[:,aa,bb,aa,bb] -= (ht[aa,:]+ht[bb,:])


        
    def td_reference_implementation(self, Na, Nt, HH, tt, gt, ll):
        return _td_reference_implementation(Na, Nt, HH, tt,
                                            gt, ll, _td_fintegral)
                                             


def _td_reference_implementation(Na, Nt, HH, tt, gt, ll, fce):
                                     
    #
    # Rates between states a and b
    # 
    KK = numpy.zeros((Nt,Na,Na), dtype=numpy.float64)
    for a in range(Na):
        for b in range(Na):
            if a != b:
                ed = HH[b,b] # donor
                ea = HH[a,a] # acceptor
                KK[:,a,b] = (HH[a,b]**2)*fce(tt, gt[a,:], gt[b,:],
                                             ed, ea, ll[b])
   
    return KK     
            
    
def _td_fintegral(tt, gtd, gta, ed, ea, ld):
        """Time dependent Foerster integral
        
        
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

        prod = numpy.exp(-gtd-gta +1j*((ed-ea)-2.0*ld)*tt)
        
        preal = numpy.real(prod)
        pimag = numpy.imag(prod)
        splr = interp.UnivariateSpline(tt,
                                   preal, s=0).antiderivative()(tt)
        spli = interp.UnivariateSpline(tt,
                                   pimag, s=0).antiderivative()(tt)
        hoft = splr + 1j*spli


        ret = 2.0*numpy.real(hoft)
        
        return ret
       