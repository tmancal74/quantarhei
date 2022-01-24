# -*- coding: utf-8 -*-
import numpy
import scipy.interpolate as interp

from ..hilbertspace.hamiltonian import Hamiltonian
from ..liouvillespace.systembathinteraction import SystemBathInteraction
from .foerstertensor import FoersterRelaxationTensor
from ..corfunctions.correlationfunctions import c2g
from ...core.managers import energy_units

from ...core.time import TimeDependent

class TDFoersterRelaxationTensor(FoersterRelaxationTensor, TimeDependent):
    """Weak resonance coupling relaxation tensor by Foerster theory
    
    
    
    """
    def __init__(self, ham, sbi, initialize=True, cutoff_time=None):
        
        super().__init__(ham, sbi, initialize, cutoff_time)

        
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
                        
            KK = _td_reference_implementation(Na, Nt, HH, tt, gt, ll)
   
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


        ret = 2.0*numpy.real(hoft)
        
        return ret

def _td_reference_implementation(Na, Nt, HH, tt, gt, ll):
                                         
        #
        # Rates between states a and b
        # 
        KK = numpy.zeros((Nt,Na,Na), dtype=numpy.float64)
        for a in range(Na):
            for b in range(Na):
                if a != b:
                    ed = HH[b,b] # donor
                    ea = HH[a,a] # acceptor
                    KK[:,a,b] = (HH[a,b]**2)*_td_fintegral(tt, gt[a,:],
                                                           gt[b,:],
                                                           ed, ea, ll[b])
   
        return KK        