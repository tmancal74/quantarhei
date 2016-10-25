# -*- coding: utf-8 -*-
import numpy

from ..hilbertspace.hamiltonian import Hamiltonian
from ..liouvillespace.systembathinteraction import SystemBathInteraction
from .relaxationtensor import RelaxationTensor
from ..corfunctions.correlationfunctions import c2g

class FoersterRelaxationTensor(RelaxationTensor):
    """Weak resonance coupling relaxation tensor by Foerster theory
    
    """
    def __init__(self, ham, sbi, initialize=True, cutoff_time=None):
        
        super().__init__()
        
        if not isinstance(ham, Hamiltonian):
            raise Exception("First argument must be a Hamiltonian")
            
        if not isinstance(sbi, SystemBathInteraction):
            raise Exception("Second argument must be of type SystemBathInteraction")
            
        self._is_initialized = False
        self._has_cutoff_time = False
        
        if cutoff_time is not None:
            self.cutoff_time = cutoff_time
            self._has_cutoff_time = True            
            
        self.Hamiltonian = ham
        self.dim = ham.dim
        self.BilinearSystemBathInteraction = sbi
        
        if initialize:
            self.__reference_implementation()
            self._is_initialized = True

    def set_system_bath_interaction(self,sbi):
        self.BilinearSystemBathInteraction = sbi
        

    def _fintegral(self,tm,gtd,gta,ed,ea,ld):
        """Foerster integral
        
        
        Parameters
        ----------
        tm : cu.oqs.time.TimeAxis
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
        fl = numpy.exp(-numpy.conj(gtd) +1j*(ed-2.0*ld)*tm.data)
        ab = numpy.exp(-numpy.conj(gta) -1j*ea*tm.data)
        prod = ab*fl
        ret = 2.0*numpy.real(numpy.sum(prod)*tm.step)
        return ret

    def __reference_implementation(self):
        Na = self.Hamiltonian.dim
        ham = self.Hamiltonian
        HH = self.Hamiltonian.data
        if ham._has_remainder_coupling: 
            JR = self.Hamiltonian.JR
        else:
            JR = numpy.zeros((Na,Na), dtype=numpy.float64)
            for i in range(Na):
                for j in range(i,Na):
                    JR[i,j] = HH[i,j]
                    JR[j,i] = HH[j,i]
                    
        sbi = self.BilinearSystemBathInteraction
        ta = sbi.TimeAxis

        self.data = numpy.zeros((Na,Na,Na,Na))

        # line shape functions
        Gt = numpy.zeros((Na,ta.length),dtype=numpy.complex64)
        gt = numpy.zeros((Na,ta.length),dtype=numpy.complex64)
        for ii in range(1,Na):
            Gt[ii,:] = c2g(ta,
                           sbi.CC.get_coft(ii-1,ii-1))
        try:
            SS = ham.SS
        except:
            SS = numpy.diag(numpy.ones(Na,dtype=numpy.float64))
            
        for aa in range(Na):
            for bb in range(Na):
                """ Here we assume no correlation between sites """
                gt[aa,:] += (SS[bb,aa]**4)*Gt[bb,:]  

        # reorganization energies
        ll = numpy.zeros(Na)
        lT = numpy.zeros(Na)
        for ii in range(1,Na):
            ll[ii] = sbi.CC.get_reorganization_energy(ii-1,ii-1)
        for aa in range(Na):
            for bb in range(Na):
                """ Here we assume no correlation between sites """
                lT[aa] += (SS[bb,aa]**4)*ll[bb]                 
            
                           
        KK = numpy.zeros((Na,Na))
        for a in range(Na):
            for b in range(Na):
                if a != b:
                    ed = HH[b,b] # donor
                    ea = HH[a,a] # acceptor
                    ld = lT[b] # donor reorganization energy
                    KK[a,b] = (JR[a,b]**2)*self._fintegral(ta,
                                                           gt[a,:],gt[b,:],
                                                           ed,ea,ld)
        #
        # Transfer rates
        #                                                          
        for aa in range(Na):
            for bb in range(Na):
                if aa != bb:
                    self.data[aa,aa,bb,bb] = KK[aa,bb]
           
        # calculates dephasing rates and depopulation rates
        self.updateStructure()
        self._data_initialized = True    
