# -*- coding: utf-8 -*-
import numpy

from .redfieldtensor import RedfieldRelaxationTensor
from ..corfunctions.correlationfunctions import CorrelationFunction
from ..corfunctions.cfmatrix import CorrelationFunctionMatrix
from .foerstertensor import FoersterRelaxationTensor
from ..corfunctions.correlationfunctions import c2g

class RedfieldFoersterRelaxationTensor(RedfieldRelaxationTensor):
    """Combination of Redfield and Foerster relaxation tensors
    
    Paramaters
    ----------
    ham : cu.oqs.hilbertspace.Hamiltonian
        Hamiltonian of the system
        
    sbi : cu.oqs.liouvillespace.SystemBathInteraction
        Object specifying system bath interaction
        
    initialize : bool
        If True, the tensor is imediately calculated
        
    cutoff_time : float
        Time after which the integration kernel of the Redfield tensor
        is assumed to be zero
        
    coupling_cutoff : float
        The smallest value of coupling which is still considered to cause
        delocalization.
        
    
    """
    def __init__(self, ham, sbi, initialize=True,
                 cutoff_time=None, coupling_cutoff=None):
                     
        # non-zero coupling cut-off requires a calculation of both Redfield
        # and Foerster contributions
        if coupling_cutoff is not None:
            self.coupling_cutoff = coupling_cutoff
            super().__init__(ham,sbi,initialize=False,cutoff_time=cutoff_time)
            if initialize:            
                self.__reference_implementation()
                self._is_initialized = True
        # non cut-off - standard Redfield calculation
        else:
            super().__init__(ham,sbi,initialize,cutoff_time) 
            
    def __reference_implementation(self):
        """ Reference all Python implementation
        
        """
        ham = self.Hamiltonian
        Na = ham.dim
        sbi = self.SystemBathInteraction
        ta = sbi.TimeAxis
        Nt = ta.length

        # Identify weak couplings        
        SS,JR = self.Hamiltonian.diagonalize(coupling_cutoff
                                             =self.coupling_cutoff)
        self.Hamiltonian.undiagonalize(with_remainder=False)
        calcFT = True
        # is the remainder coupling different from zero?
        if numpy.allclose(JR,numpy.zeros(JR.shape)):
            calcFT = False

        # calculate Redfield tensor for the strong coupling part
        if self._has_cutoff_time:
            RT = RedfieldRelaxationTensor(self.Hamiltonian,sbi,
                                      cutoff_time=self.cutoff_time)
        else:
            RT = RedfieldRelaxationTensor(self.Hamiltonian,sbi)
                                      
        self.data = RT.data

        if calcFT:

            self.Hamiltonian.diagonalize()
            
            # identify correlation functions of excitonic states
            cvals = numpy.zeros((Na,Nt),dtype=numpy.complex128)
            Gt = numpy.zeros((Na,ta.length),dtype=numpy.complex128)
            for ii in range(1,Na):
                Gt[ii,:] = c2g(ta,
                           sbi.CC.get_coft(ii-1,ii-1))
            for aa in range(Na):
                for bb in range(Na):
                    """ Here we assume no correlation between sites """
                    cvals[aa,:] += (ham.SS[bb,aa]**4)*Gt[bb,:]  
                    
            
            # calculate reorganization energies of exciton states
            lamb = numpy.zeros(Na-1)
            lamb_sites = numpy.zeros(Na-1)
            for ii in range(Na-1):
                lamb_sites[ii] = sbi.CC.get_reorganization_energy(ii,ii)

            #FIXME: Transformation of the reorganization energies
            #       and correlation functions should be done here, not in
            #       the FoersterRelaxationTensor

            cfm = CorrelationFunctionMatrix(ta,Na,Na)
            for ii in range(Na-1):
                params = dict(ftype="Value-defined",reorg=lamb[ii])
                fc = CorrelationFunction(ta,params,values=cvals[ii,:])
                cfm.set_correlation_function(fc,[(ii,ii)],ii+1)

            if self._has_cutoff_time:
                FT = FoersterRelaxationTensor(self.Hamiltonian,sbi,
                                        cutoff_time=self.cutoff_time)
            else:
                FT = FoersterRelaxationTensor(self.Hamiltonian,sbi)                        

            self.data += FT.data  
            
            self.Hamiltonian.undiagonalize(with_remainder=False)
            
        self._is_initialized = True
        
