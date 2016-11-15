# -*- coding: utf-8 -*-
import numpy

from .redfieldtensor import RedfieldRelaxationTensor
from ..corfunctions.correlationfunctions import CorrelationFunction
from ..corfunctions.cfmatrix import CorrelationFunctionMatrix
from .foerstertensor import FoersterRelaxationTensor
from ..corfunctions.correlationfunctions import c2g
from ...core.managers import Manager, eigenbasis_of, energy_units
from ..hilbertspace.hamiltonian import Hamiltonian
from ..hilbertspace.operators import Operator
from ..liouvillespace.systembathinteraction import SystemBathInteraction
from ...core.units import cm2int

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
            m = Manager()
            self.coupling_cutoff = m.convert_energy_2_internal_u(coupling_cutoff)
        else:
            self.coupling_cutoff = 0.0
            
        super().__init__(ham, sbi, initialize=False, 
                             cutoff_time=cutoff_time)
        if initialize: 
            with energy_units("int"):
                self._reference_implementation()
                
    def _reference_implementation(self):
        """ Reference all Python implementation
        
        """
        ham = self.Hamiltonian 
        sbi = self.SystemBathInteraction
        
        ta = sbi.TimeAxis
        Nt = ta.length
        Na = ham.dim

        if ham._has_remainder_coupling:
            JR = ham.JR
        else:
            JR = numpy.zeros((ham.dim,ham.dim), dtype=numpy.float64)
        
        calcRT = True
        calcFT = False

        # is the remainder coupling different from zero?
        if numpy.allclose(JR,numpy.zeros(JR.shape)):
            calcFT = False

        # create empty data
        self.data = numpy.zeros((ham.dim,ham.dim,ham.dim,ham.dim),
                                dtype=numpy.complex128)

        #
        # calculate Redfield tensor for the strong coupling part
        #
        if calcRT:
            
            if self._has_cutoff_time:
                RT = RedfieldRelaxationTensor(ham, sbi,
                                      cutoff_time=self.cutoff_time)
            else:
                RT = RedfieldRelaxationTensor(ham, sbi)
        
            self.data += RT.data
        
        
        #
        # Calculate Foerster for the remainder coupling
        #
        if calcFT:

            hD, SS = numpy.linalg.eigh(ham.data) 
           
            # identify correlation functions of excitonic states
            cvals = numpy.zeros((Na,Nt),dtype=numpy.complex128)
            Gt = numpy.zeros((Na,ta.length),dtype=numpy.complex128)
            for ii in range(1,Na):
                Gt[ii,:] = c2g(ta, sbi.CC.get_coft(ii-1,ii-1))
            for aa in range(Na):
                for bb in range(Na):
                    # Here we assume no correlation between sites 
                    cvals[aa,:] += (SS[bb,aa]**4)*Gt[bb,:]  
           
            # calculate reorganization energies of exciton states
            lamb = numpy.zeros(Na-1)
            lamb_sites = numpy.zeros(Na-1)
            for ii in range(Na-1):
                lamb_sites[ii] = sbi.CC.get_reorganization_energy(ii,ii)
            for aa in range(Na-1):
                for bb in range(Na-1):
                    # Here we assume no correlation between sites 
                    lamb[aa] += (SS[bb,aa]**4)*lamb_sites[bb]
                    
            #FIXME: Transformation of the reorganization energies
            #       and correlation functions should be done here, not in
            #       the FoersterRelaxationTensor

            #
            # create a new system-bath interaction object
            #

            # operators
            nsbi_op = []
            for aa in range(1,Na):
                op = Operator(dim=ham.dim, real=True)
                op.data[aa,aa] = 1.0
                nsbi_op.append(op)
            
            # correlation function matrix
            cfm = CorrelationFunctionMatrix(ta, Na-1,Na-1)
            for ii in range(Na-1):
                params = dict(ftype="Value-defined",reorg=lamb[ii])
                fc = CorrelationFunction(ta,params,values=cvals[ii,:])
                cfm.set_correlation_function(fc,[(ii,ii)],ii+1)
            nsbi = SystemBathInteraction(nsbi_op, cfm)

            if self._has_cutoff_time:
                FT = FoersterRelaxationTensor(ham, nsbi,
                                        cutoff_time=self.cutoff_time)
            else:
                FT = FoersterRelaxationTensor(ham, sbi)                        
                
            self.data += FT.data
            

            
        self._is_initialized = True
        
