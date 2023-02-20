# -*- coding: utf-8 -*-
import numpy
#import scipy.interpolate as interp

from ..hilbertspace.hamiltonian import Hamiltonian
from ..liouvillespace.systembathinteraction import SystemBathInteraction
from .relaxationtensor import RelaxationTensor
from .rates.modifiedredfieldrates import ModifiedRedfieldRateMatrix
from ...core.managers import energy_units
from ...core.managers import eigenbasis_of
from ... import COMPLEX

class ModRedfieldRelaxationTensor(RelaxationTensor):
    """Modified Redfield Theory
    
    
    
    """
    def __init__(self, ham, sbi, initialize=True, cutoff_time=None,
                 as_operators=False):
        
        if as_operators:
            import warnings
            warnings.warn("Modiefied Redfield Tensor does not have "
                          +"operator form: using tensor form.")
            as_operators = False
            
            
        self._initialize_basis()
        
        if not isinstance(ham, Hamiltonian):
            raise Exception("First argument must be a Hamiltonian")
            
        if not isinstance(sbi, SystemBathInteraction):
            raise Exception("Second argument must be of"
                           +" type SystemBathInteraction")
            
        self._is_initialized = False
        self._has_cutoff_time = False
        self.as_operators = False
        
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
            

    def initialize(self):
        
        #
        # Tensor data
        #
        Na = self.dim
        self.data = numpy.zeros((Na,Na,Na,Na),dtype=COMPLEX)
        
        
        with energy_units("int"):
            
            HH = self.Hamiltonian
            sbi = self.SystemBathInteraction             
            
            #HH.protect_basis()
            #with eigenbasis_of(HH):
            if True:
                if self._has_cutoff_time:
                    cft = self.cut_off_time
                else:
                    cft = None
                
                frm = ModifiedRedfieldRateMatrix(HH, sbi, # sbi.TimeAxis,
                                                 initialize=True,
                                                 cutoff_time=cft)
    
            with eigenbasis_of(HH):
                #
                # Transfer rates
                #                                                          
                for aa in range(self.dim):
                    for bb in range(self.dim):
                        if aa != bb:
                            self.data[aa,aa,bb,bb] = frm.data[aa,bb]
                    
                #  
                # calculate dephasing rates and depopulation rates
                #
                self.updateStructure()
                
            #HH.unprotect_basis()

