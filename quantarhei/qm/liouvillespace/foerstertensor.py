# -*- coding: utf-8 -*-
import numpy
#import scipy.interpolate as interp

from ..hilbertspace.hamiltonian import Hamiltonian
from ..liouvillespace.systembathinteraction import SystemBathInteraction
from .relaxationtensor import RelaxationTensor
from .rates.foersterrates import FoersterRateMatrix
#from ..corfunctions.correlationfunctions import c2g
from ...core.managers import energy_units

class FoersterRelaxationTensor(RelaxationTensor):
    """Weak resonance coupling relaxation tensor by Foerster theory
    
    
    
    """
    def __init__(self, ham, sbi, initialize=True, cutoff_time=None):
        
        #super().__init__()
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
        self.data = numpy.zeros((Na,Na,Na,Na),dtype=numpy.complex128)
        
        
        with energy_units("int"):
            
            # Hamiltonian matrix
            #HH = self.Hamiltonian.data
            sbi = self.SystemBathInteraction             
            
            if self._has_cutoff_time:
                cft = self.cut_off_time
            else:
                cft = None
                
            frm = FoersterRateMatrix(self.Hamiltonian, sbi, initialize=True,
                                     cutoff_time=cft)
 
            #tt = sbi.TimeAxis.data
            
            # line shape functions
            #gt = numpy.zeros((Na, sbi.TimeAxis.length),
            #                 dtype=numpy.complex64)
    
            # SBI is defined with "sites"
            #for ii in range(1, Na):
            #    gt[ii,:] = c2g(sbi.TimeAxis, sbi.CC.get_coft(ii-1,ii-1))
            
            # reorganization energies
            #ll = numpy.zeros(Na)
            #for ii in range(1, Na):
            #    ll[ii] = sbi.CC.get_reorganization_energy(ii-1,ii-1)
                
            #KK = _reference_implementation(Na, HH, tt, gt, ll)
   
            KK = frm.data
    
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
        
        


