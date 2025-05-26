# -*- coding: utf-8 -*-
import numpy
#import scipy.interpolate as interp

from ..hilbertspace.hamiltonian import Hamiltonian
from ..liouvillespace.systembathinteraction import SystemBathInteraction
from .relaxationtensor import RelaxationTensor
from .rates.foersterrates import FoersterRateMatrix
from ...core.managers import energy_units
from ... import COMPLEX

class FoersterRelaxationTensor(RelaxationTensor):
    """Weak resonance coupling relaxation tensor by Foerster theory
    
    
    
    """
    def __init__(self, ham, sbi, initialize=True, cutoff_time=None, 
                 pure_dephasing=False):
        
        self._initialize_basis()
        
        if not isinstance(ham, Hamiltonian):
            raise Exception("First argument must be a Hamiltonian")
            
        if not isinstance(sbi, SystemBathInteraction):
            raise Exception("Second argument must be of"
                           +" type SystemBathInteraction")
            
        self._is_initialized = False
        self._has_cutoff_time = False
        self.as_operators = False
        self.pure_dephasing = pure_dephasing
        self.has_Iterm = False
        
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
            
            if self._has_cutoff_time:
                cft = self.cut_off_time
            else:
                cft = None
                
            frm = FoersterRateMatrix(HH, sbi, initialize=True,
                                     cutoff_time=cft)
    
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

            if self.pure_dephasing:
                # additional pure dephasing 
                self.add_dephasing()
               
          
    def add_dephasing(self):
      

        # line shape function derivatives
        sbi = self.SystemBathInteraction
        Na = self.dim
        Nt = sbi.TimeAxis.length
        
        ht = numpy.zeros((Na, Nt),
                         dtype=numpy.complex64)
        
        sbi.CC.create_one_integral()
    
        for ii in range(1, Na):
           ht[ii,:] = sbi.CC.get_hoft(ii-1,ii-1)
        
                        
        for aa in range(self.dim):
            for bb in range(self.dim):
                if aa != bb:
 
                    self.data[aa,bb,aa,bb] -= (ht[aa,Nt-1]+ht[bb,Nt-1])
