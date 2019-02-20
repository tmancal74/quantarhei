# -*- coding: utf-8 -*-
"""
    Pure dephasing Tensor
    
    
    Class Details
    -------------
    
"""

#
# We create this intentionally as not basis managed. It must be applied
# only while working in the correct basis!
#
import numpy

from ...builders.aggregates import Aggregate
from ...core.managers import eigenbasis_of
from ... import REAL

class PureDephasing: #(BasisManaged):
    
    
    def __init__(self, drates=None, dtype="Lorentzian", cutoff_time=None):
        
        self.data=drates
        self.dtype=dtype
        self.cutoff_time=cutoff_time
        
        
        
        
class ElectronicPureDephasing(PureDephasing):
    """Electronic pure dephasing for one-exciton states
    
    """
    
    def __init__(self, system, drates=None, dtype="Lorentzian",
                 cutoff_time=None):
        
        super().__init__(drates, dtype, cutoff_time)
        
        self.system = system
        self.system_is_aggregate = False
        
        if isinstance(self.system, Aggregate):
            
            self.system_is_aggregate = True
            
        else:
            
            raise Exception("Non-Aggregate systems not implemented yet")
            
            
        if self.system_is_aggregate:
            
            Nstates = self.system.Ntot
            Nel = self.system.number_of_electronic_states_in_band(1)+1 
            self.data = numpy.zeros((Nstates, Nstates), dtype=REAL)
            
            widths = numpy.zeros(Nel, dtype=REAL)
            for ii in range(Nel):
                if ii > 0:
                    widths[ii] = self.system.monomers[ii
                                      -1].get_transition_width((0,1))
                
            self.system.diagonalize()
            
            Xi = self.system.Xi
            
            for aa in range(Nstates):
                for bb in range(Nstates):
                    for ii in range(Nel):
                        self.data[aa,bb] += \
                        widths[ii]*(Xi[aa,ii]**2 - Xi[bb,ii]**2)**2

            
    def eigenbasis(self):
        """Returns the context for the eigenbasis in which pure dephasing is defined
        
        
        To be used as
        
        pd = PureDephasing(sys)
        with pd.eigenbasis:
            ...
        
        """
        ham = self.system.get_Hamiltonian()
        return eigenbasis_of(ham)
            