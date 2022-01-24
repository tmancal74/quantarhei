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
from ..liouvillespace.superoperator import SuperOperator

class PureDephasing: #(BasisManaged):
    
    dtypes = ["Lorentzian", "Gaussian"]
    
    def __init__(self, drates=None, dtype="Lorentzian", cutoff_time=None):
        
        if dtype in self.dtypes:
        
            self.data=drates
            self.dtype=dtype
            self.cutoff_time=cutoff_time
            
        else:
            raise Exception("Unknown dephasing type")
        
        
    def get_SuperOperator(self):
        """Returns a superoperator representing the pure dephasing

        The return superoperator is always time-independent, i.e. in the
        case of "Gaussian" dephasing, only the constants are returned.

        """        
        dim = self.data.shape[0] 
        sup = SuperOperator(dim=dim, real=True)
        for aa in range(dim):
            for bb in range(dim):
                sup.data[aa,bb,aa,bb] = self.data[aa,bb]
        
        return sup

        
    def convert_to(self, dtype=None):
        """Converts between Lorenzian and Gaussian dephasings      
        
        The conversion is done approximatively, so that the FWHM of 
        the corresponding lineshapes are the same.
        
        Parameters
        ----------
        
        dtype: str
            Type of PureDephasing to which one should convert. Conversion
            factor is such that the FFT of the time evolution gives a curve
            with the same FWHM.
    
        """
        
        
        if dtype in self.dtypes:
            
            #factor = 2.0*numpy.sqrt(numpy.log(2.0))
            factor = numpy.sqrt(numpy.log(2.0))
            if dtype == "Lorentzian" and self.dtype == "Gaussian":
                self.data = numpy.sqrt(self.data)*factor
                self.dtype = dtype
            elif dtype == "Gaussian" and self.dtype == "Lorenzian":
                self.data = (self.data**2)/(factor**2)
                self.dtype = dtype
                
        else:
            raise Exception("Unknown dephasing type")

    
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
                    widths[ii] = (((self.system.monomers[ii
                                      -1].get_transition_width((0,1)))**2)/
                                  (8.0*numpy.log(2.0)))
                
            self.system.diagonalize()
            
            Xi = self.system.Xi
            
            for aa in range(Nstates):
                for bb in range(Nstates):
                    for ii in range(Nel):
                        self.data[aa,bb] += \
                        widths[ii]*(Xi[aa,ii]**2 - Xi[bb,ii]**2)**2

            
    def eigenbasis(self):
        """Returns the context for the eigenbasis in which pure dephasing
        is defined
        
        
        To be used as
        
        pd = PureDephasing(sys)
        with pd.eigenbasis:
            ...
        
        """
        ham = self.system.get_Hamiltonian()
        return eigenbasis_of(ham)
            