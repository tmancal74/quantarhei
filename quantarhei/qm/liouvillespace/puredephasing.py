# -*- coding: utf-8 -*-
"""
    Pure dephasing Tensor
    
    
    We create this class intentionally as not basis managed. Pure dephasing
    must be applied only while working in the correct basis!    
    
    
    Class Details
    -------------
    
    
    
    
    
    Examples
    --------
    
    >>> from quantarhei import energy_units
    >>> from quantarhei import Molecule
    
    >>> d_rates = [[0.0, 1.0/100.0], [1.0/100.0, 0.0]]
    >>> pd = PureDephasing(d_rates)
    >>> print(pd.data)
    [[ 0.    0.01]
     [ 0.01  0.  ]]
    
    >>> with energy_units("1/cm"):
    ...     mol1 = Molecule([0.0, 12000.0])
    ...     mol2 = Molecule([0.0, 12100.0])
    ...     agg = Aggregate(molecules=[mol1, mol2])
    ...     agg.set_resonance_coupling(0,1,100.0)
    >>> agg.build()
    >>> pd = PureDephasing(drates=d_rates, system=agg)
    Traceback (most recent call last):
        ...
    Exception: Incompatible dimension of the rate matrix: system has dimension = 3
    
    
    >>> d_rates = [[0.0, 0.0, 0.0],[0.0, 0.0, 1.0/100.0], [0.0, 0.0, 0.0]]
    >>> pd = PureDephasing(drates=d_rates, system=agg)
    >>> pd.data[1,2] == d_rates[1][2]
    True
    
    >>> print(pd.data[1,2])
    0.01
    
    
    >>> pd = PureDephasing(system=agg)
    Traceback (most recent call last):
        ... 
    Exception: Dephasing rates must be specified.
    
    >>> pd = PureDephasing(drates=d_rates, system=1.0)
    Traceback (most recent call last):
        ...    
    Exception: Non-Aggregate systems not implemented yet.
    
    >>> pd = PureDephasing(drates=d_rates, system=agg)
    >>> HH = agg.get_Hamiltonian()
    >>> with pd.eigenbasis:
    ...     print(numpy.max(numpy.abs(HH.data
    ...                   -numpy.diag(numpy.diag(HH.data))))<1.0e-15)
    True
    
    >>> pd.eigenbasis = None
    Traceback (most recent call last):
        ...
    Exception: The property 'eigenbasis' is protected and cannot be set.
    
"""

#

#
import numpy

from ...builders.aggregates import Aggregate
#from ...builders.molecules import Molecule
from ...core.managers import eigenbasis_of
from ... import REAL
from ..liouvillespace.superoperator import SuperOperator
from ...qm.hilbertspace.operators import Operator


def _eigenb():
    """Pointer to eigenbasis property
    
    """
    
    @property
    def prop(self):
        return self._eigenbasis()
 
    @prop.setter
    def prop(self, value):
        raise Exception("The property 'eigenbasis' is protected"+
                        " and cannot be set.")   
 
    return prop




class PureDephasing: #(BasisManaged):
    
    dtypes = ["Lorentzian", "Gaussian"]
    
    eigenbasis = _eigenb()
    
    def __init__(self, drates=None, dtype="Lorentzian", system=None, 
                 cutoff_time=None):

        if drates is None:
            raise Exception("Dephasing rates must be specified.")
           
        self.data=numpy.array(drates, dtype=REAL)    

        if system is not None:
            
            self.system = system 
            
            
            self.system_is_aggregate = False
            self.system_is_molecule = False
            
            if isinstance(self.system, Aggregate):
                
                self.system_is_aggregate = True
                
            # FIXME: we have a circular reference problem with Molecule
            #if isinstance(self.system, Molecule):
            #    
            #    self.system_is_molecule = True
                
            else:
                
                raise Exception("Non-Aggregate systems not implemented yet.")

            HH = self.system.get_Hamiltonian()
            if (HH.dim != self.data.shape[0]):
                raise Exception("Incompatible dimension of the rate matrix:"
                                +" system has dimension = "+str(HH.dim))
       
        else:
            
            self.system = None


        if dtype in self.dtypes:
        
            
            self.dim = self.data.shape[0]
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


    def _eigenbasis(self):
        """Returns the context for the eigenbasis in which pure dephasing
        is defined
        
        
        To be used as
        
        pd = PureDephasing(sys)
        with pd.eigenbasis:
            ...
        
        """
        if self.system is not None:
            
            ham = self.system.get_Hamiltonian()
            return eigenbasis_of(ham)
        
        else:
            
            op = Operator(dim=self.dim)
            return eigenbasis_of(op)
            
            


    
class ElectronicPureDephasing(PureDephasing):
    """Electronic pure dephasing for one-exciton states
    
    """
    
    def __init__(self, system, drates=None, dtype="Lorentzian",
                 cutoff_time=None):
        
        if drates is None:
            HH = system.get_Hamiltonian()
            dim = HH.dim
            rates = numpy.zeros((dim, dim), dtype=REAL)
        else:
            rates = drates
        super().__init__(rates, dtype, system, cutoff_time)
        
            
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

            

            