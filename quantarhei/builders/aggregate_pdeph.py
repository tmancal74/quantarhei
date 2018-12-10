# -*- coding: utf-8 -*-
"""Class adding calculation of PureDephasing object

    Most of the present methods are available after the aggregate is
    diagonalized by calling the ``diagonalize`` method.

    **This class should not be used directly**. Use `Aggregate` class, which
    inherits all the methods from here, instead.
    

    Class Details
    -------------

"""
import numpy

from .aggregate_excitonanalysis import AggregateExcitonAnalysis

from .. import REAL

class AggregatePureDephasing(AggregateExcitonAnalysis):
    """Class calculation of PureDephasing object
    
    
    
    """

    def get_PureDephasing(self, dtype="Lorentzian"):
        """Returns pure dephasing object of this aggregate
        
        """
        
        from ..qm.liouvillespace.puredephasing import ElectronicPureDephasing
        
        # collect site basis dephasing rates
        
        pdrates = numpy.zeros(self.nmono, dtype=REAL)
        k = 0
        if dtype == "Lorentzian":
            
            for mono in self.monomers:
                if mono.dephs is not None:
                    # electronic optical dephasing rates from monomers
                    pdrates[k] = mono.dephs[0, 1]
                else:
                    pdrates[k] = 0.0

        elif dtype == "Gaussian":
            for mono in self.monomers:
                if mono.dephs is not None:
                    # electronic optical dephasing rates from monomers
                    pdrates[k] = mono.widths[0, 1]
                else:
                    pdrates[k] = 0.0
        else:
            raise Exception("Unknown dephasing type")
            
        self.diagonalize()
        
        # Na is the number of states (vibronic origin)
        Na = self.Ntot
        # self.Nel number of electronic states
        xiai = numpy.zeros((Na, self.Nel), dtype=REAL)
        for aa in range(Na):
            st = 0 # st counts all states in the site basis
            for ii in range(self.Nel):
                xiai[aa,ii] = 0.0
                for alph_i in self.vibindices[ii]:
                    xiai[aa,ii] += self.SS[aa,st]**2
                    st += 1
        self.xi = xiai
        
        # return PureDephasing object
        return ElectronicPureDephasing(self, dtype=dtype)