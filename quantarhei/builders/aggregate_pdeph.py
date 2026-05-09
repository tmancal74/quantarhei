"""Class adding calculation of PureDephasing object

Most of the present methods are available after the aggregate is
diagonalized by calling the ``diagonalize`` method.

**This class should not be used directly**. Use `Aggregate` class, which
inherits all the methods from here, instead.


Class Details
-------------

"""

from __future__ import annotations

import numpy

from .. import REAL
from .aggregate_excitonanalysis import AggregateExcitonAnalysis


class AggregatePureDephasing(AggregateExcitonAnalysis):
    """Class calculation of PureDephasing object"""

    def get_PureDephasing(self, dtype: str = "Lorentzian") -> object:
        """Returns pure dephasing object of this aggregate"""
        from ..qm.liouvillespace.puredephasing import ElectronicPureDephasing

        if dtype not in ("Lorentzian", "Gaussian"):
            raise Exception("Unknown dephasing type")

        self.diagonalize()

        Na = self.Ntot
        xiai: numpy.ndarray = numpy.zeros((Na, self.Nel), dtype=REAL)  # type: ignore[explicit-any]
        for aa in range(Na):
            st = 0
            for ii in range(self.Nel):
                xiai[aa, ii] = 0.0
                for alph_i in self.vibindices[ii]:
                    xiai[aa, ii] += self.SS[aa, st] ** 2
                    st += 1
        self.xi = xiai

        return ElectronicPureDephasing(self, dtype=dtype)
