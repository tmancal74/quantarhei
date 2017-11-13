# -*- coding: utf-8 -*-
"""

Lindblad form implementation


"""
import numpy

from .redfieldtensor import RedfieldRelaxationTensor
from ...builders.aggregates import Aggregate

class LindbladForm(RedfieldRelaxationTensor):
    """Lindblad form of relaxation tensor

    We use the Redfield tensor class and reimplement its implementation

    """

    def __init__(self, ham, sbi, initialize=True,
                 as_operators=True, name=""):

        super().__init__(ham, sbi, initialize=initialize,
                         as_operators=as_operators, name=name)


    def _implementation(self, ham, sbi):
        """Very simple implementation of Lindblad form

        """
        # dimension of the operators
        Na = ham.dim

        # number of operators
        Nb = sbi.N

        llm = numpy.zeros((Nb, Na, Na), dtype=numpy.float64)
        lld = numpy.zeros((Nb, Na, Na), dtype=numpy.float64)
        for i in range(Nb):
            llm[i, :, :] = sbi.rates[i]*sbi.KK[i, :, :]/2.0
            lld[i, :, :] = numpy.transpose(llm[i, :, :])

        self._post_implementation(sbi.KK, llm, lld)



class ElectronicLindbladForm(LindbladForm):
    """Lindblad for corresponding to purely electronic relaxation


    Purely electronic relaxation defined in the site basis


    """

    def __init__(self, ham, sbi, initialize=True,
                 as_operators=True, name=""):

        if isinstance(sbi.system, Aggregate):
            agg = sbi.system

            if agg.Nel == agg.Ntot:
                # if the number of electronic states corresponds to the total
                # number of states, the aggregate is purely electronic. There
                # is then no difference between Lindblad form and its
                # electronic version. We call super()
                super().__init__(ham, sbi, initialize=initialize,
                                 as_operators=as_operators, name=name)
            else:
                # check the sbi for rates and (electronic) operators
                pass

        else:

            raise Exception("SystemBathInteraction has to have `system`"+
                            " attribute set to an Aggregate")
