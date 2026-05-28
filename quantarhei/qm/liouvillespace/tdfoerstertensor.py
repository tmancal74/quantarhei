from __future__ import annotations

import numpy
import numpy.typing

# import matplotlib.pyplot as plt
from ...core.managers import energy_units
from ...core.time import TimeDependent
from ..corfunctions.correlationfunctions import c2g
from ..hilbertspace.hamiltonian import Hamiltonian
from ..liouvillespace.systembathinteraction import SystemBathInteraction
from .foerstertensor import FoersterRelaxationTensor
from .rates.tdfoersterrates import _td_fintegral, _td_reference_implementation


class TDFoersterRelaxationTensor(FoersterRelaxationTensor, TimeDependent):
    """Weak resonance coupling relaxation tensor by Foerster theory

    Time-dependent version



    """

    def __init__(
        self,
        ham: Hamiltonian,
        sbi: SystemBathInteraction,
        initialize: bool = True,
        cutoff_time: float | None = None,
    ) -> None:

        super().__init__(ham, sbi, initialize, cutoff_time)

        self.is_time_dependent = True

    def initialize(self) -> None:

        sbi_ta = self.SystemBathInteraction.TimeAxis
        assert sbi_ta is not None, "SystemBathInteraction must have a TimeAxis set"
        tt = sbi_ta.data
        Nt = len(tt)
        #
        # Tensor data
        #
        Na = self.dim
        self.data = numpy.zeros((Nt, Na, Na, Na, Na), dtype=numpy.complex128)

        with energy_units("int"):
            # Hamiltonian matrix
            HH = self.Hamiltonian.data

            sbi = self.SystemBathInteraction
            Na = self.dim

            ta = sbi_ta
            cc = sbi.CC
            assert cc is not None, "SystemBathInteraction must have CC set"

            # line shape functions
            gt = numpy.zeros((Na, ta.length), dtype=numpy.complex64)

            # SBI is defined with "sites"
            for ii in range(1, Na):
                gt[ii, :] = c2g(ta, cc.get_coft(ii - 1, ii - 1))

            # reorganization energies
            ll = numpy.zeros(Na)
            for ii in range(1, Na):
                ll[ii] = cc.get_reorganization_energy(ii - 1, ii - 1)

            KK = self.td_reference_implementation(Na, Nt, HH, tt, gt, ll)

            #
            # Transfer rates
            #
            for aa in range(self.dim):
                for bb in range(self.dim):
                    if aa != bb:
                        self.data[:, aa, aa, bb, bb] = KK[:, aa, bb]

            #
            # calculate dephasing rates and depopulation rates
            #
            self.updateStructure()

            # additional pure dephasing
            self.add_dephasing()

    def add_dephasing(self) -> None:

        # line shape function derivatives
        sbi = self.SystemBathInteraction
        Na = self.dim

        ta = sbi.TimeAxis
        assert ta is not None
        cc = sbi.CC
        assert cc is not None

        ht = numpy.zeros((Na, ta.length), dtype=numpy.complex64)

        cc.create_one_integral()

        for ii in range(1, Na):
            ht[ii, :] = cc.get_hoft(ii - 1, ii - 1)

        for aa in range(self.dim):
            for bb in range(self.dim):
                if aa != bb:
                    self.data[:, aa, bb, aa, bb] -= ht[aa, :] + ht[bb, :]

    def td_reference_implementation(
        self,
        Na: int,
        Nt: int,
        HH: numpy.ndarray,
        tt: numpy.ndarray,
        gt: numpy.ndarray,
        ll: numpy.ndarray,
    ) -> numpy.ndarray:
        return _td_reference_implementation(Na, Nt, HH, tt, gt, ll, _td_fintegral)
