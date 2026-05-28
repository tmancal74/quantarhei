"""Class comprising the aggregate methods for support of spectroscopic simulations



Class Details
-------------

"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy

from ..core.managers import eigenbasis_of
from ..exceptions import ImplementationError, QuantarheiError
from ..qm.liouvillespace.supopunity import SOpUnity
from ..spectroscopy import diagramatics as diag
from .aggregate_base import AggregateBase


class AggregateSpectroscopy(AggregateBase):
    """Class comprising the aggregate methods for support of spectroscopic simulations"""

    ########################################################################
    #
    #   SPECTROSCOPY
    #
    ########################################################################

    def liouville_pathways_3(
        self,
        ptype: str | tuple | list = "R3g",
        dtol: float = 0.01,
        ptol: float = 1.0e-3,
        lab: object = None,
        verbose: int = 0,
    ) -> list:
        """Generator of Liouville pathways"""
        ham = self.get_Hamiltonian()
        self.lab = lab
        return self.liouville_pathways_3T(
            ptype,
            dtol=dtol,
            ptol=ptol,
            lab=lab,
            eUt=SOpUnity(dim=ham.dim),
            verbose=verbose,
        )

    def liouville_pathways_3T(
        self,
        ptype: str | tuple | list = "R3g",
        eUt: Any = None,
        ham: Any = None,
        t2: float = 0.0,
        dtol: float = 1.0e-12,
        ptol: float = 1.0e-3,
        etol: float = 1.0e-6,
        verbose: int = 0,
        lab: Any = None,
    ) -> list:
        """Generator of Liouville pathways with energy transfer




        Parameters
        ----------
        ptype : tuple, list, str
            List of strings or a string representing one or more
            Liouville pathway types that are to be calculated

        eUt : EvolutionSuperOperator
            Evolution superoperator representing the energy
            transfer in the system

        t2 : float
            Waiting time at which the spectrum is calculated

        dtol : float
            Minimum acceptable strength of the transition from ground
            to excited state, relative to the maximum dipole strength
            available in the system

        ptol : float
            Minimum acceptable population of the ground state (e.g. states
            not thermally populated are excluded)

        lab : LaboratorySetup
            Object representing laboratory setup - number of pulses,
            polarization etc.

        Returns
        -------
        lst : list
            List of LiouvillePathway objects


        """
        self.lab = lab

        if self._diagonalized:
            if verbose > 0:
                print("Diagonalizing aggregate")
            self.diagonalize()
            if verbose > 0:
                print("..done")

        pop_tol = ptol
        dip_tol = numpy.sqrt(self.D2_max) * dtol
        evf_tol = etol

        # Check if the ptype is a tuple
        ptype_tuple: Any
        if not isinstance(ptype, (tuple, list)):
            ptype_tuple = (ptype,)
        else:
            ptype_tuple = ptype
        lst: list[Any] = []

        if verbose > 0:
            print("Pathways", ptype_tuple)

        #
        # data of the evolution superoperator in eigenstate basis
        #

        try:
            # either the eUt is a complete evolution superoperator
            eUt2 = eUt.at(t2)

            #
            # TODO: Check this part I had before without eigenbasis of and used the evolution superoperator values directly
            # ----------------------------------------------
            eUt2_dat = numpy.zeros(eUt2.data.shape, dtype=eUt2.data.dtype)  #
            HH = eUt.get_Hamiltonian()  #
            with eigenbasis_of(HH):  #
                eUt2_dat[:, :, :, :] = eUt2.data  #
        # ----------------------------------------------

        except AttributeError:
            # or it is only a super operator at a given time t2
            # in this case 'ham' must be specified
            eUt2 = eUt
            eUt2_dat = numpy.zeros(eUt2.data.shape, dtype=eUt2.data.dtype)
            with eigenbasis_of(ham):
                eUt2_dat[:, :, :, :] = eUt2.data

        for ptp in ptype_tuple:
            if ptp == "R1g":
                generate_R1g(self, lst, eUt2_dat, pop_tol, dip_tol, evf_tol, verbose)

            elif ptp == "R2g":
                generate_R2g(self, lst, eUt2_dat, pop_tol, dip_tol, evf_tol, verbose)

            elif ptp == "R3g":
                generate_R3g(self, lst, eUt2_dat, pop_tol, dip_tol, verbose)

            elif ptp == "R4g":
                generate_R4g(self, lst, eUt2_dat, pop_tol, dip_tol, verbose)

            elif ptp == "R1f*":
                generate_R1f(self, lst, eUt2_dat, pop_tol, dip_tol, evf_tol, verbose)

            elif ptp == "R2f*":
                generate_R2f(self, lst, eUt2_dat, pop_tol, dip_tol, evf_tol, verbose)

            elif ptp == "R1gE":
                generate_R1gE(self, lst, eUt2_dat, pop_tol, dip_tol, evf_tol, verbose)

            elif ptp == "R2gE":
                generate_R2gE(self, lst, eUt2_dat, pop_tol, dip_tol, evf_tol, verbose)

            elif ptp == "R1f*E":
                generate_R1fE(self, lst, eUt2_dat, pop_tol, dip_tol, evf_tol, verbose)

            elif ptp == "R2f*E":
                generate_R2fE(self, lst, eUt2_dat, pop_tol, dip_tol, evf_tol, verbose)

            else:
                raise QuantarheiError("Unknown pythway type: " + str(ptp))

        if lab is not None:
            for l in lst:
                l.orientational_averaging(lab)

        return lst

    def liouville_pathways_1(
        self,
        eUt: object = None,
        ham: object = None,
        dtol: float = 0.01,
        ptol: float = 1.0e-3,
        etol: float = 1.0e-6,
        verbose: int = 0,
        lab: object = None,
    ) -> list:
        """Generator of the first order Liouville pathways


        Generator of the pathways for an absorption spectrum
        calculation.



        Parameters
        ----------
        eUt : EvolutionSuperOperator
            Evolution superoperator representing the evolution of optical
            coherence in the system


        dtol : float
            Minimum acceptable strength of the transition from ground
            to excited state, relative to the maximum dipole strength
            available in the system

        ptol : float
            Minimum acceptable population of the ground state (e.g. states
            not thermally populated are excluded)

        lab : LaboratorySetup
            Object representing laboratory setup - number of pulses,
            polarization etc.

        Returns
        -------
        lst : list
            List of LiouvillePathway objects


        """
        if self._diagonalized:
            if verbose > 0:
                print("Diagonalizing aggregate")
            self.diagonalize()
            if verbose > 0:
                print("..done")

        pop_tol = ptol
        dip_tol = numpy.sqrt(self.D2_max) * dtol
        evf_tol = etol

        if eUt is None:
            # secular absorption spectrum calculation
            eUt2_dat = None
            sec = True

        else:
            raise ImplementationError("Not implemented yet")

        lst: list[Any] = []

        if sec:
            generate_1orderP_sec(self, lst, pop_tol, dip_tol, verbose)
        else:
            raise ImplementationError("Not implemented yet")

        if lab is not None:
            for l in lst:
                l.orientational_averaging(lab)

        return lst


# ---------------------------------------------------------------------------
# Data-driven pathway generation (incremental refactoring of #366)
# ---------------------------------------------------------------------------

from dataclasses import dataclass


@dataclass(frozen=True)
class _GroundExcitedPathwayDesc:
    rtype: str
    pname: str
    evf_fn: Callable[[int, int], tuple[int, int, int, int]]
    inner_dip_check: Callable[..., bool]
    width1_fn: Callable[..., tuple[int, int]]
    width3_fn: Callable[..., tuple[int, int]]
    transitions: tuple[tuple[Callable[..., tuple[int, int]], int, int | None], ...]


_R3G_DESC = _GroundExcitedPathwayDesc(
    rtype="R",
    pname="R3g",
    evf_fn=lambda i1g, i3g: (i1g, i3g, i1g, i3g),
    inner_dip_check=lambda self, i4e, i1g, i3g, dip_tol: (
        self.D2[i4e, i1g] > dip_tol and self.D2[i3g, i4e] > dip_tol
    ),
    width1_fn=lambda i2e, i1g, **_: (i2e, i1g),
    width3_fn=lambda i4e, i3g, **_: (i4e, i3g),
    transitions=(
        # (idx_fn returning (a,b), side, interval)
        (lambda i2e, i1g, **_: (i2e, i1g), -1, 1),
        (lambda i3g, i2e, **_: (i3g, i2e), -1, None),
        (lambda i4e, i1g, **_: (i4e, i1g), +1, None),
        (lambda i3g, i4e, **_: (i3g, i4e), +1, 3),
    ),
)

_R4G_DESC = _GroundExcitedPathwayDesc(
    rtype="NR",
    pname="R4g",
    evf_fn=lambda i1g, i3g: (i1g, i3g, i1g, i3g),
    inner_dip_check=lambda self, i4e, i1g, i3g, dip_tol: (
        self.D2[i4e, i3g] > dip_tol and self.D2[i1g, i4e] > dip_tol
    ),
    width1_fn=lambda i2e, i1g, **_: (i2e, i1g),
    width3_fn=lambda i4e, i1g, **_: (i4e, i1g),
    transitions=(
        (lambda i2e, i1g, **_: (i2e, i1g), +1, 1),
        (lambda i3g, i2e, **_: (i3g, i2e), +1, None),
        (lambda i4e, i3g, **_: (i4e, i3g), +1, None),
        (lambda i1g, i4e, **_: (i1g, i4e), +1, 3),
    ),
)


@dataclass(frozen=True)
class _RelaxationPathwayDesc:
    rtype: str
    pname: str
    requires_band2: bool
    relax_band: int  # 1 = nes, 0 = ngs
    inner_band: int  # 0 = ngs, 1 = nes, 2 = nfs
    evf_fn: Callable[..., tuple[int, int, int, int]]
    inner_dip_check: Callable[..., bool]
    width1_fn: Callable[..., tuple[int, int]]
    width3_fn: Callable[..., tuple[int, int]]
    transfer_fn: Callable[..., tuple[tuple[int, int], tuple[int, int]]]
    transitions_before: tuple[
        tuple[Callable[..., tuple[int, int]], int, int | None], ...
    ]
    transitions_after: tuple[
        tuple[Callable[..., tuple[int, int]], int, int | None], ...
    ]


_R1G_DESC = _RelaxationPathwayDesc(
    rtype="NR",
    pname="R1g",
    requires_band2=False,
    relax_band=1,
    inner_band=0,
    evf_fn=lambda iRa, iRb, i2e, i3e: (iRa, iRb, i2e, i3e),
    inner_dip_check=lambda self, iX, iRa, iRb, i2e, i3e, dip_tol: (
        self.D2[iX, iRb] > dip_tol and self.D2[iX, iRa] > dip_tol
    ),
    width1_fn=lambda i2e, i1g, **_: (i2e, i1g),
    width3_fn=lambda iRa, iX, **_: (iRa, iX),
    transfer_fn=lambda iRa, iRb, i2e, i3e: ((iRa, iRb), (i2e, i3e)),
    transitions_before=(
        (lambda i2e, i1g, **_: (i2e, i1g), +1, 1),
        (lambda i3e, i1g, **_: (i3e, i1g), -1, None),
    ),
    transitions_after=(
        (lambda iX, iRb, **_: (iX, iRb), -1, None),
        (lambda iX, iRa, **_: (iX, iRa), +1, 3),
    ),
)

_R2G_DESC = _RelaxationPathwayDesc(
    rtype="R",
    pname="R2g",
    requires_band2=False,
    relax_band=1,
    inner_band=0,
    evf_fn=lambda iRa, iRb, i2e, i3e: (iRa, iRb, i3e, i2e),
    inner_dip_check=lambda self, iX, iRa, iRb, i2e, i3e, dip_tol: (
        self.D2[iX, i2e] > dip_tol and self.D2[iX, i3e] > dip_tol
    ),
    width1_fn=lambda i2e, i1g, **_: (i2e, i1g),
    width3_fn=lambda iRa, iX, **_: (iRa, iX),
    transfer_fn=lambda iRa, iRb, i2e, i3e: ((iRa, iRb), (i3e, i2e)),
    transitions_before=(
        (lambda i2e, i1g, **_: (i2e, i1g), -1, 1),
        (lambda i3e, i1g, **_: (i3e, i1g), +1, None),
    ),
    transitions_after=(
        (lambda iX, iRb, **_: (iX, iRb), -1, None),
        (lambda iX, iRa, **_: (iX, iRa), +1, 3),
    ),
)

_R1F_DESC = _RelaxationPathwayDesc(
    rtype="R",
    pname="R1f*",
    requires_band2=True,
    relax_band=1,
    inner_band=2,
    evf_fn=lambda iRa, iRb, i2e, i3e: (iRa, iRb, i3e, i2e),
    inner_dip_check=lambda self, iX, iRa, iRb, i2e, i3e, dip_tol: (
        self.D2[iX, iRa] > dip_tol and self.D2[iRb, iX] > dip_tol
    ),
    width1_fn=lambda i2e, i1g, **_: (i2e, i1g),
    width3_fn=lambda iX, iRb, **_: (iX, iRb),
    transfer_fn=lambda iRa, iRb, i2e, i3e: ((iRa, iRb), (i3e, i2e)),
    transitions_before=(
        (lambda i2e, i1g, **_: (i2e, i1g), -1, 1),
        (lambda i3e, i1g, **_: (i3e, i1g), +1, None),
    ),
    transitions_after=(
        (lambda iX, iRa, **_: (iX, iRa), +1, None),
        (lambda iRb, iX, **_: (iRb, iX), +1, 3),
    ),
)

_R2F_DESC = _RelaxationPathwayDesc(
    rtype="NR",
    pname="R2f*",
    requires_band2=True,
    relax_band=1,
    inner_band=2,
    evf_fn=lambda iRa, iRb, i2e, i3e: (iRa, iRb, i2e, i3e),
    inner_dip_check=lambda self, iX, iRa, iRb, i2e, i3e, dip_tol: (
        self.D2[iX, iRa] > dip_tol and self.D2[iRb, iX] > dip_tol
    ),
    width1_fn=lambda i2e, i1g, **_: (i2e, i1g),
    width3_fn=lambda iX, iRb, **_: (iX, iRb),
    transfer_fn=lambda iRa, iRb, i2e, i3e: ((iRa, iRb), (i2e, i3e)),
    transitions_before=(
        (lambda i2e, i1g, **_: (i2e, i1g), +1, 1),
        (lambda i3e, i1g, **_: (i3e, i1g), -1, None),
    ),
    transitions_after=(
        (lambda iX, iRa, **_: (iX, iRa), +1, None),
        (lambda iRb, iX, **_: (iRb, iX), +1, 3),
    ),
)

_R1GE_DESC = _RelaxationPathwayDesc(
    rtype="NR",
    pname="R1gE",
    requires_band2=False,
    relax_band=0,
    inner_band=1,
    evf_fn=lambda iRa, iRb, i2e, i3e: (iRa, iRb, i2e, i3e),
    inner_dip_check=lambda self, iX, iRa, iRb, i2e, i3e, dip_tol: (
        self.D2[iRa, iX] > dip_tol and self.D2[iRb, iX] > dip_tol
    ),
    width1_fn=lambda i2e, i1g, **_: (i2e, i1g),
    width3_fn=lambda iX, iRa, **_: (iX, iRa),
    transfer_fn=lambda iRa, iRb, i2e, i3e: ((iRa, iRb), (i2e, i3e)),
    transitions_before=(
        (lambda i2e, i1g, **_: (i2e, i1g), +1, 1),
        (lambda i3e, i1g, **_: (i3e, i1g), -1, None),
    ),
    transitions_after=(
        (lambda iX, iRa, **_: (iX, iRa), +1, None),
        (lambda iRb, iX, **_: (iRb, iX), +1, 3),
    ),
)

_R2GE_DESC = _RelaxationPathwayDesc(
    rtype="R",
    pname="R2gE",
    requires_band2=False,
    relax_band=0,
    inner_band=1,
    evf_fn=lambda iRa, iRb, i2e, i3e: (iRa, iRb, i3e, i2e),
    inner_dip_check=lambda self, iX, iRa, iRb, i2e, i3e, dip_tol: (
        self.D2[iRa, iX] > dip_tol and self.D2[iRb, iX] > dip_tol
    ),
    width1_fn=lambda i2e, i1g, **_: (i2e, i1g),
    width3_fn=lambda iX, iRa, **_: (iX, iRa),
    transfer_fn=lambda iRa, iRb, i2e, i3e: ((iRa, iRb), (i3e, i2e)),
    transitions_before=(
        (lambda i2e, i1g, **_: (i2e, i1g), -1, 1),
        (lambda i3e, i1g, **_: (i3e, i1g), +1, None),
    ),
    transitions_after=(
        (lambda iX, iRa, **_: (iX, iRa), +1, None),
        (lambda iRb, iX, **_: (iRb, iX), +1, 3),
    ),
)

_R1FE_DESC = _RelaxationPathwayDesc(
    rtype="R",
    pname="R1f*E",
    requires_band2=False,
    relax_band=0,
    inner_band=1,
    evf_fn=lambda iRa, iRb, i2e, i3e: (iRa, iRb, i3e, i2e),
    inner_dip_check=lambda self, iX, iRa, iRb, i2e, i3e, dip_tol: (
        self.D2[iX, iRa] > dip_tol and self.D2[iRb, iX] > dip_tol
    ),
    width1_fn=lambda i2e, i1g, **_: (i2e, i1g),
    width3_fn=lambda iX, iRb, **_: (iX, iRb),
    transfer_fn=lambda iRa, iRb, i2e, i3e: ((iRa, iRb), (i3e, i2e)),
    transitions_before=(
        (lambda i2e, i1g, **_: (i2e, i1g), -1, 1),
        (lambda i3e, i1g, **_: (i3e, i1g), +1, None),
    ),
    transitions_after=(
        (lambda iX, iRa, **_: (iX, iRa), +1, None),
        (lambda iRb, iX, **_: (iRb, iX), +1, 3),
    ),
)

_R2FE_DESC = _RelaxationPathwayDesc(
    rtype="NR",
    pname="R2f*E",
    requires_band2=False,
    relax_band=0,
    inner_band=1,
    evf_fn=lambda iRa, iRb, i2e, i3e: (iRa, iRb, i2e, i3e),
    inner_dip_check=lambda self, iX, iRa, iRb, i2e, i3e, dip_tol: (
        self.D2[iX, iRa] > dip_tol and self.D2[iRb, iX] > dip_tol
    ),
    width1_fn=lambda i2e, i1g, **_: (i2e, i1g),
    width3_fn=lambda iX, iRb, **_: (iX, iRb),
    transfer_fn=lambda iRa, iRb, i2e, i3e: ((iRa, iRb), (i2e, i3e)),
    transitions_before=(
        (lambda i2e, i1g, **_: (i2e, i1g), +1, 1),
        (lambda i3e, i1g, **_: (i3e, i1g), -1, None),
    ),
    transitions_after=(
        (lambda iX, iRa, **_: (iX, iRa), +1, None),
        (lambda iRb, iX, **_: (iRb, iX), +1, 3),
    ),
)


def _generate_relaxation_pathway(
    self: AggregateSpectroscopy,
    desc: _RelaxationPathwayDesc,
    lst: list,
    eUt2: numpy.ndarray,
    pop_tol: float,
    dip_tol: float,
    evf_tol: float,
    verbose: int = 0,
) -> None:
    ngs = self.get_electronic_groundstate()
    nes = self.get_excitonic_band(band=1)

    if desc.requires_band2:
        try:
            nfs = self.get_excitonic_band(band=2)
        except Exception:
            raise Exception(
                f"Band-2 states not available for {desc.pname} pathway generation"
            )

    relax_set = nes if desc.relax_band == 1 else ngs
    if desc.inner_band == 0:
        inner_set = ngs
    elif desc.inner_band == 1:
        inner_set = nes
    else:
        inner_set = nfs  # type: ignore[possibly-undefined]

    if verbose > 0:
        print(f"Liouville pathway {desc.pname}")
        print("Population tolerance: ", pop_tol)
        print("Dipole tolerance:     ", dip_tol)
        print("Evolution amplitude:  ", evf_tol)

    for i1g in ngs:
        if verbose > 0:
            print("Ground state: ", i1g, "of", len(ngs))

        if self.rho0[i1g, i1g] > pop_tol:
            for i2e in nes:
                if verbose > 1:
                    print("Excited state: ", i2e, "of", len(nes))

                if self.D2[i2e, i1g] > dip_tol:
                    for i3e in nes:
                        if self.D2[i3e, i1g] > dip_tol:
                            for iRa in relax_set:
                                for iRb in relax_set:
                                    evf_idx = desc.evf_fn(
                                        iRa=iRa, iRb=iRb, i2e=i2e, i3e=i3e
                                    )
                                    evf = eUt2[evf_idx]

                                    if abs(evf) > evf_tol:
                                        for iX in inner_set:
                                            if desc.inner_dip_check(
                                                self,
                                                iX,
                                                iRa,
                                                iRb,
                                                i2e,
                                                i3e,
                                                dip_tol,
                                            ):
                                                try:
                                                    if verbose > 5:
                                                        print(
                                                            f" * Generating {desc.pname}",
                                                            i1g,
                                                            i2e,
                                                            i3e,
                                                        )

                                                    lp = diag.liouville_pathway(
                                                        desc.rtype,
                                                        i1g,
                                                        aggregate=self,
                                                        order=3,
                                                        pname=desc.pname,
                                                        popt_band=1,
                                                        relax_order=1,
                                                    )

                                                    w1_pair = desc.width1_fn(
                                                        i2e=i2e,
                                                        i1g=i1g,
                                                        iRa=iRa,
                                                        iRb=iRb,
                                                        iX=iX,
                                                    )
                                                    width1 = self.get_transition_width(
                                                        w1_pair
                                                    )
                                                    deph1 = (
                                                        self.get_transition_dephasing(
                                                            w1_pair
                                                        )
                                                    )

                                                    w3_pair = desc.width3_fn(
                                                        i2e=i2e,
                                                        i1g=i1g,
                                                        iRa=iRa,
                                                        iRb=iRb,
                                                        iX=iX,
                                                    )
                                                    width3 = self.get_transition_width(
                                                        w3_pair
                                                    )
                                                    deph3 = (
                                                        self.get_transition_dephasing(
                                                            w3_pair
                                                        )
                                                    )

                                                    kw = dict(
                                                        i1g=i1g,
                                                        i2e=i2e,
                                                        i3e=i3e,
                                                        iRa=iRa,
                                                        iRb=iRb,
                                                        iX=iX,
                                                    )
                                                    for (
                                                        idx_fn,
                                                        side,
                                                        interval,
                                                    ) in desc.transitions_before:
                                                        pair = idx_fn(**kw)
                                                        if interval == 1:
                                                            lp.add_transition(
                                                                pair,
                                                                side,
                                                                interval=1,
                                                                width=width1,
                                                                deph=deph1,
                                                            )
                                                        else:
                                                            lp.add_transition(
                                                                pair, side
                                                            )

                                                    transfer_target, transfer_source = (
                                                        desc.transfer_fn(
                                                            iRa=iRa,
                                                            iRb=iRb,
                                                            i2e=i2e,
                                                            i3e=i3e,
                                                        )
                                                    )
                                                    lp.add_transfer(
                                                        transfer_target,
                                                        transfer_source,
                                                    )
                                                    lp.set_evolution_factor(evf)

                                                    for (
                                                        idx_fn,
                                                        side,
                                                        interval,
                                                    ) in desc.transitions_after:
                                                        pair = idx_fn(**kw)
                                                        if interval == 3:
                                                            lp.add_transition(
                                                                pair,
                                                                side,
                                                                interval=3,
                                                                width=width3,
                                                                deph=deph3,
                                                            )
                                                        else:
                                                            lp.add_transition(
                                                                pair, side
                                                            )

                                                except Exception:
                                                    raise Exception(
                                                        f"Generation of {desc.pname}"
                                                        " pathway failed"
                                                    )

                                                lp.build()
                                                lst.append(lp)


def _generate_ground_excited_pathway(
    self: AggregateSpectroscopy,
    desc: _GroundExcitedPathwayDesc,
    lst: list,
    eUt2: numpy.ndarray,
    pop_tol: float,
    dip_tol: float,
    verbose: int = 0,
) -> None:
    ngs = self.get_electronic_groundstate()
    nes = self.get_excitonic_band(band=1)

    if verbose > 0:
        print(f"Liouville pathway {desc.pname}")
        print("Population tolerance: ", pop_tol)
        print("Dipole tolerance:     ", dip_tol)

    for i1g in ngs:
        if verbose > 0:
            print("Ground state: ", i1g, "of", len(ngs))

        if self.rho0[i1g, i1g] > pop_tol:
            for i2e in nes:
                if verbose > 1:
                    print("Excited state: ", i2e, "of", len(nes))

                if self.D2[i2e, i1g] > dip_tol:
                    for i3g in ngs:
                        if self.D2[i3g, i2e] > dip_tol:
                            evf_idx = desc.evf_fn(i1g, i3g)
                            evf = eUt2[evf_idx]

                            for i4e in nes:
                                if desc.inner_dip_check(self, i4e, i1g, i3g, dip_tol):
                                    try:
                                        if verbose > 5:
                                            print(
                                                f" * Generating {desc.pname}",
                                                i1g,
                                                i2e,
                                            )

                                        lp = diag.liouville_pathway(
                                            desc.rtype,
                                            i1g,
                                            aggregate=self,
                                            order=3,
                                            pname=desc.pname,
                                        )

                                        w1_pair = desc.width1_fn(
                                            i2e=i2e, i1g=i1g, i3g=i3g, i4e=i4e
                                        )
                                        width1 = self.get_transition_width(w1_pair)
                                        deph1 = self.get_transition_dephasing(w1_pair)

                                        w3_pair = desc.width3_fn(
                                            i2e=i2e, i1g=i1g, i3g=i3g, i4e=i4e
                                        )
                                        width3 = self.get_transition_width(w3_pair)
                                        deph3 = self.get_transition_dephasing(w3_pair)

                                        kw = dict(i1g=i1g, i2e=i2e, i3g=i3g, i4e=i4e)
                                        for i, (idx_fn, side, interval) in enumerate(
                                            desc.transitions
                                        ):
                                            pair = idx_fn(**kw)
                                            kwargs: dict[str, Any] = {}
                                            if interval == 1:
                                                kwargs.update(
                                                    interval=1,
                                                    width=width1,
                                                    deph=deph1,
                                                )
                                            elif interval == 3:
                                                kwargs.update(
                                                    interval=3,
                                                    width=width3,
                                                    deph=deph3,
                                                )
                                            lp.add_transition(pair, side, **kwargs)

                                        lp.set_evolution_factor(evf)

                                    except Exception:
                                        raise Exception(
                                            f"Generation of {desc.pname} pathway failed"
                                        )

                                    lp.build()
                                    lst.append(lp)


def generate_R1g(
    self: AggregateSpectroscopy,
    lst: list,
    eUt2: numpy.ndarray,
    pop_tol: float,
    dip_tol: float,
    evf_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self, _R1G_DESC, lst, eUt2, pop_tol, dip_tol, evf_tol, verbose
    )


def generate_R1gE(
    self: AggregateSpectroscopy,
    lst: list,
    eUt2: numpy.ndarray,
    pop_tol: float,
    dip_tol: float,
    evf_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self, _R1GE_DESC, lst, eUt2, pop_tol, dip_tol, evf_tol, verbose
    )


def generate_R2g(
    self: AggregateSpectroscopy,
    lst: list,
    eUt2: numpy.ndarray,
    pop_tol: float,
    dip_tol: float,
    evf_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self, _R2G_DESC, lst, eUt2, pop_tol, dip_tol, evf_tol, verbose
    )


def generate_R2gE(
    self: AggregateSpectroscopy,
    lst: list,
    eUt2: numpy.ndarray,
    pop_tol: float,
    dip_tol: float,
    evf_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self, _R2GE_DESC, lst, eUt2, pop_tol, dip_tol, evf_tol, verbose
    )


def generate_R3g(
    self: AggregateSpectroscopy,
    lst: list,
    eUt2: numpy.ndarray,
    pop_tol: float,
    dip_tol: float,
    verbose: int = 0,
) -> None:
    _generate_ground_excited_pathway(
        self, _R3G_DESC, lst, eUt2, pop_tol, dip_tol, verbose
    )


def generate_R4g(
    self: AggregateSpectroscopy,
    lst: list,
    eUt2: numpy.ndarray,
    pop_tol: float,
    dip_tol: float,
    verbose: int = 0,
) -> None:
    _generate_ground_excited_pathway(
        self, _R4G_DESC, lst, eUt2, pop_tol, dip_tol, verbose
    )


def generate_R1f(
    self: AggregateSpectroscopy,
    lst: list,
    eUt2: numpy.ndarray,
    pop_tol: float,
    dip_tol: float,
    evf_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self, _R1F_DESC, lst, eUt2, pop_tol, dip_tol, evf_tol, verbose
    )


def generate_R2f(
    self: AggregateSpectroscopy,
    lst: list,
    eUt2: numpy.ndarray,
    pop_tol: float,
    dip_tol: float,
    evf_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self, _R2F_DESC, lst, eUt2, pop_tol, dip_tol, evf_tol, verbose
    )


def generate_R1fE(
    self: AggregateSpectroscopy,
    lst: list,
    eUt2: numpy.ndarray,
    pop_tol: float,
    dip_tol: float,
    evf_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self, _R1FE_DESC, lst, eUt2, pop_tol, dip_tol, evf_tol, verbose
    )


def generate_R2fE(
    self: AggregateSpectroscopy,
    lst: list,
    eUt2: numpy.ndarray,
    pop_tol: float,
    dip_tol: float,
    evf_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self, _R2FE_DESC, lst, eUt2, pop_tol, dip_tol, evf_tol, verbose
    )


def generate_1orderP_sec(
    self: Any, lst: list, pop_tol: float, dip_tol: float, verbose: int = 0
) -> None:

    ngs = self.get_electronic_groundstate()
    nes = self.get_excitonic_band(band=1)

    if verbose > 0:
        print("Liouville pathway of first order")
        print("Population tolerance: ", pop_tol)
        print("Dipole tolerance:     ", dip_tol)

    k = 0
    l = 0
    for i1g in ngs:
        if verbose > 0:
            print("Ground state: ", i1g, "of", len(ngs))

        # Only thermally allowed starting states are considered
        if self.rho0[i1g, i1g] > pop_tol:
            for i2e in nes:
                if self.D2[i2e, i1g] > dip_tol:
                    l += 1

                    #      Diagram P1
                    #
                    #
                    #      |g_i1> <g_i1|
                    # <----|-----------|
                    #      |e_i2> <g_i1|
                    # ---->|-----------|
                    #      |g_i1> <g_i1|

                    try:
                        if verbose > 5:
                            print(" * Generating P1", i1g, i2e)

                        lp = diag.liouville_pathway(
                            "NR",
                            i1g,
                            aggregate=self,
                            order=1,
                            pname="P1",
                            popt_band=1,
                            relax_order=1,
                        )

                        # first transition lineshape
                        width1 = self.get_transition_width((i2e, i1g))
                        deph1 = self.get_transition_dephasing((i2e, i1g))

                        #      |g_i1> <g_i1|
                        lp.add_transition(
                            (i2e, i1g), +1, interval=1, width=width1, deph=deph1
                        )
                        #      |e_i2> <g_i1|
                        lp.add_transition(
                            (i1g, i2e), +1, interval=1, width=width1, deph=deph1
                        )
                        #      |g_i1> <g_i1|

                    except Exception:
                        break

                    lp.build()
                    lst.append(lp)
                    k += 1
