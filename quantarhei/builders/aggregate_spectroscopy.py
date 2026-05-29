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
        #
        # Rest is ignored for now (may be valuabel in the future)
        #

        population_tol = ptol
        dipole_tol = numpy.sqrt(self.D2_max) * dtol

        # Check if the ptype is a tuple
        ptype_tuple: Any
        if not isinstance(ptype, (tuple, list)):
            ptype_tuple = (ptype,)
        else:
            ptype_tuple = ptype
        pathways: list[Any] = []

        for ptp in ptype_tuple:
            if ptp == "R3g":
                ground_states = self.get_electronic_groundstate()
                excited_states = self.get_excitonic_band(band=1)

                k = 0
                l = 0
                for ground in ground_states:
                    # Only thermally allowed starting states are considered
                    if self.rho0[ground, ground] > population_tol:
                        for exc_ket in excited_states:
                            if self.D2[exc_ket, ground] < dipole_tol:
                                break

                            for mid_ground in ground_states:
                                if self.D2[mid_ground, exc_ket] < dipole_tol:
                                    break

                                for final_exc in excited_states:
                                    if (self.D2[final_exc, ground] < dipole_tol) and (
                                        self.D2[mid_ground, final_exc] < dipole_tol
                                    ):
                                        break

                                    l += 1

                                    #      Diagram R3g
                                    #
                                    #
                                    #      |g_i3> <g_i3|
                                    # <----|-----------|
                                    #      |e_i4> <g_i3|
                                    # ---->|-----------|
                                    #      |g_i1> <g_i3|
                                    #      |-----------|---->
                                    #      |g_i1> <e_i2|
                                    #      |-----------|<----
                                    #      |g_i1> <g_i1|

                                    try:
                                        pathway = diag.liouville_pathway(
                                            "R",
                                            ground,
                                            aggregate=self,
                                            order=3,
                                            pname=ptp,
                                        )
                                        # |g_i1> <g_i1|
                                        pathway.add_transition((exc_ket, ground), -1)
                                        # |g_i1> <e_i2|
                                        pathway.add_transition(
                                            (mid_ground, exc_ket), -1
                                        )
                                        # |g_i1> <g_i3|
                                        pathway.add_transition((final_exc, ground), +1)
                                        # |e_i5> <g_i3|
                                        pathway.add_transition(
                                            (mid_ground, final_exc), +1
                                        )
                                        # |g_i3> <g_i3|

                                    except Exception:
                                        break

                                    pathway.build()
                                    pathways.append(pathway)
                                    k += 1

            if ptp == "R2g":
                ground_states = self.get_electronic_groundstate()
                excited_states = self.get_excitonic_band(band=1)

                k = 0
                l = 0
                for ground in ground_states:
                    # Only thermally allowed starting states are considered
                    if self.rho0[ground, ground] > population_tol:
                        for exc_ket in excited_states:
                            if self.D2[exc_ket, ground] < dipole_tol:
                                break

                            for exc_bra in excited_states:
                                if self.D2[exc_bra, ground] < dipole_tol:
                                    break

                                for final_ground in ground_states:
                                    if (
                                        self.D2[final_ground, exc_ket] < dipole_tol
                                    ) or (self.D2[final_ground, exc_bra] < dipole_tol):
                                        break

                                    l += 1

                                    #      Diagram R2g
                                    #
                                    #
                                    #      |g_i4> <g_i4|
                                    # <----|-----------|
                                    #      |e_i3> <g_i4|
                                    #      |-----------|---->
                                    #      |e_i3> <e_i2|
                                    # ---->|-----------|
                                    #      |g_i1> <e_i2|
                                    #      |-----------|<----
                                    #      |g_i1> <g_i1|

                                    try:
                                        pathway = diag.liouville_pathway(
                                            "R",
                                            ground,
                                            aggregate=self,
                                            order=3,
                                            pname=ptp,
                                            popt_band=1,
                                        )
                                        #      |g_i1> <g_i1|
                                        pathway.add_transition((exc_ket, ground), -1)
                                        #      |g_i1> <e_i2|
                                        pathway.add_transition((exc_bra, ground), +1)
                                        #      |e_i3> <e_i2|
                                        pathway.add_transition(
                                            (final_ground, exc_ket), -1
                                        )
                                        #      |e_i3> <g_i4|
                                        pathway.add_transition(
                                            (final_ground, exc_bra), +1
                                        )
                                        #      |g_i4> <g_i4|

                                    except Exception:
                                        break

                                    pathway.build()
                                    pathways.append(pathway)
                                    k += 1

            if ptp == "R1g":
                ground_states = self.get_electronic_groundstate()
                excited_states = self.get_excitonic_band(band=1)

                # nrg = len(ground_states)
                # nre = len(excited_states)

                # print("Ground state : ", nrg)
                # print("Excited state: ", nre)
                # print("R1g: ",nrg*nre*nre*nrg)

                k = 0
                l = 0
                for ground in ground_states:
                    # Only thermally allowed starting states are considered
                    if self.rho0[ground, ground] > population_tol:
                        for exc_ket in excited_states:
                            if self.D2[exc_ket, ground] < dipole_tol:
                                break

                            for exc_bra in excited_states:
                                if self.D2[exc_bra, ground] < dipole_tol:
                                    break

                                for final_ground in ground_states:
                                    if (
                                        self.D2[final_ground, exc_bra] < dipole_tol
                                    ) or (self.D2[final_ground, exc_ket] < dipole_tol):
                                        break

                                    l += 1

                                    #      Diagram R1g
                                    #
                                    #
                                    #      |g_i4> <g_i4|
                                    # <----|-----------|
                                    #      |e_i2> <g_i4|
                                    #      |-----------|---->
                                    #      |e_i2> <e_i3|
                                    #      |-----------|<----
                                    #      |e_i2> <g_i1|
                                    # ---->|-----------|
                                    #      |g_i1> <g_i1|

                                    try:
                                        pathway = diag.liouville_pathway(
                                            "NR",
                                            ground,
                                            aggregate=self,
                                            order=3,
                                            pname=ptp,
                                            popt_band=1,
                                        )
                                        #      |g_i1> <g_i1|
                                        pathway.add_transition((exc_ket, ground), +1)
                                        #      |e_i2> <g_i1|
                                        pathway.add_transition((exc_bra, ground), -1)
                                        #      |e_i2> <e_i3|
                                        pathway.add_transition(
                                            (final_ground, exc_bra), -1
                                        )
                                        #      |e_i2> <g_i4|
                                        pathway.add_transition(
                                            (final_ground, exc_ket), +1
                                        )
                                        #      |g_i4> <g_i4|

                                    except Exception:
                                        break

                                    pathway.build()
                                    pathways.append(pathway)
                                    k += 1

            if ptp == "R4g":
                ground_states = self.get_electronic_groundstate()
                excited_states = self.get_excitonic_band(band=1)

                # nrg = len(ground_states)
                # nre = len(excited_states)

                # print("Ground state : ", nrg)
                # print("Excited state: ", nre)
                # print("R4g: ",nrg*nre*nrg*nrg*nre)

                k = 0
                l = 0
                for ground in ground_states:
                    # Only thermally allowed starting states are considered
                    if self.rho0[ground, ground] > population_tol:
                        for exc_ket in excited_states:
                            if self.D2[exc_ket, ground] < dipole_tol:
                                break

                            for mid_ground in ground_states:
                                if self.D2[mid_ground, exc_ket] < dipole_tol:
                                    break

                                for final_exc in excited_states:
                                    if (
                                        self.D2[final_exc, mid_ground] < dipole_tol
                                    ) or (self.D2[ground, final_exc] < dipole_tol):
                                        break

                                    l += 1

                                    #      Diagram R4g
                                    #
                                    #
                                    #      |g_i1> <g_i1|
                                    # <----|-----------|
                                    #      |e_i4> <g_i1|
                                    # ---->|-----------|
                                    #      |g_i3> <g_i1|
                                    # <----|-----------|
                                    #      |e_i2> <g_i1|
                                    # ---->|-----------|
                                    #      |g_i1> <g_i1|

                                    try:
                                        pathway = diag.liouville_pathway(
                                            "NR",
                                            ground,
                                            aggregate=self,
                                            order=3,
                                            pname=ptp,
                                        )
                                        #      |g_i1> <g_i1|
                                        pathway.add_transition((exc_ket, ground), +1)
                                        #      |e_i2> <g_i1|
                                        pathway.add_transition(
                                            (mid_ground, exc_ket), +1
                                        )
                                        #      |g_i3> <g_i1|
                                        pathway.add_transition(
                                            (final_exc, mid_ground), +1
                                        )
                                        #      |e_i4> <g_i1|
                                        pathway.add_transition((ground, final_exc), +1)
                                        #      |g_i1> <g_i1|

                                    except Exception:
                                        break

                                    pathway.build()
                                    pathways.append(pathway)
                                    k += 1

            if ptp == "R1f*":
                ground_states = self.get_electronic_groundstate()
                excited_states = self.get_excitonic_band(band=1)
                try:
                    two_exciton_states = self.get_excitonic_band(band=2)
                except Exception:
                    break

                #                print(ground_states)
                #                print(excited_states)
                #                print(two_exciton_states)
                #                for a in excited_states:
                #                    for b in two_exciton_states:
                #                        print(a,b," : ",self.D2[a,b],self.D2[b,a])

                k = 0
                l = 0
                for ground in ground_states:
                    # Only thermally allowed starting states are considered
                    if self.rho0[ground, ground] > population_tol:
                        for exc_ket in excited_states:
                            if self.D2[exc_ket, ground] < dipole_tol:
                                break

                            for exc_bra in excited_states:
                                if self.D2[exc_bra, ground] < dipole_tol:
                                    break

                                for two_exc in two_exciton_states:
                                    if (self.D2[two_exc, exc_bra] < dipole_tol) or (
                                        self.D2[exc_ket, two_exc] < dipole_tol
                                    ):
                                        # print("Breaking")
                                        # print(self.D2[two_exc,exc_bra],self.D2[exc_ket,two_exc])
                                        break

                                    l += 1

                                    #      Diagram R4g
                                    #
                                    #
                                    #      |e_i2> <e_i2|
                                    # <----|-----------|
                                    #      |f_i4> <e_i2|
                                    # ---->|-----------|
                                    #      |e_i3> <e_i2|
                                    # ---->|-----------|
                                    #      |g_i1> <e_i2|
                                    #      |-----------|<----
                                    #      |g_i1> <g_i1|

                                    try:
                                        pathway = diag.liouville_pathway(
                                            "R",
                                            ground,
                                            aggregate=self,
                                            order=3,
                                            pname=ptp,
                                            popt_band=1,
                                        )
                                        #      |g_i1> <g_i1|
                                        pathway.add_transition((exc_ket, ground), -1)
                                        #      |g_i1> <e_i2|
                                        pathway.add_transition((exc_bra, ground), +1)
                                        #      |e_i3> <e_i2|
                                        pathway.add_transition((two_exc, exc_bra), +1)
                                        #      |f_i4> <e_i2|
                                        pathway.add_transition((exc_ket, two_exc), +1)
                                        #      |e_i2> <e_i2|

                                    except Exception:
                                        break

                                    pathway.build()
                                    pathways.append(pathway)
                                    k += 1

            if ptp == "R2f*":
                ground_states = self.get_electronic_groundstate()
                excited_states = self.get_excitonic_band(band=1)

                try:
                    two_exciton_states = self.get_excitonic_band(band=2)
                except Exception:
                    break

                k = 0
                l = 0
                for ground in ground_states:
                    # Only thermally allowed starting states are considered
                    if self.rho0[ground, ground] > population_tol:
                        for exc_ket in excited_states:
                            if self.D2[exc_ket, ground] < dipole_tol:
                                break

                            for exc_bra in excited_states:
                                if self.D2[exc_bra, ground] < dipole_tol:
                                    break

                                for two_exc in two_exciton_states:
                                    if (self.D2[two_exc, exc_ket] < dipole_tol) or (
                                        self.D2[exc_bra, two_exc] < dipole_tol
                                    ):
                                        break

                                    l += 1

                                    #      Diagram R4g
                                    #
                                    #
                                    #      |e_i3> <e_i3|
                                    # <----|-----------|
                                    #      |f_i4> <e_i3|
                                    # ---->|-----------|
                                    #      |e_i2> <e_i3|
                                    #      |-----------|<----
                                    #      |e_i2> <g_i1|
                                    # ---->|-----------|
                                    #      |g_i1> <g_i1|

                                    try:
                                        pathway = diag.liouville_pathway(
                                            "NR",
                                            ground,
                                            aggregate=self,
                                            order=3,
                                            pname=ptp,
                                            popt_band=1,
                                        )
                                        #      |g_i1> <g_i1|
                                        pathway.add_transition((exc_ket, ground), +1)
                                        #      |e_i2> <g_i1|
                                        pathway.add_transition((exc_bra, ground), -1)
                                        #      |e_i2> <e_i3|
                                        pathway.add_transition((two_exc, exc_ket), +1)
                                        #      |f_i4> <e_i3|
                                        pathway.add_transition((exc_bra, two_exc), +1)
                                        #      |e_i3> <e_i3|

                                    except Exception:
                                        break

                                    pathway.build()
                                    pathways.append(pathway)
                                    k += 1

            if ptp == "R2g->3g":
                ground_states = self.get_electronic_groundstate()
                excited_states = self.get_excitonic_band(band=1)

                k = 0
                l = 0
                for ground in ground_states:
                    # Only thermally allowed starting states are considered
                    if self.rho0[ground, ground] > population_tol:
                        for exc_ket in excited_states:
                            if self.D2[exc_ket, ground] < dipole_tol:
                                break

                            for exc_bra in excited_states:
                                if self.D2[exc_bra, ground] < dipole_tol:
                                    break

                                # relaxation
                                for final_ground in ground_states:
                                    for mid_ground_b in ground_states:
                                        for detect_exc in excited_states:
                                            if (
                                                self.D2[detect_exc, final_ground]
                                                < dipole_tol
                                            ) or (
                                                self.D2[mid_ground_b, detect_exc]
                                                < dipole_tol
                                            ):
                                                break

                                            l += 1

                                            #      Diagram R2g_ETICS
                                            #      (Compensates R3g)
                                            #
                                            #
                                            #      |g_i5> <g_i5|
                                            # <----|-----------|
                                            #      |e_i6> <g_i5|
                                            # ---->|-----------|
                                            #      |g_i4> <g_i5|
                                            #      |***********|
                                            #      |e_i3> <e_i2|
                                            # ---->|-----------|
                                            #      |g_i1> <e_i2|
                                            #      |-----------|<----
                                            #      |g_i1> <g_i1|

                                            pathway = diag.liouville_pathway(
                                                "R_E",
                                                ground,
                                                aggregate=self,
                                                order=3,
                                                relax_order=1,
                                                pname=ptp,
                                            )
                                            #      |g_i1> <g_i1|
                                            pathway.add_transition(
                                                (exc_ket, ground), -1
                                            )
                                            #      |g_i1> <e_i2|
                                            pathway.add_transition(
                                                (exc_bra, ground), +1
                                            )
                                            #      |e_i3> <e_i2|
                                            pathway.add_transfer(
                                                (final_ground, mid_ground_b),
                                                (exc_bra, exc_ket),
                                            )
                                            #      |g_i4> <g_i5|
                                            pathway.add_transition(
                                                (detect_exc, final_ground), +1
                                            )
                                            #      |e_i6> <g_i5|
                                            pathway.add_transition(
                                                (mid_ground_b, detect_exc), +1
                                            )
                                            #      |g_i5> <g_i5|

                                            pathway.build()
                                            pathways.append(pathway)
                                            k += 1

            if ptp == "R1g->4g":
                ground_states = self.get_electronic_groundstate()
                excited_states = self.get_excitonic_band(band=1)

                k = 0
                l = 0
                for ground in ground_states:
                    # Only thermally allowed starting states are considered
                    if self.rho0[ground, ground] > population_tol:
                        for exc_ket in excited_states:
                            if self.D2[exc_ket, ground] < dipole_tol:
                                break

                            for exc_bra in excited_states:
                                if self.D2[exc_bra, ground] < dipole_tol:
                                    break

                                # relaxation
                                for final_ground in ground_states:
                                    for mid_ground_b in ground_states:
                                        for detect_exc in excited_states:
                                            if (
                                                self.D2[detect_exc, final_ground]
                                                < dipole_tol
                                            ) or (
                                                self.D2[mid_ground_b, detect_exc]
                                                < dipole_tol
                                            ):
                                                break

                                            l += 1

                                            #      Diagram R2g_ETICS
                                            #      (Compensates R3g)
                                            #
                                            #
                                            #      |g_i5> <g_i5|
                                            # <----|-----------|
                                            #      |e_i6> <g_i5|
                                            # ---->|-----------|
                                            #      |g_i4> <g_i5|
                                            #      |***********|
                                            #      |e_i2> <e_i3|
                                            #      |-----------|<----
                                            #      |e_i2> <g_i1|
                                            # ---->|-----------|
                                            #      |g_i1> <g_i1|

                                            # if True:
                                            try:
                                                pathway = diag.liouville_pathway(
                                                    "NR_E",
                                                    ground,
                                                    aggregate=self,
                                                    order=3,
                                                    relax_order=1,
                                                    pname=ptp,
                                                )
                                                #      |g_i1> <g_i1|
                                                pathway.add_transition(
                                                    (exc_ket, ground), +1
                                                )
                                                #      |e_i2> <g_i1|
                                                pathway.add_transition(
                                                    (exc_bra, ground), -1
                                                )
                                                #      |e_i2> <e_i3|
                                                pathway.add_transfer(
                                                    (final_ground, mid_ground_b),
                                                    (exc_ket, exc_bra),
                                                )
                                                #      |g_i4> <g_i5|
                                                pathway.add_transition(
                                                    (detect_exc, final_ground), +1
                                                )
                                                #      |e_i6> <g_i5|
                                                pathway.add_transition(
                                                    (mid_ground_b, detect_exc), +1
                                                )
                                                #      |g_i5> <g_i5|

                                            except Exception:
                                                break

                                            pathway.build()
                                            pathways.append(pathway)
                                            k += 1

        if lab is not None:
            for l in pathways:
                l.orientational_averaging(lab)

        return pathways

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
        pathways : list
            List of LiouvillePathway objects


        """
        self.lab = lab

        if self._diagonalized:
            if verbose > 0:
                print("Diagonalizing aggregate")
            self.diagonalize()
            if verbose > 0:
                print("..done")

        population_tol = ptol
        dipole_tol = numpy.sqrt(self.D2_max) * dtol
        evolution_tol = etol

        # Check if the ptype is a tuple
        ptype_tuple: Any
        if not isinstance(ptype, (tuple, list)):
            ptype_tuple = (ptype,)
        else:
            ptype_tuple = ptype
        pathways: list[Any] = []

        if verbose > 0:
            print("Pathways", ptype_tuple)

        #
        # data of the evolution superoperator in eigenstate basis
        #

        try:
            # either the eUt is a complete evolution superoperator
            evolution_superop = eUt.at(t2)

            #
            # TODO: Check this part I had before without eigenbasis of and used the evolution superoperator values directly
            # ----------------------------------------------
            eUt2_dat = numpy.zeros(
                evolution_superop.data.shape, dtype=evolution_superop.data.dtype
            )  #
            HH = eUt.get_Hamiltonian()  #
            with eigenbasis_of(HH):  #
                eUt2_dat[:, :, :, :] = evolution_superop.data  #
        # ----------------------------------------------

        except AttributeError:
            # or it is only a super operator at a given time t2
            # in this case 'ham' must be specified
            evolution_superop = eUt
            eUt2_dat = numpy.zeros(
                evolution_superop.data.shape, dtype=evolution_superop.data.dtype
            )
            with eigenbasis_of(ham):
                eUt2_dat[:, :, :, :] = evolution_superop.data

        for ptp in ptype_tuple:
            if ptp == "R1g":
                generate_R1g(
                    self,
                    pathways,
                    eUt2_dat,
                    population_tol,
                    dipole_tol,
                    evolution_tol,
                    verbose,
                )

            elif ptp == "R2g":
                generate_R2g(
                    self,
                    pathways,
                    eUt2_dat,
                    population_tol,
                    dipole_tol,
                    evolution_tol,
                    verbose,
                )

            elif ptp == "R3g":
                generate_R3g(
                    self, pathways, eUt2_dat, population_tol, dipole_tol, verbose
                )

            elif ptp == "R4g":
                generate_R4g(
                    self, pathways, eUt2_dat, population_tol, dipole_tol, verbose
                )

            elif ptp == "R1f*":
                generate_R1f(
                    self,
                    pathways,
                    eUt2_dat,
                    population_tol,
                    dipole_tol,
                    evolution_tol,
                    verbose,
                )

            elif ptp == "R2f*":
                generate_R2f(
                    self,
                    pathways,
                    eUt2_dat,
                    population_tol,
                    dipole_tol,
                    evolution_tol,
                    verbose,
                )

            elif ptp == "R1gE":
                generate_R1gE(
                    self,
                    pathways,
                    eUt2_dat,
                    population_tol,
                    dipole_tol,
                    evolution_tol,
                    verbose,
                )

            elif ptp == "R2gE":
                generate_R2gE(
                    self,
                    pathways,
                    eUt2_dat,
                    population_tol,
                    dipole_tol,
                    evolution_tol,
                    verbose,
                )

            elif ptp == "R1f*E":
                generate_R1fE(
                    self,
                    pathways,
                    eUt2_dat,
                    population_tol,
                    dipole_tol,
                    evolution_tol,
                    verbose,
                )

            elif ptp == "R2f*E":
                generate_R2fE(
                    self,
                    pathways,
                    eUt2_dat,
                    population_tol,
                    dipole_tol,
                    evolution_tol,
                    verbose,
                )

            else:
                raise QuantarheiError("Unknown pythway type: " + str(ptp))

        if lab is not None:
            for l in pathways:
                l.orientational_averaging(lab)

        return pathways

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
        pathways : list
            List of LiouvillePathway objects


        """
        if self._diagonalized:
            if verbose > 0:
                print("Diagonalizing aggregate")
            self.diagonalize()
            if verbose > 0:
                print("..done")

        population_tol = ptol
        dipole_tol = numpy.sqrt(self.D2_max) * dtol
        evolution_tol = etol

        if eUt is None:
            # secular absorption spectrum calculation
            eUt2_dat = None
            sec = True

        else:
            raise ImplementationError("Not implemented yet")

        pathways: list[Any] = []

        if sec:
            generate_1orderP_sec(self, pathways, population_tol, dipole_tol, verbose)
        else:
            raise ImplementationError("Not implemented yet")

        if lab is not None:
            for l in pathways:
                l.orientational_averaging(lab)

        return pathways


# ---------------------------------------------------------------------------
# Data-driven pathway generation (incremental refactoring of #366)
# ---------------------------------------------------------------------------

from dataclasses import dataclass


@dataclass(frozen=True)
class _GroundExcitedPathwayDesc:
    pathway_type: str
    pathway_name: str
    evolution_factor_indices: Callable[[int, int], tuple[int, int, int, int]]
    detection_dipole_check: Callable[..., bool]
    first_transition_pair: Callable[..., tuple[int, int]]
    third_transition_pair: Callable[..., tuple[int, int]]
    transitions: tuple[tuple[Callable[..., tuple[int, int]], int, int | None], ...]


_R3G_DESC = _GroundExcitedPathwayDesc(
    pathway_type="R",
    pathway_name="R3g",
    evolution_factor_indices=lambda ground, mid_ground: (
        ground,
        mid_ground,
        ground,
        mid_ground,
    ),
    detection_dipole_check=lambda self, final_exc, ground, mid_ground, dipole_tol: (
        self.D2[final_exc, ground] > dipole_tol
        and self.D2[mid_ground, final_exc] > dipole_tol
    ),
    first_transition_pair=lambda exc_ket, ground, **_: (exc_ket, ground),
    third_transition_pair=lambda final_exc, mid_ground, **_: (final_exc, mid_ground),
    transitions=(
        # (pair_fn returning (a,b), side, interval)
        (lambda exc_ket, ground, **_: (exc_ket, ground), -1, 1),
        (lambda mid_ground, exc_ket, **_: (mid_ground, exc_ket), -1, None),
        (lambda final_exc, ground, **_: (final_exc, ground), +1, None),
        (lambda mid_ground, final_exc, **_: (mid_ground, final_exc), +1, 3),
    ),
)

_R4G_DESC = _GroundExcitedPathwayDesc(
    pathway_type="NR",
    pathway_name="R4g",
    evolution_factor_indices=lambda ground, mid_ground: (
        ground,
        mid_ground,
        ground,
        mid_ground,
    ),
    detection_dipole_check=lambda self, final_exc, ground, mid_ground, dipole_tol: (
        self.D2[final_exc, mid_ground] > dipole_tol
        and self.D2[ground, final_exc] > dipole_tol
    ),
    first_transition_pair=lambda exc_ket, ground, **_: (exc_ket, ground),
    third_transition_pair=lambda final_exc, ground, **_: (final_exc, ground),
    transitions=(
        (lambda exc_ket, ground, **_: (exc_ket, ground), +1, 1),
        (lambda mid_ground, exc_ket, **_: (mid_ground, exc_ket), +1, None),
        (lambda final_exc, mid_ground, **_: (final_exc, mid_ground), +1, None),
        (lambda ground, final_exc, **_: (ground, final_exc), +1, 3),
    ),
)


@dataclass(frozen=True)
class _RelaxationPathwayDesc:
    pathway_type: str
    pathway_name: str
    requires_two_exciton_band: bool
    relaxation_band: int  # 1 = excited_states, 0 = ground_states
    detection_band: int  # 0 = ground_states, 1 = excited_states, 2 = two_exciton_states
    evolution_factor_indices: Callable[..., tuple[int, int, int, int]]
    detection_dipole_check: Callable[..., bool]
    first_transition_pair: Callable[..., tuple[int, int]]
    third_transition_pair: Callable[..., tuple[int, int]]
    transfer_indices: Callable[..., tuple[tuple[int, int], tuple[int, int]]]
    transitions_before_transfer: tuple[
        tuple[Callable[..., tuple[int, int]], int, int | None], ...
    ]
    transitions_after_transfer: tuple[
        tuple[Callable[..., tuple[int, int]], int, int | None], ...
    ]


_R1G_DESC = _RelaxationPathwayDesc(
    pathway_type="NR",
    pathway_name="R1g",
    requires_two_exciton_band=False,
    relaxation_band=1,
    detection_band=0,
    evolution_factor_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        relax_a,
        relax_b,
        exc_ket,
        exc_bra,
    ),
    detection_dipole_check=lambda self, detect, relax_a, relax_b, exc_ket, exc_bra, dipole_tol: (
        self.D2[detect, relax_b] > dipole_tol and self.D2[detect, relax_a] > dipole_tol
    ),
    first_transition_pair=lambda exc_ket, ground, **_: (exc_ket, ground),
    third_transition_pair=lambda relax_a, detect, **_: (relax_a, detect),
    transfer_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        (relax_a, relax_b),
        (exc_ket, exc_bra),
    ),
    transitions_before_transfer=(
        (lambda exc_ket, ground, **_: (exc_ket, ground), +1, 1),
        (lambda exc_bra, ground, **_: (exc_bra, ground), -1, None),
    ),
    transitions_after_transfer=(
        (lambda detect, relax_b, **_: (detect, relax_b), -1, None),
        (lambda detect, relax_a, **_: (detect, relax_a), +1, 3),
    ),
)

_R2G_DESC = _RelaxationPathwayDesc(
    pathway_type="R",
    pathway_name="R2g",
    requires_two_exciton_band=False,
    relaxation_band=1,
    detection_band=0,
    evolution_factor_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        relax_a,
        relax_b,
        exc_bra,
        exc_ket,
    ),
    detection_dipole_check=lambda self, detect, relax_a, relax_b, exc_ket, exc_bra, dipole_tol: (
        self.D2[detect, exc_ket] > dipole_tol and self.D2[detect, exc_bra] > dipole_tol
    ),
    first_transition_pair=lambda exc_ket, ground, **_: (exc_ket, ground),
    third_transition_pair=lambda relax_a, detect, **_: (relax_a, detect),
    transfer_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        (relax_a, relax_b),
        (exc_bra, exc_ket),
    ),
    transitions_before_transfer=(
        (lambda exc_ket, ground, **_: (exc_ket, ground), -1, 1),
        (lambda exc_bra, ground, **_: (exc_bra, ground), +1, None),
    ),
    transitions_after_transfer=(
        (lambda detect, relax_b, **_: (detect, relax_b), -1, None),
        (lambda detect, relax_a, **_: (detect, relax_a), +1, 3),
    ),
)

_R1F_DESC = _RelaxationPathwayDesc(
    pathway_type="R",
    pathway_name="R1f*",
    requires_two_exciton_band=True,
    relaxation_band=1,
    detection_band=2,
    evolution_factor_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        relax_a,
        relax_b,
        exc_bra,
        exc_ket,
    ),
    detection_dipole_check=lambda self, detect, relax_a, relax_b, exc_ket, exc_bra, dipole_tol: (
        self.D2[detect, relax_a] > dipole_tol and self.D2[relax_b, detect] > dipole_tol
    ),
    first_transition_pair=lambda exc_ket, ground, **_: (exc_ket, ground),
    third_transition_pair=lambda detect, relax_b, **_: (detect, relax_b),
    transfer_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        (relax_a, relax_b),
        (exc_bra, exc_ket),
    ),
    transitions_before_transfer=(
        (lambda exc_ket, ground, **_: (exc_ket, ground), -1, 1),
        (lambda exc_bra, ground, **_: (exc_bra, ground), +1, None),
    ),
    transitions_after_transfer=(
        (lambda detect, relax_a, **_: (detect, relax_a), +1, None),
        (lambda relax_b, detect, **_: (relax_b, detect), +1, 3),
    ),
)

_R2F_DESC = _RelaxationPathwayDesc(
    pathway_type="NR",
    pathway_name="R2f*",
    requires_two_exciton_band=True,
    relaxation_band=1,
    detection_band=2,
    evolution_factor_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        relax_a,
        relax_b,
        exc_ket,
        exc_bra,
    ),
    detection_dipole_check=lambda self, detect, relax_a, relax_b, exc_ket, exc_bra, dipole_tol: (
        self.D2[detect, relax_a] > dipole_tol and self.D2[relax_b, detect] > dipole_tol
    ),
    first_transition_pair=lambda exc_ket, ground, **_: (exc_ket, ground),
    third_transition_pair=lambda detect, relax_b, **_: (detect, relax_b),
    transfer_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        (relax_a, relax_b),
        (exc_ket, exc_bra),
    ),
    transitions_before_transfer=(
        (lambda exc_ket, ground, **_: (exc_ket, ground), +1, 1),
        (lambda exc_bra, ground, **_: (exc_bra, ground), -1, None),
    ),
    transitions_after_transfer=(
        (lambda detect, relax_a, **_: (detect, relax_a), +1, None),
        (lambda relax_b, detect, **_: (relax_b, detect), +1, 3),
    ),
)

_R1GE_DESC = _RelaxationPathwayDesc(
    pathway_type="NR",
    pathway_name="R1gE",
    requires_two_exciton_band=False,
    relaxation_band=0,
    detection_band=1,
    evolution_factor_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        relax_a,
        relax_b,
        exc_ket,
        exc_bra,
    ),
    detection_dipole_check=lambda self, detect, relax_a, relax_b, exc_ket, exc_bra, dipole_tol: (
        self.D2[relax_a, detect] > dipole_tol and self.D2[relax_b, detect] > dipole_tol
    ),
    first_transition_pair=lambda exc_ket, ground, **_: (exc_ket, ground),
    third_transition_pair=lambda detect, relax_a, **_: (detect, relax_a),
    transfer_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        (relax_a, relax_b),
        (exc_ket, exc_bra),
    ),
    transitions_before_transfer=(
        (lambda exc_ket, ground, **_: (exc_ket, ground), +1, 1),
        (lambda exc_bra, ground, **_: (exc_bra, ground), -1, None),
    ),
    transitions_after_transfer=(
        (lambda detect, relax_a, **_: (detect, relax_a), +1, None),
        (lambda relax_b, detect, **_: (relax_b, detect), +1, 3),
    ),
)

_R2GE_DESC = _RelaxationPathwayDesc(
    pathway_type="R",
    pathway_name="R2gE",
    requires_two_exciton_band=False,
    relaxation_band=0,
    detection_band=1,
    evolution_factor_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        relax_a,
        relax_b,
        exc_bra,
        exc_ket,
    ),
    detection_dipole_check=lambda self, detect, relax_a, relax_b, exc_ket, exc_bra, dipole_tol: (
        self.D2[relax_a, detect] > dipole_tol and self.D2[relax_b, detect] > dipole_tol
    ),
    first_transition_pair=lambda exc_ket, ground, **_: (exc_ket, ground),
    third_transition_pair=lambda detect, relax_a, **_: (detect, relax_a),
    transfer_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        (relax_a, relax_b),
        (exc_bra, exc_ket),
    ),
    transitions_before_transfer=(
        (lambda exc_ket, ground, **_: (exc_ket, ground), -1, 1),
        (lambda exc_bra, ground, **_: (exc_bra, ground), +1, None),
    ),
    transitions_after_transfer=(
        (lambda detect, relax_a, **_: (detect, relax_a), +1, None),
        (lambda relax_b, detect, **_: (relax_b, detect), +1, 3),
    ),
)

_R1FE_DESC = _RelaxationPathwayDesc(
    pathway_type="R",
    pathway_name="R1f*E",
    requires_two_exciton_band=False,
    relaxation_band=0,
    detection_band=1,
    evolution_factor_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        relax_a,
        relax_b,
        exc_bra,
        exc_ket,
    ),
    detection_dipole_check=lambda self, detect, relax_a, relax_b, exc_ket, exc_bra, dipole_tol: (
        self.D2[detect, relax_a] > dipole_tol and self.D2[relax_b, detect] > dipole_tol
    ),
    first_transition_pair=lambda exc_ket, ground, **_: (exc_ket, ground),
    third_transition_pair=lambda detect, relax_b, **_: (detect, relax_b),
    transfer_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        (relax_a, relax_b),
        (exc_bra, exc_ket),
    ),
    transitions_before_transfer=(
        (lambda exc_ket, ground, **_: (exc_ket, ground), -1, 1),
        (lambda exc_bra, ground, **_: (exc_bra, ground), +1, None),
    ),
    transitions_after_transfer=(
        (lambda detect, relax_a, **_: (detect, relax_a), +1, None),
        (lambda relax_b, detect, **_: (relax_b, detect), +1, 3),
    ),
)

_R2FE_DESC = _RelaxationPathwayDesc(
    pathway_type="NR",
    pathway_name="R2f*E",
    requires_two_exciton_band=False,
    relaxation_band=0,
    detection_band=1,
    evolution_factor_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        relax_a,
        relax_b,
        exc_ket,
        exc_bra,
    ),
    detection_dipole_check=lambda self, detect, relax_a, relax_b, exc_ket, exc_bra, dipole_tol: (
        self.D2[detect, relax_a] > dipole_tol and self.D2[relax_b, detect] > dipole_tol
    ),
    first_transition_pair=lambda exc_ket, ground, **_: (exc_ket, ground),
    third_transition_pair=lambda detect, relax_b, **_: (detect, relax_b),
    transfer_indices=lambda relax_a, relax_b, exc_ket, exc_bra: (
        (relax_a, relax_b),
        (exc_ket, exc_bra),
    ),
    transitions_before_transfer=(
        (lambda exc_ket, ground, **_: (exc_ket, ground), +1, 1),
        (lambda exc_bra, ground, **_: (exc_bra, ground), -1, None),
    ),
    transitions_after_transfer=(
        (lambda detect, relax_a, **_: (detect, relax_a), +1, None),
        (lambda relax_b, detect, **_: (relax_b, detect), +1, 3),
    ),
)


def _generate_relaxation_pathway(
    self: AggregateSpectroscopy,
    desc: _RelaxationPathwayDesc,
    pathways: list,
    evolution_superop: numpy.ndarray,
    population_tol: float,
    dipole_tol: float,
    evolution_tol: float,
    verbose: int = 0,
) -> None:
    ground_states = self.get_electronic_groundstate()
    excited_states = self.get_excitonic_band(band=1)

    if desc.requires_two_exciton_band:
        try:
            two_exciton_states = self.get_excitonic_band(band=2)
        except Exception:
            raise Exception(
                f"Band-2 states not available for {desc.pathway_name} pathway generation"
            )

    relaxation_states = excited_states if desc.relaxation_band == 1 else ground_states
    if desc.detection_band == 0:
        detection_states = ground_states
    elif desc.detection_band == 1:
        detection_states = excited_states
    else:
        detection_states = two_exciton_states  # type: ignore[possibly-undefined]

    if verbose > 0:
        print(f"Liouville pathway {desc.pathway_name}")
        print("Population tolerance: ", population_tol)
        print("Dipole tolerance:     ", dipole_tol)
        print("Evolution amplitude:  ", evolution_tol)

    for ground in ground_states:
        if verbose > 0:
            print("Ground state: ", ground, "of", len(ground_states))

        if self.rho0[ground, ground] > population_tol:
            for exc_ket in excited_states:
                if verbose > 1:
                    print("Excited state: ", exc_ket, "of", len(excited_states))

                if self.D2[exc_ket, ground] > dipole_tol:
                    for exc_bra in excited_states:
                        if self.D2[exc_bra, ground] > dipole_tol:
                            for relax_a in relaxation_states:
                                for relax_b in relaxation_states:
                                    evolution_indices = desc.evolution_factor_indices(
                                        relax_a=relax_a,
                                        relax_b=relax_b,
                                        exc_ket=exc_ket,
                                        exc_bra=exc_bra,
                                    )
                                    evolution_factor = evolution_superop[
                                        evolution_indices
                                    ]

                                    if abs(evolution_factor) > evolution_tol:
                                        for detect in detection_states:
                                            if desc.detection_dipole_check(
                                                self,
                                                detect,
                                                relax_a,
                                                relax_b,
                                                exc_ket,
                                                exc_bra,
                                                dipole_tol,
                                            ):
                                                try:
                                                    if verbose > 5:
                                                        print(
                                                            f" * Generating {desc.pathway_name}",
                                                            ground,
                                                            exc_ket,
                                                            exc_bra,
                                                        )

                                                    pathway = diag.liouville_pathway(
                                                        desc.pathway_type,
                                                        ground,
                                                        aggregate=self,
                                                        order=3,
                                                        pname=desc.pathway_name,
                                                        popt_band=1,
                                                        relax_order=1,
                                                    )

                                                    first_lineshape_pair = (
                                                        desc.first_transition_pair(
                                                            exc_ket=exc_ket,
                                                            ground=ground,
                                                            relax_a=relax_a,
                                                            relax_b=relax_b,
                                                            detect=detect,
                                                        )
                                                    )
                                                    width1 = self.get_transition_width(
                                                        first_lineshape_pair
                                                    )
                                                    deph1 = (
                                                        self.get_transition_dephasing(
                                                            first_lineshape_pair
                                                        )
                                                    )

                                                    third_lineshape_pair = (
                                                        desc.third_transition_pair(
                                                            exc_ket=exc_ket,
                                                            ground=ground,
                                                            relax_a=relax_a,
                                                            relax_b=relax_b,
                                                            detect=detect,
                                                        )
                                                    )
                                                    width3 = self.get_transition_width(
                                                        third_lineshape_pair
                                                    )
                                                    deph3 = (
                                                        self.get_transition_dephasing(
                                                            third_lineshape_pair
                                                        )
                                                    )

                                                    state_indices = dict(
                                                        ground=ground,
                                                        exc_ket=exc_ket,
                                                        exc_bra=exc_bra,
                                                        relax_a=relax_a,
                                                        relax_b=relax_b,
                                                        detect=detect,
                                                    )
                                                    for (
                                                        pair_fn,
                                                        side,
                                                        interval,
                                                    ) in (
                                                        desc.transitions_before_transfer
                                                    ):
                                                        pair = pair_fn(**state_indices)
                                                        if interval == 1:
                                                            pathway.add_transition(
                                                                pair,
                                                                side,
                                                                interval=1,
                                                                width=width1,
                                                                deph=deph1,
                                                            )
                                                        else:
                                                            pathway.add_transition(
                                                                pair, side
                                                            )

                                                    transfer_target, transfer_source = (
                                                        desc.transfer_indices(
                                                            relax_a=relax_a,
                                                            relax_b=relax_b,
                                                            exc_ket=exc_ket,
                                                            exc_bra=exc_bra,
                                                        )
                                                    )
                                                    pathway.add_transfer(
                                                        transfer_target,
                                                        transfer_source,
                                                    )
                                                    pathway.set_evolution_factor(
                                                        evolution_factor
                                                    )

                                                    for (
                                                        pair_fn,
                                                        side,
                                                        interval,
                                                    ) in (
                                                        desc.transitions_after_transfer
                                                    ):
                                                        pair = pair_fn(**state_indices)
                                                        if interval == 3:
                                                            pathway.add_transition(
                                                                pair,
                                                                side,
                                                                interval=3,
                                                                width=width3,
                                                                deph=deph3,
                                                            )
                                                        else:
                                                            pathway.add_transition(
                                                                pair, side
                                                            )

                                                except Exception:
                                                    raise Exception(
                                                        f"Generation of {desc.pathway_name}"
                                                        " pathway failed"
                                                    )

                                                pathway.build()
                                                pathways.append(pathway)


def _generate_ground_excited_pathway(
    self: AggregateSpectroscopy,
    desc: _GroundExcitedPathwayDesc,
    pathways: list,
    evolution_superop: numpy.ndarray,
    population_tol: float,
    dipole_tol: float,
    verbose: int = 0,
) -> None:
    ground_states = self.get_electronic_groundstate()
    excited_states = self.get_excitonic_band(band=1)

    if verbose > 0:
        print(f"Liouville pathway {desc.pathway_name}")
        print("Population tolerance: ", population_tol)
        print("Dipole tolerance:     ", dipole_tol)

    for ground in ground_states:
        if verbose > 0:
            print("Ground state: ", ground, "of", len(ground_states))

        if self.rho0[ground, ground] > population_tol:
            for exc_ket in excited_states:
                if verbose > 1:
                    print("Excited state: ", exc_ket, "of", len(excited_states))

                if self.D2[exc_ket, ground] > dipole_tol:
                    for mid_ground in ground_states:
                        if self.D2[mid_ground, exc_ket] > dipole_tol:
                            evolution_indices = desc.evolution_factor_indices(
                                ground, mid_ground
                            )
                            evolution_factor = evolution_superop[evolution_indices]

                            for final_exc in excited_states:
                                if desc.detection_dipole_check(
                                    self, final_exc, ground, mid_ground, dipole_tol
                                ):
                                    try:
                                        if verbose > 5:
                                            print(
                                                f" * Generating {desc.pathway_name}",
                                                ground,
                                                exc_ket,
                                            )

                                        pathway = diag.liouville_pathway(
                                            desc.pathway_type,
                                            ground,
                                            aggregate=self,
                                            order=3,
                                            pname=desc.pathway_name,
                                        )

                                        first_lineshape_pair = (
                                            desc.first_transition_pair(
                                                exc_ket=exc_ket,
                                                ground=ground,
                                                mid_ground=mid_ground,
                                                final_exc=final_exc,
                                            )
                                        )
                                        width1 = self.get_transition_width(
                                            first_lineshape_pair
                                        )
                                        deph1 = self.get_transition_dephasing(
                                            first_lineshape_pair
                                        )

                                        third_lineshape_pair = (
                                            desc.third_transition_pair(
                                                exc_ket=exc_ket,
                                                ground=ground,
                                                mid_ground=mid_ground,
                                                final_exc=final_exc,
                                            )
                                        )
                                        width3 = self.get_transition_width(
                                            third_lineshape_pair
                                        )
                                        deph3 = self.get_transition_dephasing(
                                            third_lineshape_pair
                                        )

                                        state_indices = dict(
                                            ground=ground,
                                            exc_ket=exc_ket,
                                            mid_ground=mid_ground,
                                            final_exc=final_exc,
                                        )
                                        for i, (pair_fn, side, interval) in enumerate(
                                            desc.transitions
                                        ):
                                            pair = pair_fn(**state_indices)
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
                                            pathway.add_transition(pair, side, **kwargs)

                                        pathway.set_evolution_factor(evolution_factor)

                                    except Exception:
                                        raise Exception(
                                            f"Generation of {desc.pathway_name} pathway failed"
                                        )

                                    pathway.build()
                                    pathways.append(pathway)


def generate_R1g(
    self: AggregateSpectroscopy,
    pathways: list,
    evolution_superop: numpy.ndarray,
    population_tol: float,
    dipole_tol: float,
    evolution_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self,
        _R1G_DESC,
        pathways,
        evolution_superop,
        population_tol,
        dipole_tol,
        evolution_tol,
        verbose,
    )


def generate_R1gE(
    self: AggregateSpectroscopy,
    pathways: list,
    evolution_superop: numpy.ndarray,
    population_tol: float,
    dipole_tol: float,
    evolution_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self,
        _R1GE_DESC,
        pathways,
        evolution_superop,
        population_tol,
        dipole_tol,
        evolution_tol,
        verbose,
    )


def generate_R2g(
    self: AggregateSpectroscopy,
    pathways: list,
    evolution_superop: numpy.ndarray,
    population_tol: float,
    dipole_tol: float,
    evolution_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self,
        _R2G_DESC,
        pathways,
        evolution_superop,
        population_tol,
        dipole_tol,
        evolution_tol,
        verbose,
    )


def generate_R2gE(
    self: AggregateSpectroscopy,
    pathways: list,
    evolution_superop: numpy.ndarray,
    population_tol: float,
    dipole_tol: float,
    evolution_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self,
        _R2GE_DESC,
        pathways,
        evolution_superop,
        population_tol,
        dipole_tol,
        evolution_tol,
        verbose,
    )


def generate_R3g(
    self: AggregateSpectroscopy,
    pathways: list,
    evolution_superop: numpy.ndarray,
    population_tol: float,
    dipole_tol: float,
    verbose: int = 0,
) -> None:
    _generate_ground_excited_pathway(
        self,
        _R3G_DESC,
        pathways,
        evolution_superop,
        population_tol,
        dipole_tol,
        verbose,
    )


def generate_R4g(
    self: AggregateSpectroscopy,
    pathways: list,
    evolution_superop: numpy.ndarray,
    population_tol: float,
    dipole_tol: float,
    verbose: int = 0,
) -> None:
    _generate_ground_excited_pathway(
        self,
        _R4G_DESC,
        pathways,
        evolution_superop,
        population_tol,
        dipole_tol,
        verbose,
    )


def generate_R1f(
    self: AggregateSpectroscopy,
    pathways: list,
    evolution_superop: numpy.ndarray,
    population_tol: float,
    dipole_tol: float,
    evolution_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self,
        _R1F_DESC,
        pathways,
        evolution_superop,
        population_tol,
        dipole_tol,
        evolution_tol,
        verbose,
    )


def generate_R2f(
    self: AggregateSpectroscopy,
    pathways: list,
    evolution_superop: numpy.ndarray,
    population_tol: float,
    dipole_tol: float,
    evolution_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self,
        _R2F_DESC,
        pathways,
        evolution_superop,
        population_tol,
        dipole_tol,
        evolution_tol,
        verbose,
    )


def generate_R1fE(
    self: AggregateSpectroscopy,
    pathways: list,
    evolution_superop: numpy.ndarray,
    population_tol: float,
    dipole_tol: float,
    evolution_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self,
        _R1FE_DESC,
        pathways,
        evolution_superop,
        population_tol,
        dipole_tol,
        evolution_tol,
        verbose,
    )


def generate_R2fE(
    self: AggregateSpectroscopy,
    pathways: list,
    evolution_superop: numpy.ndarray,
    population_tol: float,
    dipole_tol: float,
    evolution_tol: float,
    verbose: int = 0,
) -> None:
    _generate_relaxation_pathway(
        self,
        _R2FE_DESC,
        pathways,
        evolution_superop,
        population_tol,
        dipole_tol,
        evolution_tol,
        verbose,
    )


def generate_1orderP_sec(
    self: Any,
    pathways: list,
    population_tol: float,
    dipole_tol: float,
    verbose: int = 0,
) -> None:

    ground_states = self.get_electronic_groundstate()
    excited_states = self.get_excitonic_band(band=1)

    if verbose > 0:
        print("Liouville pathway of first order")
        print("Population tolerance: ", population_tol)
        print("Dipole tolerance:     ", dipole_tol)

    k = 0
    l = 0
    for ground in ground_states:
        if verbose > 0:
            print("Ground state: ", ground, "of", len(ground_states))

        # Only thermally allowed starting states are considered
        if self.rho0[ground, ground] > population_tol:
            for exc_ket in excited_states:
                if self.D2[exc_ket, ground] > dipole_tol:
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
                            print(" * Generating P1", ground, exc_ket)

                        pathway = diag.liouville_pathway(
                            "NR",
                            ground,
                            aggregate=self,
                            order=1,
                            pname="P1",
                            popt_band=1,
                            relax_order=1,
                        )

                        # first transition lineshape
                        width1 = self.get_transition_width((exc_ket, ground))
                        deph1 = self.get_transition_dephasing((exc_ket, ground))

                        #      |g_i1> <g_i1|
                        pathway.add_transition(
                            (exc_ket, ground), +1, interval=1, width=width1, deph=deph1
                        )
                        #      |e_i2> <g_i1|
                        pathway.add_transition(
                            (ground, exc_ket), +1, interval=1, width=width1, deph=deph1
                        )
                        #      |g_i1> <g_i1|

                    except Exception:
                        break

                    pathway.build()
                    pathways.append(pathway)
                    k += 1
