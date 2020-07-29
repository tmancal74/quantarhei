# -*- coding: utf-8 -*-
"""
Class representing aggregates of molecules.

The class enables building of complicated objects from objects of the Molecule
type, their mutual interactions and system-bath interactions. It also provides
an interface to various methods of open quantum systems theory.


"""

import numpy
#import h5py

from ..core.managers import UnitsManaged
#from ..core.units import cm2int
from .interactions import dipole_dipole_interaction

from ..qm.oscillators.ho import fcstorage
from ..qm.oscillators.ho import operator_factory

from ..qm.hilbertspace.operators import Operator
from ..qm.hilbertspace.operators import DensityMatrix
from ..qm.hilbertspace.operators import ReducedDensityMatrix
from ..qm.hilbertspace.statevector import StateVector
from ..qm.propagators.dmevolution import DensityMatrixEvolution
from ..qm.propagators.dmevolution import ReducedDensityMatrixEvolution
from ..qm.propagators.statevectorevolution import StateVectorEvolution
from ..qm.liouvillespace.systembathinteraction import SystemBathInteraction
from ..qm.hilbertspace.hamiltonian import Hamiltonian
from ..qm.hilbertspace.dmoment import TransitionDipoleMoment

from ..qm.corfunctions import CorrelationFunctionMatrix

#from .aggregate_states import aggregate_state
from .aggregate_states import ElectronicState
from .aggregate_states import VibronicState

#from ..core.managers import energy_units
#from .molecules import Molecule
from ..core.managers import Manager
from ..core.saveable import Saveable

import quantarhei as qr


class AggregateBase(UnitsManaged, Saveable):
    """ Molecular aggregate


    Parameters
    ----------

    name : str
        Specifies the name of the aggregate

    molecules : list or tuple
        List of molecules out of which the aggregate is built

    """

    def __init__(self, molecules=None, name=""):

        self.mnames = {}    #
        self.monomers = []
        self.nmono = 0         #
        self.name = name       #
        self.mult = 0          #

        self._has_egcf_matrix = False #
        self.egcf_matrix = None

        self._has_system_bath_interaction = False #

        self._has_lindich_axes = False

        self.coupling_initiated = False #
        self.resonance_coupling = None

        if molecules is not None:
            for m in molecules:
                self.add_Molecule(m)

        self._init_me()


    def _init_me(self):
        """Initializes all built attributes of the aggregate

        This should put the object into a pre-build state


        """

        self.FC = fcstorage()
        self.ops = operator_factory()

        self._has_relaxation_tensor = False #

        self._relaxation_theory = "" #

        self._built = False     #
        self._diagonalized = False

        self.mult = 0                #
        self.sbi_mult = 0            #
        self.Nel = 0                 #
        self.Ntot = 0                #
        self.Nb = 0                  #

        self.vibindices = []
        self.which_band = None
        self.elsigs = None

        self.HH = None
        self.HamOp = None
        self.DD = None
        self.Wd = None
        self.Dr = None
        self.D2 = None
        self.D2_max = 0
        self.sbi = None


    def clean(self):
        """Cleans the aggregate object of anything built

        This operation leaves the molecules of the aggregate intact and keeps
        few more pieces of information it it. E. g. coupling matrix is not
        deleted. You call build again after this.


        """
        self._init_me()


    def wipe_out(self):
        """Removes everything except of name attribute

        You have to set molecules and recalculate interactions before you can
        build


        """
        self.mnames = {}   #
        self.monomers = []
        self.nmono = 0      #
        self.mult = 0       #
        self.sbi_mult = 0    #

        self._has_egcf_matrix = False   #
        self.egcf_matrix = None

        self._has_system_bath_interaction = False  #

        self.coupling_initiated = False     #
        self.resonance_coupling = None    #

        self._init_me()


    ########################################################################
    #
    #    BUILDING METHODS
    #
    ########################################################################

    def init_coupling_matrix(self):
        """Nullifies coupling matrix


        """
        self.resonance_coupling = numpy.zeros((self.nmono,self.nmono),
                                              dtype=numpy.float64)
        self.coupling_initiated = True
        #
        # TESTED


    def set_resonance_coupling(self, i, j, coupling):
        """Sets resonance coupling value between two sites

        """
        if not self.coupling_initiated:
            self.init_coupling_matrix()

        coup = self.convert_energy_2_internal_u(coupling)

        self.resonance_coupling[i,j] = coup
        self.resonance_coupling[j,i] = coup
        #
        # TESTED


    def get_resonance_coupling(self, i, j):
        """Returns resonance coupling value between two sites

        """
        coupling = self.resonance_coupling[i,j]
        return self.convert_energy_2_current_u(coupling)
        #
        # TESTED


    def set_resonance_coupling_matrix(self, coupmat):
        """Sets resonance coupling values from a matrix

        """

        if type(coupmat) in (list, tuple):
            coupmat = numpy.array(coupmat)

        coup = self.convert_energy_2_internal_u(coupmat)
        self.resonance_coupling = coup
        if not self.coupling_initiated:
            self.coupling_initiated = True
        #
        # TESTED


    def dipole_dipole_coupling(self, kk, ll, epsr=1.0):
        """Calculates dipole-dipole coupling

        """
        if kk == ll:
            raise Exception("Only coupling between different molecules \
            can be calculated")

        #FIXME: this works only for first excited states of two-level molecules
        d1 = self.monomers[kk].dmoments[0,1,:]
        r1 = self.monomers[kk].position
        d2 = self.monomers[ll].dmoments[0,1,:]
        r2 = self.monomers[ll].position

        val =  dipole_dipole_interaction(r1, r2, d1, d2, epsr)
        return self.convert_energy_2_current_u(val)
        #
        # TESTED


    def set_coupling_by_dipole_dipole(self, epsr=1.0):
        """Sets resonance coupling by dipole-dipole interaction

        """

        if not self.coupling_initiated:
            self.init_coupling_matrix()
        for kk in range(self.nmono):
            for ll in range(kk+1,self.nmono):
                cc = self.dipole_dipole_coupling(kk,ll,epsr=epsr)
                c1 = self.convert_energy_2_internal_u(cc)
                self.resonance_coupling[kk,ll] = c1
                self.resonance_coupling[ll,kk] = c1
        #
        # TESTED


    def calculate_resonance_coupling(self, method="dipole-dipole",
                               params=dict(epsr=1.0)):
        """ Sets resonance coupling calculated by a given method

        Parameters
        ----------

        method: string
            Method to be used for calculation of resonance coupling

        """

        if method == "dipole-dipole":
            epsr = params["epsr"]
            self.set_coupling_by_dipole_dipole(epsr=epsr)
        else:
            raise Exception("Unknown method for calculation"+
                            " of resonance coupling")
        #
        # TESTED


    # FIXME: This must be set in coordination with objects describing laboratory
    def set_lindich_axes(self, axis_orthog_membrane):
        """ Creates a coordinate system with one axis supplied by the user
        (typically an axis orthogonal to the membrane), and two other axes, all
        of which are orthonormal.
        """

        qr = numpy.vstack((axis_orthog_membrane, numpy.array([1,0,0]), numpy.array([0,1,0]))).T
        self.q, r = numpy.linalg.qr(qr)
        self._has_lindich_axes = True

    # FIXME: This must be set in coordination with objects describing laboratory
    def get_lindich_axes(self):
        if self._has_lindich_axes:
            return self.q
        else:
            raise Exception("No linear dichroism coordinate system supplied")


    # FIXME: This should be delegated SystemBathInteraction
    def set_egcf_matrix(self,cm):
        """Sets a matrix describing system bath interaction

        """
        self.egcf_matrix = cm
        self._has_egcf_matrix = True


    #
    # Molecues
    #
    def add_Molecule(self, mono):
        """Adds monomer to the aggregate


        """

        # If at least one monomer has energy gap correlation function
        # we will try to build system bath interaction for a the aggregate.
        # Exception will be thrown if not all monomers have the same egcf
        if mono._has_egcf:
            self._has_system_bath_interaction = True

        self.monomers.append(mono)
        self.mnames[mono.name] = len(self.monomers)-1
        self.nmono += 1
        #
        # TESTED


    def get_Molecule_by_name(self, name):
        try:
            im = self.mnames[name]
            return self.monomers[im]
        except:
            raise Exception()


    def get_Molecule_index(self, name):
        try:
            im = self.mnames[name]
            return im
        except:
            raise Exception()


    def remove_Molecule(self, mono):
        self.monomers.remove(mono)
        self.nmono -= 1


    def get_nearest_Molecule(self,molecule):
        """Returns a molecule nearest in the aggregate to a given molecule

        Parameters
        ----------

        molecule : Molecule
            Molecule whose neighbor we look for

        """
        tol = 1.0e-3
        rmin = 1.0e20
        r1 = molecule.position
        mmin = None
        for m in self.monomers:
            r2 = m.position
            r = r1 - r2
            dist = numpy.sqrt(numpy.dot(r,r))
            if (dist > tol) and (dist < rmin):
                mmin = m
                rmin = dist

        return mmin, rmin

    #
    # Vibrational modes
    #
    def add_Mode_by_name(self,name,mode):
        try:
            im = self.mnames[name]
            mn = self.monomers[im]
            mn.add_mode(mode)
        except:
            raise Exception()

    def get_Mode_by_name(self,name,N):
        try:
            im = self.mnames[name]
            mn = self.monomers[im]
            return mn.get_mode(N)
        except:
            raise Exception("Mode not found")

    #
    # Transition dipole
    #
    def get_dipole_by_name(self,name,N,M):
        try:
            im = self.mnames[name]
            mn = self.monomers[im]
            return mn.get_dipole(N,M)
        except:
            raise Exception()

    def get_dipole(self, n, N, M):
        nm = self.monomers[n]
        return nm.get_dipole(N,M)


    #
    # Various info
    #
    def get_width(self, n, N, M):
        nm = self.monomers[n]
        return nm.get_transition_width((N,M))


    def get_max_excitations(self):
        """Returns a list of maximum number of excitations on each monomer

        """
        omax = []
        for nm in self.monomers:
            omax.append(nm.nel-1)
        return omax


    def get_energy_by_name(self, name, N):
        """ Electronic energy """
        try:
            im = self.mnames[name]
            mn = self.monomers[im]
            return mn.get_energy(N)
        except:
            raise Exception()


    def fc_factor(self, state1, state2):
        """Franck-Condon factors between two vibrational states


        Calculates Franck-Condon factor between two aggregate_states
        regardless of their electronic parts


        """

        inx1 = state1.vsig
        inx2 = state2.vsig
        sta1 = state1.elstate.vibmodes
        sta2 = state2.elstate.vibmodes

        if not (len(sta1)==len(sta2)):
            raise Exception("Incompatible states")

        res = 1.0
        for kk in range(len(sta1)):
            smod1 = sta1[kk]
            smod2 = sta2[kk]

            # difference in shifts
            shft = smod1.shift - smod2.shift
            # quantum numbers
            qn1 = inx1[kk]
            qn2 = inx2[kk]

            #calculate FC factors
            #
            #Best implementation would be a table look-up. First we calculate
            #a table of FC factors from known omegas and shifts and here we
            #just consult the table.

            if not self.FC.lookup(shft):
                fc = self.ops.shift_operator(shft)[:20,:20]
                self.FC.add(shft,fc)

            ii = self.FC.index(shft)
            rs = self.FC.get(ii)[qn1,qn2]

            res = res*rs

        return res


    def get_transition_width(self, state1, state2=None):
        """Returns phenomenological width of a given transition


        Parameters
        ----------

        state1 : {ElectroniState/VibronicState, tuple}
            If both state1 and state2 are specified, it is assumed they are
            of the type of Electronic of Vibronic state. Otherwise, if state2
            is None, it is assumed that it is a tuple representing
            a transition

        state2 : {ElectroniState/VibronicState, None}
        If not None it is of the type of Electronic of Vibronic state

        """
        if state2 is not None:

            b1 = state1.elstate.band
            b2 = state2.elstate.band

            if abs(b2-b1) == 1:

                # index of a monomer on which the transition occurs
                exindx = self._get_exindx(state1, state2)
                width = self.monomers[exindx].get_transition_width((0,1))
                #print(exindx, width)
                return width

            elif abs(b2-b1) == 2:
                sig1 = state1.elstate.get_signature()
                sig2 = state2.elstate.get_signature()
                if (2 in sig1) or (2 in sig2):
                    exindx = self._get_exindx(state1, state2)
                    #print("Two-ex:", exindx)
                    # FIXME: The factor of 2 needs to be checked and justified
                    width = \
                    2.0*self.monomers[exindx].get_transition_width((0,1))
                    return width
                else:
                    (indx1, indx2) = self._get_twoexindx(state1, state2)
                    #print(state1.elstate.elsignature,
                    #      state2.elstate.elsignature, indx1, indx2)
                    width = self.monomers[indx1].get_transition_width((0,1))
                    width += self.monomers[indx2].get_transition_width((0,1))
                    #print(indx1, indx2, width)
                    return width

            else:
                return -1.0


        else:

            transition = state1

            Nf = transition[0]
            Ni = transition[1]

            eli = self.elinds[Ni]
            elf = self.elinds[Nf]

            # g -> 1 exciton band transitions
            if (self.which_band[eli] == 0) and (self.which_band[elf] == 1):
                # this simulates bath correlation function
                #print("0->1 :", self.Wd[Nf, Nf]**2)
                return self.Wd[Nf, Nf]**2

            elif (self.which_band[eli] == 1) and (self.which_band[elf] == 0):
                # this simulates bath correlation function
                #print("0->1 :", self.Wd[Nf, Nf]**2)
                return self.Wd[Ni, Ni]**2

            # 1 exciton -> 2 exciton transitions
            elif (self.which_band[eli] == 1) and (self.which_band[elf] == 2):
                # this simulates the term  g_ff + g_ee - 2Re g_fe
                ret =  (self.Wd[Ni, Ni]**2 + self.Wd[Nf, Nf]**2
                        - 2.0*(self.Wd[Nf, Ni]**2))
                #print("1->2 (", eli, elf,") :", ret, self.Wd[Nf, Ni]**2)
                return ret

            elif (self.which_band[eli] == 2) and (self.which_band[elf] == 1):
                # this simulates the term  g_ff + g_ee - 2Re g_fe
                ret =  (self.Wd[Ni, Ni]**2 + self.Wd[Nf, Nf]**2
                        - 2.0*(self.Wd[Nf, Ni]**2))
                #print("1->2 (", eli, elf,") :", ret, self.Wd[Nf, Ni]**2)
                return ret


            else:
                print("This should not be used")
                return 0.0


    def get_transition_dephasing(self, state1, state2=None):
        """Returns phenomenological dephasing of a given transition



        Parameters
        ----------

        state1 : {ElectroniState/VibronicState, tuple}
            If both state1 and state2 are specified, it is assumed they are
            of the type of Electronic of Vibronic state. Otherwise, if state2
            is None, it is assumed that it is a tuple representing
            a transition

        state2 : {ElectroniState/VibronicState, None}
        If not None it is of the type of Electronic of Vibronic state


        """
        if state2 is not None:

            # index of a monomer on which the transition occurs
            exindx = self._get_exindx(state1, state2)
            if exindx < 0:
                return 0.0

            deph = self.monomers[exindx].get_transition_dephasing((0,1))
            return deph

        else:

            transition = state1

            Nf = transition[0]
            Ni = transition[1]

            eli = self.elinds[Ni]
            elf = self.elinds[Nf]

            # g -> 1 exciton band transitions
            if (self.which_band[eli] == 0) and (self.which_band[elf] == 1):
                return self.Dr[Nf, Nf]**2

            elif (self.which_band[eli] == 1) and (self.which_band[elf] == 0):
                return self.Dr[Ni, Ni]**2

            # 1 exciton -> 2 exciton band transitions
            elif (self.which_band[eli] == 1) and (self.which_band[elf] == 2):
                # this simulates the term  g_ff + g_ee - 2Re g_fe
                return (self.Dr[Ni, Ni]**2 + self.Dr[Nf, Nf]
                        - 2.0*self.Dr[Nf, Ni])

            elif (self.which_band[eli] == 2) and (self.which_band[elf] == 1):
                # this simulates the term  g_ff + g_ee - 2Re g_fe
                return (self.Dr[Ni, Ni]**2 + self.Dr[Nf, Nf]
                        - 2.0*self.Dr[Nf, Ni])

            else:
                return -1.0


    def transition_dipole(self, state1, state2):
        """ Transition dipole moment between two states

        Parameters
        ----------
        state1 : class VibronicState
            state 1

        state2 : class VibronicState
            state 2

        """
        exindx = self._get_exindx(state1, state2)

        if (exindx < 0):
            return 0.0

        eldip = self.get_dipole(exindx, 0, 1)

        # Franck-Condon factor between the two states
        fcfac = self.fc_factor(state1,state2)

        return eldip*fcfac

    def _get_twoexindx(self, state1, state2):
        """ Indices of two molecule with transitions or negative number
        if not found

        Parameters
        ----------
        state1 : class VibronicState
            state 1

        state2 : class VibronicState
            state 2

        """
        # get electronic signatures
        els1 = state1.elstate.elsignature
        els2 = state2.elstate.elsignature

        # only states in neighboring bands can be connected by dipole moment
        b1 = state1.elstate.band
        b2 = state2.elstate.band
        if (abs(b1-b2) != 2):
            return -1

        # count the number of differences
        l = 0
        count = 0
        for kk in els1:
            if kk != els2[l]:
                count += 1
            l += 1


        if count != 2:
            return -1

        # now that we know that the states differ by two excitations, let
        # us find on which molecule they are
        exstates = []
        exindxs = []
        l = -1
        for kk in els1: # signature is just a tuple; iterate over it
            l += 1
            if kk != els2[l]: # this is the index where they differ
                # which of them is excited
                if kk > els2[l]:
                    exstates.append(els1)
                else:
                    exstates.append(els2)
                exindxs.append(l)

        if len(exstates) == 0:
            raise Exception()

        return exindxs[0], exindxs[1]


    def _get_exindx(self, state1, state2):
        """ Index of molecule with transition or negative number if not found

        Parameters
        ----------
        state1 : class VibronicState
            state 1

        state2 : class VibronicState
            state 2

        """

        # get electronic signatures
        els1 = state1.elstate.elsignature
        els2 = state2.elstate.elsignature

        # only states in neighboring bands can be connected by dipole moment
        b1 = state1.elstate.band
        b2 = state2.elstate.band
        if (abs(b1-b2) != 1) and (abs(b1-b2) != 2):
            return -1

        # count the number of differences
        l = 0
        count = 0
        for kk in els1:
            if kk != els2[l]:
                count += 1
            l += 1


        if count != 1:
            return -1

        # now that we know that the states differ by one excitation, let
        # us find on which molecule it is
        exstate = None
        l = -1
        for kk in els1: # signature is just a tuple; iterate over it
            l += 1
            if kk != els2[l]: # this is the index where they differ
                # which of them is excited
                if kk > els2[l]:
                    exstate = els1
                else:
                    exstate = els2
                exindx = l

        if exstate is None:
            raise Exception()

        return exindx


    def total_number_of_states(self, mult=1, vibgen_approx=None, Nvib=None,
                               vibenergy_cutoff=None, save_indices=False):
        """ Total number of states in the aggregate

        Counts all states of the aggregate by iterating through them. States
        are generated with a set of constraints.

        """

        nret = 0

        for state in self.allstates(mult=mult,
                                    save_indices=save_indices,
                                    vibgen_approx=vibgen_approx,
                                    Nvib=Nvib,
                                    vibenergy_cutoff=vibenergy_cutoff):
            nret += 1

        return nret


    def total_number_of_electronic_states(self, mult=1):
        """ Total number of electronic states in the aggregate"""

        nret = 0

        for elsig in self.elsignatures(mult=mult):
            nret += 1

        return nret


    def number_of_states_in_band(self, band=1, vibgen_approx=None,
                                 Nvib=None, vibenergy_cutoff=None):
        """ Number of states in a given excitonic band """

        nret = 0

        for state in self.allstates(mult=band, mode="EQ", save_indices=False,
                                    vibgen_approx=vibgen_approx, Nvib=Nvib,
                                    vibenergy_cutoff=vibenergy_cutoff):
            nret += 1

        return nret


    def number_of_electronic_states_in_band(self, band=1):
        """ Number of states in a given excitonic band """

        nret = 0

        for elsig in self.elsignatures(mult=band, mode="EQ"):
            nv = 1
            nret += nv

        return nret


    def get_ElectronicState(self, sig, index=None):
        """Returns electronic state corresponding to this aggregate

        Parameters
        ----------

        sig : tuple
            Tuple defining the electronic state of the aggregate

        index : integer or None
            If integer is specified, this number is recorded as an index
            of this state in the aggregate. It is recorded only internally
            in the state object. Aggregate keeps its own record which is
            created during the build.

        """
        return ElectronicState(self, sig, index)


    def get_VibronicState(self, esig, vsig):
        """Returns vibronic state corresponding to the two specified signatures

        """
        elstate = self.get_ElectronicState(sig=esig)
        return VibronicState(elstate, vsig)


    def coupling(self, state1, state2, full=False):
        """Coupling between two aggregate states


        Parameters
        ----------

        state1 : {ElectronicState, VibronicState}
            States for which coupling should be calculated


        """

        #
        # Coupling between two purely electronic states
        #
        if (isinstance(state1, ElectronicState)
           and isinstance(state2, ElectronicState)):

            if self.nmono > 1:
                # coupling within the bands
                if state1.band == state2.band:
                    #print("Band:", state1.band)

                    if state1.band == 1:

                        kk = state1.index - 1
                        ll = state2.index - 1

                        if (kk >= 0) and (ll >= 0):
                            coup = self.resonance_coupling[kk,ll]
                        else:
                            coup = 0.0

                    else:

                        els1 = state1.elsignature
                        els2 = state2.elsignature
                        Ns = len(els1)
                        sites = [0,0]
                        k = 0
                        # count differences
                        for i in range(Ns):
                            if els1[i] != els2[i]:
                                if (k == 0) or (k == 1):
                                    sites[k] = i
                                k += 1
                        # if there are exactly 2 differences, the differing
                        # two molecules are those coupled; sites[k] contains
                        # indiced those coupled molecules
                        if k == 2:
                            kk = sites[0]
                            ll = sites[1]
                            coup = self.resonance_coupling[kk,ll]
                        else:
                            coup = 0.0

                elif state1.band + 2 == state2.band:

                    #print(state1.elsignature, state2.elsignature)
                    coup = 0.0

                else:
                    coup = 0.0

            else:
                coup = 0.0

        #
        # Coupling between two general states
        #
        elif (isinstance(state1, VibronicState)
          and isinstance(state2, VibronicState)):

            es1 = state1.elstate
            es2 = state2.elstate

            fc = self.fc_factor(state1, state2)

            # it make sense to calculate coupling only when the number
            # of molecules is larger than 1
            if self.nmono > 1:

                # coupling within the bands
                if es1.band == es2.band:

                    # single exciton band
                    if es1.band == 1:

                        kk = es1.index - 1
                        ll = es2.index - 1

                        if (kk >= 0) and (ll >= 0):
                            coup = self.resonance_coupling[kk,ll]*fc
                        else:
                            coup = 0.0

                    else:

                        els1 = es1.elsignature
                        els2 = es2.elsignature
                        Ns = len(els1)
                        sites = [0,0]
                        k = 0
                        # count differences
                        for i in range(Ns):
                            if els1[i] != els2[i]:
                                if (k == 0) or (k == 1):
                                    sites[k] = i
                                k += 1
                        # if there are exactly 2 differences, the differing
                        # two molecules are those coupled; sites[k] contains
                        # indiced those coupled molecules
                        if k == 2:
                            kk = sites[0]
                            ll = sites[1]
                            #print(kk,ll,els1,els2)
                            ar1 = numpy.array(els1)
                            ar2 = numpy.array(els2)
                            df = numpy.abs(ar1-ar2)
                            sdf = numpy.sum(df)
                            if (sdf == 2):
                                mx1 = numpy.max([ar1[kk],ar2[kk]])
                                mx2 = numpy.max([ar1[ll],ar2[ll]])
                                #print("max:",mx1,mx2)
                                harm_fc = numpy.sqrt(numpy.real(mx1))
                                harm_fc = harm_fc*numpy.sqrt(numpy.real(mx2))
                                fc = fc*harm_fc
                                #print(harm_fc)
                                coup = self.resonance_coupling[kk,ll]*fc
                            else:
                                coup = 0.0
                        else:
                            coup = 0.0

                elif (numpy.abs(es1.band - es2.band) == 2) and full:

                    #print(es1.elsignature, es2.elsignature)
                    #print("Here we calculate coupling between bands")
                    els1 = es1.elsignature
                    els2 = es2.elsignature
                    Ns = len(els1)
                    sites = [0,0]
                    k = 0
                    # count differences
                    for i in range(Ns):
                        if els1[i] != els2[i]:
                            if (k == 0) or (k == 1):
                                sites[k] = i
                            k += 1
                    # if there are exactly 2 differences, the differing
                    # two molecules are those coupled; sites[k] contains
                    # indiced those coupled molecules
                    if k == 2:
                        kk = sites[0]
                        ll = sites[1]
                        coup = self.resonance_coupling[kk,ll]*fc
                    else:
                        coup = 0.0

                else:
                    coup = 0.0
            else:
                coup = 0.0

        return self.convert_energy_2_current_u(coup)



    #######################################################################
    #
    # Generators
    #
    #######################################################################

    def elsignatures(self, mult=1, mode="LQ", emax=None):
        """ Generator of electronic signatures

        Here we create signature tuples of electronic states. The signature
        is a tuple with as many integer numbers as the members of
        the aggregate. Each integer represents the state in which the
        member of the aggregate is, e.g. 0 for ground state, 1 for the first
        excited state etc.


        Parameters
        ----------
        mult : int
            multiplicity of excitons

        mode : str {"LQ", "EQ"}
            mode of the functions.

            mode="LQ" returns all signatures of states with
            multiplicity less than or equal to the `mult`

            mode="EQ" returns signatures of states with a multiplicity
            given by `mult`

        """

        if mode not in ["LQ","EQ"]:
            raise Exception("Unknown mode")

        l = len(self.monomers)

        # list of maximum numbers of excitations on each sites
        if emax is None:
            omax = self.get_max_excitations()
        else:
            omax = emax

        if mult < 0:
            raise Exception("mult must be larger than or equal to zero")

        mlt = 0
        # iterate over all excition multiplicities
        while mlt <= mult:
            # no excitations (ground state)
            out = [0 for k in range(l)]
            # if this is the multiplicity 0, yield the ground state
            if (((mlt == 0) and (mode == "LQ")) or (mult==0)):
                yield tuple(out)
            else:
                k = 1
                # first we have only ground state signature
                ins = [out]
                strt = [0]
                while k <= mlt:
                    nins = []
                    nstr = []
                    # take all signatures in "ins" and add one excitation
                    for out_added, last in self._add_excitation(ins,strt,omax):
                        # if mlt excitation was added yield
                        if (((k == mlt) and (mode == "LQ"))
                          or((mult == k) and (mult == mlt))):
                            yield tuple(out_added)
                        else:
                            # make a list of all new signatures
                            nins.append(out_added)
                            # for each signature, save the index
                            # on which an excitation was added last
                            nstr.append(last)
                    # set the new signatures for processing in the iteration
                    ins = nins
                    strt = nstr
                    k += 1
            mlt += 1


    def _add_excitation(self, inlists, strt, omax):
        """Adds one excitation to all submitted electronic signatures"""

        k = 0
        # go through all signatures
        for inlist in inlists:
            l = len(inlist)
            if len(omax) != l:
                raise Exception("arg omax has to be a list of the same \
                length as arg inlist")
            # go through all positions from the last index on (in order
            # to create unique signatures)
            for i in range(strt[k],l):
                # if it is possible to add an excitation
                # make a new list and add
                if inlist[i] < omax[i]:
                    out = inlist.copy()
                    out[i] += 1
                    # yield the list and the index of the last added exitation
                    yield out, i
            k += 1


    def vibsignatures(self, elsignature, approx=None):
        """ Generator of vibrational signatures

        Parameters
        ----------

        approx : None or str
            Approximation used in generation of vibrational states
            Allowed values are None or "SPA"

        """
        cs = ElectronicState(self, elsignature)
        return cs.vsignatures(approx=approx)


    def allstates(self, mult=1, mode="LQ", all_vibronic=True,
                  save_indices=False, vibgen_approx=None, Nvib=None,
                  vibenergy_cutoff=None):
        """ Generator of all aggregate states

        Iterator generating all states of the aggregate given a set
        of constraints.


        Parameters
        ----------

        mult : integer {0, 1, 2}
            Exciton multiplicity (ground state, single and double excitons).
            All excitons with the multiplicity smaller or equal to ``mult``
            are generated by default

        mode : str {"LQ", "EQ"}
            If set to "LQ" generates excitons with smaller or equal
            multiplicity than specified. If "EQ" is specified, generates only
            excitons with given multiplicity

        save_indices : bool
            If True, saves indices of all generated states, so that they can
            be later used.

        all_vibronic : bool
            If True, all generated states are of the type ``VibronicState``,
            even if no vibrational modes are specified. If False,
            ``ElectronicState`` is returned for pure electronic states

        vibgen_approx : str {"ZPA", "SPA", "TPA", "NPA", "SPPMA", "TPPMA", "NPPMA"}
            Type of approximation in generating vibrational states

        Nvib : integer
            Number of vibrational states that goes into "NPA" and "NPPMA"
            approximations

        vibenergy_cutoff: float
            Maximum vibrational energy allowed in generation of vibrational
            states


        """
        ast = 0  # index counting all states
        ist = 0  # index counting electronic states

        # run over all electronic signatures
        for ess1 in self.elsignatures(mult=mult, mode=mode):

            # generate electronic state
            es1 = self.get_ElectronicState(ess1, ist)

            # loop over all vibrational signatures in electronic states
            nsig = 0
            for vsig1 in es1.vsignatures(approx=vibgen_approx, N=Nvib,
                                         vibenergy_cutoff=vibenergy_cutoff):

                # create vibronic state with a given signature
                s1 = VibronicState(es1, vsig1)

                if save_indices:
                    # save indices corresponding to vibrational sublevels
                    # of a given electronic state
                    self.vibindices[ist].append(ast)
                    self.vibsigs[ast] = (ess1, vsig1)
                    self.elinds[ast] = ist

                yield ast ,s1

                ast += 1 # count all states
                nsig += 1 # count the number of vibrational signatures

            # if no vibrational signatures
            if nsig == 0:
                # if True return vibronic states even
                # for purely electronic state
                if all_vibronic:
                    s1 = VibronicState(es1, None)
                else:
                    s1 = es1

            if save_indices:
                # save electronic signatures to be searchable later
                self.elsigs[ist] = ess1
                # save the band to which this electronic index corresponds
                self.which_band[ist] = numpy.sum(ess1)

            ist += 1 # count electronic states


    def elstates(self, mult=1, mode="LQ", save_indices=False):
        """ Generator of electronic states

        """
        a = 0
        for ess1 in self.elsignatures(mult=mult, mode=mode):
            es1 = self.get_ElectronicState(ess1, a)
            yield a,es1
            a += 1



    def __str__(self):
        out  = "\nquantarhei.Aggregate object"
        out += "\n==========================="
        out += "\nname = %s" % self.name
        out += "\nnumber of molecules = %i " % self.nmono
        count = 0
        for nm in self.monomers:
            out += "\n\nMonomer %i" % count
            out += str(nm)
            count += 1

        out += "\n\nResonance coupling matrix: "
        out +=   "\n-------------------------- "
        out += "\n"+str(self.resonance_coupling)

        out += "\n\nAggregate built = "+str(self._built)

        out +="\n\nSelected attributes"
        out +="\n--------------------"
        out +="\nmult = "+str(self.mult)
        out +="\nNel  = "+str(self.Nel)
        out +="\nNtot = "+str(self.Ntot)


        return out


    ###########################################################################
    #
    #    BUILDING
    #
    ###########################################################################

    def build(self, mult=1, sbi_for_higher_ex=False,
              vibgen_approx=None, Nvib=None, vibenergy_cutoff=None,
              fem_full=False):
        """Builds aggregate properties

        Calculates Hamiltonian and transition dipole moment matrices and
        sets up system-bath interaction

        Parameters
        ----------

        mult : int
            exciton multiplicity

        sbi_for_higher_ex: bool
            If set True, system-bath information is explicitely created for
            higher exciton states (consistent with the specified parameters
            `mult`). If set False, it is expected that if system-bath
            interaction for higher excitons is needed, it will be reconstructed
            from the single exciton part of this object

        vibge_approx:
            Approximation used in the generation of vibrational state.

        """
        manager = Manager()
        manager.set_current_units("energy", "int")

        # maximum multiplicity of excitons handled by this aggregate
        self.mult = mult
        if sbi_for_higher_ex:
            self.sbi_mult = mult
        else:
            self.sbi_mult = 1

        #######################################################################
        #
        # Electronic and vibrational states
        #
        #######################################################################

        # total number of electronic states
        self.Nel = self.total_number_of_electronic_states(mult=mult)

        # storage for indices of vibrational states
        self.vibindices = []
        # there are as many lists of indices as there are electronic states
        for i in range(self.Nel):
            self.vibindices.append([])

        # number of states in the aggregate (taking into account
        # approximations in generation of vibrational states)
        Ntot = self.total_number_of_states(mult=mult,
                                           vibgen_approx=vibgen_approx,
                                           Nvib=Nvib, save_indices=False,
                                           vibenergy_cutoff=vibenergy_cutoff)
        # save total number of states (including vibrational)
        self.Ntot = Ntot
        # information about the band to which a state belongs
        self.which_band = numpy.zeros(self.Ntot, dtype=numpy.int)
        # electronic signature for every state
        self.elsigs = [None]*self.Nel
        # vibrational signature for each state
        self.vibsigs = [None]*self.Ntot
        # FIXME: what is this???
        self.elinds = numpy.zeros(self.Ntot, dtype=numpy.int)
        # Hamiltonian matrix
        HH = numpy.zeros((Ntot, Ntot), dtype=numpy.float64)
        # Transition dipole moment matrix
        DD = numpy.zeros((Ntot, Ntot, 3),dtype=numpy.float64)
        # Matrix of Franck-Condon factors
        FC = numpy.zeros((Ntot, Ntot), dtype=numpy.float64)
        # Matrix of the transition widths (their square roots)
        Wd = numpy.zeros((Ntot, Ntot), dtype=qr.REAL)
        # Matrix of dephasing rates
        Dr = numpy.zeros((Ntot, Ntot), dtype=qr.REAL)
        # Matrix of dephasing transformation coefficients
        self.Xi = numpy.zeros((Ntot, self.Nel), dtype=qr.REAL)

        # electronic indices if twice excited state (zero for all other states)
        twoex_indx = numpy.zeros((Ntot, 2), dtype=numpy.int)

        # Initialization of the matrix of couplings between states
        if not self.coupling_initiated:
            self.init_coupling_matrix()

        Ntot = self.total_number_of_states(mult=mult,
                                           vibgen_approx=vibgen_approx,
                                           Nvib=Nvib, save_indices=True,
                                           vibenergy_cutoff=vibenergy_cutoff)

        self.all_states = []

        for a, s1 in self.allstates(mult=self.mult,
                                    vibgen_approx=vibgen_approx, Nvib=Nvib,
                                    vibenergy_cutoff=vibenergy_cutoff):
            self.all_states.append((a, s1))


        # Set up Hamiltonian and Transition dipole moment matrices
        for a, s1 in self.all_states: #self.allstates(mult=self.mult,
                     #               vibgen_approx=vibgen_approx, Nvib=Nvib,
                     #               vibenergy_cutoff=vibenergy_cutoff):

            if a == 0:
                s0 = s1

            # diagonal Hamiltonian elements
            HH[a,a] = s1.energy()

            # get dephasing and width from the ground-state
            # for each excited state
            elind = self.elinds[a]
            if (self.which_band[elind] == 1) or (self.which_band[elind] == 2):
                
                Wd[a,a] = numpy.sqrt(self.get_transition_width(s1, s0)) 
                
                Dr[a,a] = numpy.sqrt(self.get_transition_dephasing(s1, s0))

            # save composition of twice excited states
            if self.which_band[elind] == 2:
                # k_s counts excited molecules in the doubly exc. state
                # there are molecules 0 and 1 in diad (n,m)
                k_s = 0
                # counts positons in the electronic signature
                # i.e. it counts molecular index
                sig_position = 0
                for i_s in s1.elstate.elsignature:
                    if i_s == 1:
                        # we save indices of electronic states and
                        # 0 is taken by the ground state
                        twoex_indx[a, k_s] = sig_position + 1
                        k_s += 1
                    sig_position += 1


            for b, s2 in self.all_states: #self.allstates(mult=self.mult,
                                    #vibgen_approx=vibgen_approx, Nvib=Nvib,
                                    #vibenergy_cutoff=vibenergy_cutoff):

                DD[a,b,:] = numpy.real(self.transition_dipole(s1, s2))
                FC[a,b] = numpy.real(self.fc_factor(s1, s2))

                if a != b:
                    HH[a,b] = numpy.real(self.coupling(s1, s2, full=fem_full))

        #
        ###3
        #
        #print("Hamiltonian:")
        #print(HH)        

        # Storing Hamiltonian and dipole moment matrices
        self.HH = HH
        # Hamiltonian operator
        self.HamOp = Hamiltonian(data=HH)
        # dipole moments
        self.DD = DD

        # FIXME: make this on-demand (if poissible)
        trdata = numpy.zeros((DD.shape[0],DD.shape[1],DD.shape[2]),dtype=qr.REAL)
        trdata[:,:,:] = DD[:,:,:]
        self.TrDMOp = TransitionDipoleMoment(data=trdata)

        # Franck-Condon factors
        self.FCf = FC
        # widths
        self.Wd = Wd
        # dephasings
        self.Dr = Dr

        #
        ###3
        #
        #print("Wd:")
        #print(self.Wd) 

        # composition of two-ex states
        # first index of state a is twoex_indx[0, a]
        self.twoex_indx = twoex_indx

        # squares of transition dipoles
        dd2 = numpy.zeros((Ntot, Ntot),dtype=numpy.float64)
        for a in range(Ntot):
            for b in range(Ntot):
                dd2[a,b] = numpy.dot(self.DD[a,b,:],self.DD[a,b,:])
        self.D2 = dd2
        # FIXME: do I need this??? Is it even corrrect??? (maybe amax?)
        # maximum of transition dipole moment elements
        self.D2_max = numpy.max(dd2)

        # Number of states in individual bands
        self.Nb = numpy.zeros(self.mult+1, dtype=numpy.int)
        for ii in range(self.mult+1):
            self.Nb[ii] = self.number_of_states_in_band(band=ii,
            vibgen_approx=vibgen_approx, Nvib=Nvib,
            vibenergy_cutoff=vibenergy_cutoff)

        # Number of electronic states in individual bands
        self.Nbe = numpy.zeros(self.mult+1, dtype=numpy.int)
        for ii in range(self.mult+1):
            self.Nbe[ii] = self.number_of_electronic_states_in_band(band=ii)

        # prepare RWA indices and set info for Rotating Wave Approximation
        rwa_indices = numpy.zeros(self.mult+1, numpy.int)
        for ii in range(self.mult):
            rwa_indices[ii+1] = rwa_indices[ii]+self.Nb[ii]
        self.HamOp.set_rwa(rwa_indices)


        #######################################################################
        #
        # System-bath interaction
        #
        #######################################################################

        #
        #  There are two methods to set system-bath interaction
        #      1) Each molecule gets its bath correlation function
        #      2) Global energy gap correlation function matrix is set
        #

        # is energy gap correlation function matrix present?
        if self._has_egcf_matrix:

            # Check the consistency of the energy gap correlation matrix
            if self.egcf_matrix.nob != self.nmono:
                raise Exception("Correlation matrix has a size different" +
                                " from the number of monomers")
            #FIXME The aggregate having a egcf matrix does not mean the monomers
            #have egcf matrices. They could just have correlation funtions.
            for i in range(self.nmono):
                if self.monomers[i]._is_mapped_on_egcf_matrix and \
                not (self.monomers[i].egcf_matrix is self.egcf_matrix):
                    raise Exception("Correlation matrix in the monomer" +
                                    " has to be the same as the one of" +
                                    " the aggregate.")
            # seems like everything is consistent -> we can calculate system-
            # -bath interaction
            self._has_system_bath_interaction = True

        # if not, try to get one from monomers later
        else:

            self._has_system_bath_interaction = False

        # try to set energy gap correlation matrix from monomers
        if not self._has_system_bath_interaction:

            # let's assume we can calculate EGCF matrix from monomers
            egcf_ok = True

            # get correlation function from a monomer
            try:
                egcf1 = self.monomers[0].get_transition_environment((0,1))
            except:
                # we cannot calculate EGCF matrix, there is no system-bath
                # interaction, or it is not based on correlation functions
                egcf_ok = False

            # if we have correlation functions for nonomers, let's construct
            # the EGCF matrix
            if egcf_ok:
                # time axis of the first monomer
                time = egcf1.axis
                # Number of correlation functions is the number of electronic
                # states minus ground state (this assumes that only electronic
                # states are coupled to the bath)
                Nelg = 1  # ASSUMPTION: here we assume a single electronic
                          # ground state
                if sbi_for_higher_ex:
                    # except for ground state, all electronic states have EGCF
                    Ncf = self.Nel - Nelg
                else:
                    # in single exciton, two-level molecule picture, there is
                    # a single correlation function per monomer
                    # ASSUMPTION: Two-level molecules
                    Ncf = self.nmono

                # instantiate the EGCF matrix object
                self.egcf_matrix = CorrelationFunctionMatrix(time, Ncf)

                # run over all electronic states
                for i in range(self.Nel):

                    # in single exciton band
                    if self.which_band[i] == 1:
                        j = i - Nelg
                        mon = self.monomers[j]
                        # get correlation for a monomer
                        # ASSUMPTION: Two-level molecule
                        cfce = mon.get_transition_environment((0,1))
                        # set correlation function into the diagonal of the
                        # EGCF matrix. Index corresponds to the monomer
                        mapi = self.egcf_matrix.set_correlation_function(cfce,
                                                                     [(j,j)])
                        # FIXME: what is returned?
                        if mapi <= 0:
                            raise Exception("Something's wrong")

                    # in two-exciton band
                    elif (self.which_band[i] == 2) and sbi_for_higher_ex:
                        l = i - Nelg
                        # monomers of a two-exciton state are obtaines
                        # FIXME: is this correct???
                        j = self.elsigs[i][0]
                        k = self.elsigs[i][1]
                        mon1 = self.monomers[j]
                        mon2 = self.monomers[k]
                        # we get correlation functions of the two monomers
                        # ASSUMPTION: Two-level molecules
                        cfce1 = mon1.get_transition_environment((0,1))
                        cfce2 = mon2.get_transition_environment((0,1))
                        # correlation functions are added to form a two-exciton
                        # correlation function
                        cfce = cfce1 + cfce2
                        # Two-exciton correlation function is set into
                        # EGCF matrix
                        mapi = self.egcf_matrix.set_correlation_function(cfce,
                                                                     [(l,l)])

                        # FIXME: cross-correlation between double excitons
                        # needs to be handled.

                        if mapi <= 0:
                            raise Exception("Something's wrong")

                    # some effective theory here
                    # FIXME: make sure we know what the ``sbi_for_higher_ex``
                    #        switch actually means
                    elif (self.which_band[i] == 2) and (not sbi_for_higher_ex):
                        # this should be handled by
                        # a map between double excitons and site cor. functions
                        pass

                    # no theory for higher bands so-far
                    elif (self.which_band[i] > 2) and sbi_for_higher_ex:
                        pass

                self._has_system_bath_interaction = True
                self._has_egcf_matrix = True

        # if all needed for system-bath interaction is present
        # we can construct the SystemBathInteraction object
        if self._has_system_bath_interaction:

            # system interaction operators
            iops = []

            # how many operators should be created
            if sbi_for_higher_ex:
                Nop = self.Nel-1 # all electronic states
            else:
                Nop = self.Nbe[1] # we count only single excited states

            # if there are more states in the single exciton block
            # than the number of sites, it means we have vibrational states
            if self.nmono != self.Nb[1]:
                # create a projection operator for each monomer
                # a monomer corresponds to one single excited state starting
                # with electronic index 1 (0 is the ground state)
                # ASSUMPTION: Two-level molecules
                for i in range(1, Nop+1):
                    op1 = Operator(dim=self.HH.shape[0],real=True)
                    # here we make a projector on a given electronic state |i>
                    # ASSUMPTION: Oscillator is represented by its eigenstates
                    for j in self.vibindices[i]:
                        op1.data[j,j] = 1.0
                    iops.append(op1)

            # standard case with only electronic states
            else:
                # create a projection operator for each monomer
                # a monomer corresponds to one single excited state starting
                # with electronic index 1 (0 is the ground state)
                # ASSUMPTION: Two-level molecules
                for i in range(1, Nop+1):
                    op1 = Operator(dim=self.HH.shape[0],real=True)
                    op1.data[i,i] = 1.0
                    iops.append(op1)

            # we create SystemBathInteraction object
            self.sbi = SystemBathInteraction(iops,
                                self.egcf_matrix, system=self)

        else:
            # system-bath interaction is not present
            pass

        self._built = True

        manager.unset_current_units("energy")


    def rebuild(self, mult=1, sbi_for_higher_ex=False,
              vibgen_approx=None, Nvib=None, vibenergy_cutoff=None):
        """Cleans the object and rebuilds it

        """
        self.clean()
        self.build(mult=mult, sbi_for_higher_ex=sbi_for_higher_ex,
              vibgen_approx=vibgen_approx, Nvib=Nvib,
              vibenergy_cutoff=vibenergy_cutoff)



    ###########################################################################
    #
    #    POST BUILDING METHODS
    #
    ###########################################################################

    def trace_over_vibrations(self, operator, Nt=None):
        """Average an operator over vibrational degrees of freedom

        Average MUST be done in site basis. Only in site basis
        we can distinguish the vibrational states properly

        """
        n_indices = 2
        evolution = False
        whole = False

        if operator.dim == self.Ntot:

            if isinstance(operator, ReducedDensityMatrix) or \
               isinstance(operator, DensityMatrix):

                nop = ReducedDensityMatrix(dim=self.Nel)


            elif isinstance(operator, ReducedDensityMatrixEvolution) or \
               isinstance(operator, DensityMatrixEvolution):

                if Nt is not None:
                    nop = ReducedDensityMatrix(dim=self.Nel)
                    evolution = True
                    whole = False
                else:
                    nop = ReducedDensityMatrixEvolution(operator.TimeAxis)
                    rhoi = ReducedDensityMatrix(dim=self.Nel)
                    # we set zero initial condition, because this initialized
                    # the data storage
                    nop.set_initial_condition(rhoi)
                    evolution = True
                    whole = True

            else:
                raise Exception("Operation not implemented for this type: "+
                                operator.__class__.__name__)

            if n_indices == 2:

                # convert to representation by ground-state oscillator

                # FIXME: This limitation might not be necessary
                # in the ground states of all monomers, there must be the same
                # or greater number of levels than in the excited state

                # over all monomers
                for k in range(self.nmono):
                    mono = self.monomers[k]
                    # over all modes
                    n_mod = mono.get_number_of_modes()
                    for i in range(n_mod):
                        mod = mono.get_Mode(i)
                        n_g = mod.get_nmax(0)
                        # over all excited states
                        # FIXME: this should be mono.Nel as in Aggregate
                        for j in range(mono.nel):
                            if (j > 0):
                                n_e = mod.get_nmax(j)
                                if n_e > n_g:
                                    raise Exception("Number of levels"+
                        " in the excited state of a molecule has to be \n"+
                        "the same or smaller than in the ground state")


                # do the conversion

                #
                # ground state vibrational states
                #
                stgs = []
                for i_g in self.vibindices[0]:
                    vs_g = self.vibsigs[i_g]
                    stg = self.get_VibronicState(vs_g[0],
                                                vs_g[1])
                    stgs.append(stg)

                FcProd = numpy.zeros_like(self.FCf)
                for i in range(FcProd.shape[0]):
                    for j in range(FcProd.shape[1]):
                        for i_g in range(self.Nb[0]):
                            FcProd[i, j] += self.FCf[i_g, i]*self.FCf[j, i_g]

                if evolution:
                    if whole:
                        # loop over electronic states n, m
                        for n in range(self.Nel):
                            for i_n in self.vibindices[n]:
                                for m in range(self.Nel):
                                    for i_m in self.vibindices[m]:
                                        nop._data[:, n, m] += \
                                            operator._data[:, i_n, i_m]*FcProd[i_n, i_m]

                    else:
                        # loop over electronic states n, m
                        for n in range(self.Nel):
                            for i_n in self.vibindices[n]:
                                for m in range(self.Nel):
                                    for i_m in self.vibindices[m]:
                                        nop._data[n,m] += \
                                            operator._data[Nt, i_n, i_m]*FcProd[i_n, i_m]

                else:
                    # loop over electronic states n, m
                    for n in range(self.Nel):
                        for i_n in self.vibindices[n]:
                            for m in range(self.Nel):
                                for i_m in self.vibindices[m]:
                                    nop._data[n,m] += \
                                        operator._data[i_n, i_m]*FcProd[i_n, i_m]

            else:
                raise Exception("Cannot trace over this object: "+
                                "wrong number of indices")

            return nop

        else:
            raise Exception("Incompatible operator")


    def convert_to_DensityMatrix(self, psi, trace_over_vibrations=True):
        """Converts StateVector into DensityMatrix (possibly reduced one)

        """

        if trace_over_vibrations:

            if isinstance(psi, StateVector):
                rho = psi.get_DensityMatrix()
                rho = self.trace_over_vibrations()
            elif isinstance(psi, StateVectorEvolution):
                # FIXME: Implement direct conversion
                rho = psi.get_DensityMatrixEvolution()
                rho = self.trace_over_vibrations()

        else:

            if isinstance(psi, StateVector):
                rho = psi.get_DensityMatrix()
            elif isinstance(psi, StateVectorEvolution):
                rho = psi.get_DensityMatrixEvolution()

        return rho


    def get_RWA_suggestion(self):
        """Returns average transition energy

        Average transition energy of the monomer as a suggestion for
        RWA frequency

        """

        Nn = self.Nb[1]  # number of monomers
        esum = 0.0
        for i in range(Nn):
            mn = self.monomers[i]
            omeg = mn.get_energy(1) - mn.get_energy(0)
            esum += omeg

        return esum/Nn


    def get_RelaxationTensor(self, timeaxis,
                       relaxation_theory=None,
                       time_dependent=False,
                       secular_relaxation=False,
                       relaxation_cutoff_time=None,
                       coupling_cutoff=None,
                       recalculate=True,
                       as_operators=False):
        """Returns a relaxation tensor corresponding to the aggregate


        Parameters
        ----------

        timeaxis : TimeAxis
            Time axis of the relaxation tensor calculation. It has to be
            compatible with the time axis of the correlation functions

        relaxation_theory: str
            One of the available relaxation theories

        time_dependent : boolean
            Should the relaxation tensor time dependent?

        secular_relaxation :
            Should the tensor be secular?


        Returns
        -------

        RR : RelaxationTensor
            Relaxation tensor of the aggregate

        ham : Hamiltonian
            Hamiltonian corresponding to the aggregate, renormalized by
            the system-bath interaction


        """

        from ..qm import RedfieldRelaxationTensor
        from ..qm import TDRedfieldRelaxationTensor
        from ..qm import FoersterRelaxationTensor
        from ..qm import TDFoersterRelaxationTensor
        from ..qm import RedfieldFoersterRelaxationTensor
        from ..qm import TDRedfieldFoersterRelaxationTensor
        from ..qm import LindbladForm

        from ..core.managers import eigenbasis_of

        if self._built:
            ham = self.get_Hamiltonian()
            sbi = self.get_SystemBathInteraction()
        else:
            raise Exception()

        #
        # Dictionary of available theories
        #
        theories = dict()
        theories[""] = [""]
        theories["standard_Redfield"] = ["standard_Redfield","stR","Redfield",
                                         "CLME2","QME"]
        theories["standard_Foerster"] = ["standard_Foerster","stF","Foerster"]
        theories["combined_RedfieldFoerster"] = ["combined_RedfieldFoerster",
                                                 "cRF","Redfield-Foerster"]
        #
        # Future
        #
        theories["modified_Redfield"] = ["modifield_Redfield", "mR"]
        theories["noneq_modified_Redfield"] = ["noneq_modified_Redfield",
                                               "nemR"]
        theories["generalized_Foerster"] = ["generalized_Foerster", "gF",
                                            "multichromophoric_Foerster"]
        theories["noneq_Foerster"] = ["noneq_Foerster", "neF"]
        theories["combined_WeakStrong"] = ["combined_WeakStrong", "cWS"]
        theories["Lindblad_form"] = ["Lindblad_form", "Lf"]
        theories["electronic_Lindblad"] = ["electronic_Lindblad", "eLf"]

        #if ((not recalculate) and
        #    (relaxation_theory in theories[self._relaxation_theory])):
        #    return self.RelaxationTensor, self.RelaxationHamiltonian


        if relaxation_theory in theories["standard_Redfield"]:

            if time_dependent:

                # Time dependent standard Refield

                ham.protect_basis()
                with eigenbasis_of(ham):
                    relaxT = TDRedfieldRelaxationTensor(ham, sbi,
                                        cutoff_time=relaxation_cutoff_time,
                                        as_operators=as_operators)
                    if secular_relaxation:
                        relaxT.secularize()
                ham.unprotect_basis()

            else:

                # Time independent standard Refield


                ham.protect_basis()

                with eigenbasis_of(ham):
                    relaxT = RedfieldRelaxationTensor(ham, sbi,
                                                    as_operators=as_operators)

                    if secular_relaxation:
                        relaxT.secularize()

                ham.unprotect_basis()


            self.RelaxationTensor = relaxT
            self.RelaxationHamiltonian = ham
            self._has_relaxation_tensor = True
            self._relaxation_theory = "standard_Redfield"

            return relaxT, ham

        elif relaxation_theory in theories["standard_Foerster"]:

            if time_dependent:

                # Time dependent standard Foerster
                relaxT = TDFoersterRelaxationTensor(ham, sbi)
                dat = numpy.zeros((ham.dim,ham.dim),dtype=numpy.float64)
                for i in range(ham.dim):
                    dat[i,i] = ham._data[i,i]
                ham_0 = Hamiltonian(data=dat)

            else:

                # Time independent standard Foerster

                #
                # This is done strictly in site basis
                #

                relaxT = FoersterRelaxationTensor(ham, sbi)
                dat = numpy.zeros((ham.dim,ham.dim),dtype=numpy.float64)
                for i in range(ham.dim):
                    dat[i,i] = ham._data[i,i]
                ham_0 = Hamiltonian(data=dat)

            # The Hamiltonian for propagation is the one without
            # resonance coupling
            self.RelaxationTensor = relaxT
            self.RelaxationHamiltonian = ham_0
            self._has_relaxation_tensor = True
            self._relaxation_theory = "standard_Foerster"

            return relaxT, ham_0

        elif relaxation_theory in theories["combined_RedfieldFoerster"]:

            if time_dependent:

                # Time dependent combined tensor
                ham.subtract_cutoff_coupling(coupling_cutoff)
                ham.protect_basis()
                with eigenbasis_of(ham):
                    relaxT = \
                             TDRedfieldFoersterRelaxationTensor(ham, sbi,
                                            coupling_cutoff=coupling_cutoff,
                                            cutoff_time=relaxation_cutoff_time)
                    if secular_relaxation:
                        relaxT.secularize()
                ham.unprotect_basis()
                ham.recover_cutoff_coupling()

            else:

                # Time independent combined tensor
                ham.subtract_cutoff_coupling(coupling_cutoff)
                ham.protect_basis()
                with eigenbasis_of(ham):
                    relaxT = \
                             RedfieldFoersterRelaxationTensor(ham, sbi,
                                            coupling_cutoff=coupling_cutoff,
                                            cutoff_time=relaxation_cutoff_time)
                    if secular_relaxation:
                        relaxT.secularize()

                    #print("Last line of the context", Manager().get_current_basis())
                #print("Left context", Manager().get_current_basis())
                ham.unprotect_basis()
                ham.recover_cutoff_coupling()

            #
            # create a corresponding propagator
            #
            ham1 = Hamiltonian(data=ham.data.copy())
            #ham1.subtract_cutoff_coupling(coupling_cutoff)
            ham1.remove_cutoff_coupling(coupling_cutoff)

            self.RelaxationTensor = relaxT
            self.RelaxationHamiltonian = ham1
            self._has_relaxation_tensor = True
            self._relaxation_theory = "combined_RedfieldFoerster"

            return relaxT, ham1

        elif relaxation_theory in theories["combined_WeakStrong"]:

            pass

        elif relaxation_theory in theories["Lindblad_form"]:

            if time_dependent:

                # Time dependent standard Refield
                raise Exception("Time dependent Lindblad not implemented yet")

            else:

                # Linblad form

                #ham.protect_basis()
                #with eigenbasis_of(ham):
                relaxT = LindbladForm(ham, sbi)
                if secular_relaxation:
                    relaxT.convert_2_tensor()
                    relaxT.secularize()
                #ham.unprotect_basis()

            self.RelaxationTensor = relaxT
            self.RelaxationHamiltonian = ham
            self._has_relaxation_tensor = True
            self._relaxation_theory = "Lindblad_form"

            return relaxT, ham

        elif relaxation_theory in theories["electronic_Lindblad"]:

            if time_dependent:

                # Time dependent standard Refield
                raise Exception("Time dependent Lindblad not implemented yet")

            else:

                # For purely electronic system, calculate normal Lindblad form
                if self.Ntot == self.Nel:
                    relaxT, ham = self.get_RelaxationTensor(timeaxis,
                                relaxation_theory="Lindblad_form",
                                time_dependent=time_dependent,
                                secular_relaxation=secular_relaxation,
                                relaxation_cutoff_time=relaxation_cutoff_time,
                                coupling_cutoff=coupling_cutoff,
                                recalculate=recalculate)
                # if vibrational states are present, we create a new SBI
                else:
                    # we assume that we have only electronic sbi
                    # FIXME: make sure that Molecule also has Nel
                    if sbi.system.Nel == sbi.KK.shape[1]:
                        # upgrade sbi to vibrational levels

                        eKK = sbi.KK
                        vKK = numpy.zeros((sbi.KK.shape[0], ham.dim, ham.dim),
                                          dtype=numpy.float64)

                        # use eKK to calculate vKK

                        sbi.KK = vKK
                    else:
                        raise Exception("SystemBathInteraction object has to"+
                                        " purely electronic")

                    relaxT = LindbladForm(ham, sbi)

                if secular_relaxation:
                    relaxT.convert_2_tensor()
                    relaxT.secularize()

            self.RelaxationTensor = relaxT
            self.RelaxationHamiltonian = ham
            self._has_relaxation_tensor = True
            self._relaxation_theory = "Lindblad_form"

            return relaxT, ham


        else:

            raise Exception("Theory not implemented")


    def get_ReducedDensityMatrixPropagator(self, timeaxis,
                       relaxation_theory=None,
                       time_dependent=False,
                       secular_relaxation=False,
                       relaxation_cutoff_time=None,
                       coupling_cutoff=None,
                       as_operators=False,
                       recalculate=True):
        """Returns propagator of the density matrix



        """


        from ..qm import ReducedDensityMatrixPropagator
        from ..core.managers import eigenbasis_of


        relaxT, ham = self.get_RelaxationTensor(timeaxis,
                       relaxation_theory=relaxation_theory,
                       time_dependent=time_dependent,
                       secular_relaxation=secular_relaxation,
                       relaxation_cutoff_time=relaxation_cutoff_time,
                       coupling_cutoff=coupling_cutoff,
                       recalculate=recalculate,
                       as_operators=as_operators)

        with eigenbasis_of(ham):
            prop = ReducedDensityMatrixPropagator(timeaxis, ham, relaxT)


        return prop


    #FIXME: There must be a general theory here
    def get_RedfieldRateMatrix(self):

        from ..qm import RedfieldRateMatrix
        from ..core.managers import eigenbasis_of

        if self._built:
            ham = self.get_Hamiltonian()
            sbi = self.get_SystemBathInteraction()
        else:
            raise Exception()

        ham.protect_basis()
        with eigenbasis_of(ham):
            RR = RedfieldRateMatrix(ham, sbi)
        ham.unprotect_basis()

        return RR
    
    
    def get_FoersterRateMatrix(self):
        
        from ..qm import FoersterRateMatrix
        
        if self._built:        
            ham = self.get_Hamiltonian()
            sbi = self.get_SystemBathInteraction()
        else:
            raise Exception()

        return FoersterRateMatrix(ham, sbi)


    def diagonalize(self):
        """Transforms some internal quantities into diagonal basis

        """

        if self._diagonalized:
            return

        #
        ###3
        #
        #for i in range(self.HH.shape[0]):
        #    print(self.HH[i,i])

        ee,SS = numpy.linalg.eigh(self.HH)
        
        self.Hs = self.HH.copy()

        self.HD = ee
        self.SS = SS
        self.S1 = numpy.linalg.inv(SS)

        self.HH = numpy.dot(self.S1,numpy.dot(self.HH,self.SS))

       
        have_vibs = False
        if len(self.vibindices[0]) > 1:
             have_vibs = True
 
        # Kronecker delta over all states
        delta = operator_factory(self.Ntot).unity_operator()  
          
        #have_vibs = True
        
        if not have_vibs: 

            ######################################################################
            #  CASE OF NO VIBRATIONAL MODES
            ######################################################################
    
            #
            # some quantities to be precalculated for two-ex lineshape
            # 1->2 has to be trasformed first because we need untransformed 0->1
            # for such a transformation
            #
            N1b = self.Nb[0]+self.Nb[1]
    
    
            #  \kappa_{nA} =
            #  \sum_{K}(\delta_{nk}+\delta_{nl})*|\langle A | K\rangle|^2
            #
            #  where K is a two-exc. state K = (k,l), A is a two-ex. state
            #  and n is a single exciton state
            #
            #  Below aa1 = A, aa2 = n, aa3 = K, st_k = k and st_l = l
            #
            kappa = numpy.zeros((self.Ntot, self.Ntot), dtype=qr.REAL)

    
            if self.mult >= 2:
    
                N2b = self.Nb[0]+self.Nb[1]+self.Nb[2]
    
                # all states (and 2-ex band selected)
                for el1 in range(self.Nel):
                    if self.which_band[el1] == 2:
                        # all states corresponding to electronic two-exc. state kk
                        for aa1 in self.vibindices[el1]:
    
                            # all states and (1-ex band selected)
                            for el2 in range(self.Nel):
                                if self.which_band[el2] == 1:
                                    for aa2 in self.vibindices[el2]:
    
                                        # all states and (2-ex band selected)
                                        for el3 in range(self.Nel):
                                            if self.which_band[el3] == 2:
                                                for aa3 in self.vibindices[el3]:
    
                                                    st_k = self.twoex_indx[aa3,0]
                                                    st_l = self.twoex_indx[aa3,1]
    
                                                    kappa[aa2, aa1] += (
                                                         (delta[aa2, st_k]
                                                        + delta[aa2, st_l])*
                                                         (SS[aa3, aa1]**2))
    
                #
                # Cross terms
                #
                for aa_2x in range(N1b, N2b):
                    for alpha in range(N1b):
                        self.Wd[aa_2x, alpha] = 0.0
                        for nn_2x in range(N1b, N2b):
                            for k_1x in range(N1b):
                                st_n = self.twoex_indx[nn_2x, 0]
                                st_m = self.twoex_indx[nn_2x, 1]
                                #print(st_n, st_m)
                                self.Wd[aa_2x, alpha] += \
                                    ((self.Wd[st_n, st_n]**2)*delta[st_n, k_1x] +
                                     (self.Wd[st_m, st_m]**2)*delta[st_m, k_1x])*\
                                     (SS[nn_2x, aa_2x]**2)*(SS[k_1x, alpha]**2)
    
                self.Wd[N1b:N2b,0:N1b] = numpy.sqrt(self.Wd[N1b:N2b,0:N1b])
                self.Wd[0:N1b,N1b:N2b] = numpy.transpose(self.Wd[N1b:N2b,0:N1b])
    
                #
                # Transform line shapes for 1->2 transitions
                #
                Wd_a = numpy.zeros(N2b, dtype=qr.REAL)
                Dr_a = numpy.zeros(N2b, dtype=qr.REAL)

    
                for aa in range(N1b, N2b):
                    for nn in range(N1b, N2b):
                        st_n = self.twoex_indx[nn,0]
                        st_m = self.twoex_indx[nn,1]
                        Wd_a[aa] += (SS[nn, aa]**2)*\
                                    ((self.Wd[st_n, st_n]**2)*kappa[st_n, aa]
                                    +(self.Wd[st_m, st_m]**2)*kappa[st_m, aa])
                             
                        
                W_aux = numpy.diag(numpy.sqrt(Wd_a))
                #W_aux = numpy.diag(numpy.sqrt(Wd_b))
                self.Wd[N1b:N2b,N1b:N2b] = W_aux[N1b:N2b,N1b:N2b]
    
            #
            # Transform line shapes for 0->1 transitions
            #
            Wd_a = numpy.zeros(N1b, dtype=qr.REAL)
            Dr_a = numpy.zeros(N1b, dtype=qr.REAL)
            for ii in range(N1b):
                for nn in range(N1b):
                    Wd_a[ii] += (self.Wd[nn,nn]**2)*abs(SS[ii,nn])**4
                    Dr_a[ii] += (self.Dr[nn,nn]**2)*abs(SS[ii,nn])**4
            Wd_a = numpy.sqrt(Wd_a)
            Dr_a = numpy.sqrt(Dr_a)
    
            self.Wd[0:N1b,0:N1b] = numpy.diag(Wd_a)
            self.Dr[0:N1b,0:N1b] = numpy.diag(Dr_a)
            
            
            #
            ###3
            #
            #print("========")
            #for k in range(N1b):
            #    print(Wd_a[k])
    
            #raise Exception()

            #print(self.Wd)

        else:


            ######################################################################
            # CASE OF VIBRATIONAL MODES
            ######################################################################

            N1b = self.Nb[0]+self.Nb[1]            


            #
            # Transform line shapes for 0->1 transitions
            #
            Wd_a = numpy.zeros(N1b, dtype=qr.REAL)
            Dr_a = numpy.zeros(N1b, dtype=qr.REAL)
            
            Nel = 1 + self.nmono
            
            Wd_in = numpy.zeros(Nel, dtype=qr.REAL)
            Dr_in = numpy.zeros(Nel, dtype=qr.REAL)
            
            Wd_ini = self.Wd.copy()
            Dr_ini = self.Dr.copy()
            
            for ii in range(Nel):
                for k in self.vibindices[ii]:
                    Wd_in[ii] = Wd_ini[k,k] #self.Wd[k,k]
                    #print("***", self.Wd[k,k], self.Hs[k,k])
                    Dr_in[ii] = Dr_ini[k,k]
            
            #Nvib1el = len(self.vibindices[0])
            kap = numpy.zeros((N1b, Nel), dtype=qr.REAL)
            
            # loop over all states in the 1-ex band
            for aa in range(N1b):
                ela = self.elinds[aa]
                if self.which_band[ela] == 1:
                    
                    # loop over electronic states in the 1-ex band
                    st = 0  # counts the total index of the state
                    for ii in range(Nel):
                        #if self.which_band[ii] == 1:
                         if True:   
                            # loop over substructure of vib states
                            for ialph in self.vibindices[ii]:
                                #print(aa, st, "(", ii, ialph,")")
                                kap[aa, ii] += numpy.abs(SS[st,aa])**2
                                st += 1
                        #else:
                        #    for ialph in self.vibindices[ii]:
                        #        st += 1
            
            # loop over all states in the 1-ex band

            for aa in range(N1b):                    
                for nn in range(Nel):
                    
                    Wd_a[aa] += (kap[aa,nn]**2)*(Wd_in[nn]**2)    
                    Dr_a[aa] += (kap[aa,nn]**2)*(Dr_in[nn]**2)
                    
                    
            Wd_a = numpy.sqrt(Wd_a)
            Dr_a = numpy.sqrt(Dr_a)
    
            self.Wd[0:N1b,0:N1b] = numpy.diag(Wd_a)
            self.Dr[0:N1b,0:N1b] = numpy.diag(Dr_a)
            
            #print("First version")
            #print(self.Wd)

            if self.mult >= 2:
                
                Nel = self.Nel
                N2b = self.Ntot
                
                Wd_b = numpy.zeros(N1b, dtype=qr.REAL)
                
                Dr_b = numpy.zeros(N1b, dtype=qr.REAL)
                Wd_in = numpy.zeros(Nel)
                Dr_in = numpy.zeros(Nel)
                
                
                for ii in range(Nel):
                    for k in self.vibindices[ii]:
                        Wd_in[ii] = Wd_ini[k,k] #self.Wd[k,k]
                        Dr_in[ii] = Dr_ini[k,k]
                    
                kap2 = numpy.zeros((N2b, Nel), dtype=qr.REAL)
                # loop over all states in the 1-ex band
                for aa in range(N2b):
                    ela = self.elinds[aa]
                    #if self.which_band[ela] == 1:
                    if True: 
                        # loop over electronic states in the 1-ex band
                        st = 0  # counts the total index of the state
                        for ii in range(Nel):
                            #if self.which_band[ii] == 1:
                            #print("el. state: ", ii)
                            if True:   
                                # loop over substructure of vib states
                                for ialph in self.vibindices[ii]:
                                    #print(aa, st, "(", ii, ialph,")")
                                    kap2[aa, ii] += numpy.abs(SS[st,aa])**2
                                    st += 1
                            
                for aa in range(N1b):                    
                    for nn in range(Nel):
                    
                        Wd_b[aa] += (kap2[aa,nn]**2)*(Wd_in[nn]**2)    
                        Dr_b[aa] += (kap2[aa,nn]**2)*(Dr_in[nn]**2)               
                
                Wd_b = numpy.sqrt(Wd_b)
                Dr_b = numpy.sqrt(Dr_b)
                
                #
                # Single exciton band
                #
                self.Wd[0:N1b,0:N1b] = numpy.diag(Wd_b)
                self.Dr[0:N1b,0:N1b] = numpy.diag(Dr_b)
                
                
                Wd_c = numpy.zeros((self.Ntot, self.Ntot), dtype=qr.REAL)

                for aa in range(N1b, N2b):
                    for bb in range(N1b, N2b):
                        
                        for nn in range(Nel):
                            vind = self.vibindices[nn]
                            nni = vind[0]
                            n = self.twoex_indx[nni,0]
                            m = self.twoex_indx[nni,1]
                            for mm in range(Nel):
                                vind = self.vibindices[mm]
                                mmi = vind[0]
                                k = self.twoex_indx[mmi,0]
                                l = self.twoex_indx[mmi,1]
                        
                                Wd_c[aa,bb] += ((Wd_in[n]**2)*(delta[n,k]+delta[n,l])
                                               +(Wd_in[m]**2)*(delta[m,k]+delta[m,l]))\
                                              *kap2[aa,nn]*kap2[bb,mm]

                W_cc = numpy.zeros(Wd_c.shape[0], dtype=qr.REAL)
                for k in range(N2b):
                    W_cc[k] = Wd_c[k,k]
                W_aux = numpy.diag(numpy.sqrt(W_cc))
                #W_aux = numpy.diag(numpy.sqrt(Wd_b))
                
                #
                #  Two-exciton band
                #
                self.Wd[N1b:N2b,N1b:N2b] = W_aux[N1b:N2b,N1b:N2b]

                #print("Second version")
                #print(self.Wd)
                #raise Exception()
                
                Wd_c = numpy.zeros((self.Ntot, self.Ntot), dtype=qr.REAL)
                
                for aa in range(N1b, N2b):
                    for bb in range(N1b):
                        
                        for nn in range(Nel):
                            vind = self.vibindices[nn]
                            nni = vind[0]
                            n = self.twoex_indx[nni,0]
                            m = self.twoex_indx[nni,1]
                            for k in range(1+self.nmono):
                                Wd_c[aa,bb] += ((Wd_in[n]**2)*delta[n,k] 
                                               +(Wd_in[m]**2)*delta[m,k]) \
                                              *kap2[aa,nn]*kap2[bb,k] 
                    
                W_aux = numpy.sqrt(Wd_c)
                
                self.Wd[N1b:N2b,0:N1b] = W_aux[N1b:N2b,0:N1b]
                self.Wd[0:N1b,N1b:N2b] = numpy.transpose(self.Wd[N1b:N2b,0:N1b])
                

        #
        #
        # Coefficients xi_{ai} to transform pure dephasing of electronic coherence
        #
        #
        # all states and (1-ex band selected)
        for aa in range(N1b):
            
            for ii in range(self.Nel):
                #el2 = self.elinds[ii]
                if self.which_band[ii] == 1:
                    for ialph in self.vibindices[ii]:
                        self.Xi[aa, ii] += SS[ialph,aa]**2
                        


        #
        # Transform transition dipole moments
        #
        for n in range(3):
            self.DD[:,:,n] = numpy.dot(self.S1,
                               numpy.dot(self.DD[:,:,n],self.SS))

        Ntot = self.HH.shape[0]
        dd2 = numpy.zeros((Ntot,Ntot),dtype=numpy.float64)
        for a in range(Ntot):
            for b in range(Ntot):
                dd2[a,b] = numpy.dot(self.DD[a,b,:],self.DD[a,b,:])

        self.D2 = dd2
        self.D2_max = numpy.max(dd2)

        self.rho0 = numpy.zeros(self.HH.shape, dtype=qr.COMPLEX)
        self.rho0[0,0] = 1.0

        self._diagonalized = True


    def _thermal_population(self, temp=0.0, subtract=None,
                            relaxation_hamiltonian=None, start=0):
        """Thermal populations at temperature temp

        Thermal populations calculated from the diagonal elements
        of the Hamiltonian.

        Parameters
        ----------

        temp : float
            Temperature in Kelvins

        subtract : list like
            Reoreganization energies to subtract from the Hamiltonian

        relaxation_hamiltonian: array
            Hamiltonian according to which we form thermal equilibrium

        """

        from ..core.units import kB_intK

        kBT = kB_intK*temp

        #if not relaxation_hamiltonian:
        #    HH = self.get_Hamiltonian()
        #else:
        #    HH = relaxation_hamiltonian
        HH = relaxation_hamiltonian

        # This is all done with arrays, not with Qrhei objects
        #HH = HH.data
        dim = HH.shape[0]

        if subtract is None:
            subtract = numpy.zeros(dim, dtype=numpy.float64)

        rho0 = numpy.zeros((dim, dim),dtype=numpy.complex128)


        if temp == 0.0:
            rho0[start,start] = 1.0

        else:
            # FIXME: we assume only single exciton band

            ens = numpy.zeros(dim-start, dtype=numpy.float64)

            # we specify the basis from outside. This allows to choose
            # canonical equilibrium in arbitrary basis
            for i in range(start, dim):
                ens[i-start] = HH[i,i] - subtract[i-start]

            ne = numpy.exp(-ens/kBT)
            sne = numpy.sum(ne)
            rho0_diag = ne/sne
            rho0[start:,start:] = numpy.diag(rho0_diag)


        return rho0


    def _impulsive_population(self, relaxation_theory_limit="weak_coupling",
                              temperature=0.0):
        """Impulsive excitation of the density matrix from ground state

        """

        rho = self.get_DensityMatrix(condition_type="thermal",
                            relaxation_theory_limit=relaxation_theory_limit,
                            temperature=temperature)
        rho0 = rho.data

        DD = self.TrDMOp.data

        # abs value of the transition dipole moment
        dabs = numpy.sqrt(DD[:,:,0]**2 + \
                          DD[:,:,1]**2 + DD[:,:,2]**2)
        # excitation from bra and ket
        rho0 = numpy.dot(dabs, numpy.dot(rho0,dabs))

        return rho0


    def get_DensityMatrix(self, condition_type=None,
                                relaxation_theory_limit="weak_coupling",
                                temperature=None,
                                relaxation_hamiltonian=None):
        """Returns density matrix according to specified condition

        Returs density matrix to be used e.g. as initial condition for
        propagation.

        Parameters
        ----------

        condition_type : str
            Type of the initial condition. If None, the property rho0, which
            was presumably calculated in the past, is returned.

        relaxation_theory_limits : str {weak_coupling, strong_coupling}
            Type of the relaxation theory limits;
            We mean the system bath coupling. When `weak_coupling` is chosen,
            the density matrix is returned in form of a canonical equilibrium
            in terms of the exciton basis. For `strong_coupling`,
            the canonical equilibrium is calculated in site basis with site
            energies striped of reorganization energies.

        temperature : float
            Temperature in Kelvin

        relaxation_hamiltonian :
            Hamiltonian according to which we form thermal equilibrium. In case
            of `strong_coupling`, no reorganization energies are subtracted -
            we assume that the supplied energies are already void of them.

        Condition types
        ---------------

        thermal
            Thermally equilibriated population of the whole density matrix

        thermal_excited_state
            Thermally equilibriuated excited state

        impulsive_excitation
            Excitation by ultrabroad laser pulse

        """

        # aggregate must be built before we call this method
        if not self._built:
            raise Exception("Aggregate must be built before"
                            +" get_DensityMatrix can be invoked.")

        # if Aggregate has interaction with the bath, temperature
        # is already defined
        if temperature is None:
            if self.sbi is None:
                temperature = 0.0
            elif self.sbi.has_temperature():
                temperature = self.sbi.get_temperature()
            else:
                temperature = 0.0

        # if no condition is specified, it is understood that we return
        # internal rho0, which was calculated sometime in the past
        if condition_type is None:
            return DensityMatrix(data=self.rho0)


        # impulsive excitation from a thermal ground state
        elif condition_type == "impulsive_excitation":
            rho0 = self._impulsive_population(
                              relaxation_theory_limit=relaxation_theory_limit,
                              temperature=temperature)
            self.rho0 = rho0
            return DensityMatrix(data=self.rho0)


        # thermal population based on the total Hamiltonian
        elif condition_type == "thermal":

            if not relaxation_hamiltonian:
                Ham = self.get_Hamiltonian()
            else:
                Ham = relaxation_hamiltonian

            # FIXME: weak and strong limits not distinguished
            rho0 = self._thermal_population(temperature,
                                            relaxation_hamiltonian=Ham.data)

            self.rho0 = rho0
            return DensityMatrix(data=self.rho0)

        elif condition_type == "thermal_excited_state":

            if relaxation_theory_limit == "strong_coupling":

                start = self.Nb[0] # this is where excited state starts
                n1ex= self.Nb[1] # number of excited states in one-ex band

                if not relaxation_hamiltonian:
                    HH = self.get_Hamiltonian()
                    Ndim = HH.dim
                    re = numpy.zeros(Ndim-start, dtype=numpy.float64)
                    # we need to subtract reorganization energies
                    for i in range(n1ex):
                        re[i] = \
                        self.sbi.get_reorganization_energy(i)
                else:
                    HH = relaxation_hamiltonian
                    Ndim = HH.dim
                    re = numpy.zeros(Ndim-start, dtype=numpy.float64)
                    # here we assume that reorganizaton energies are already
                    # removed


                # we get this in SITE BASIS
                ham = HH.data

                rho0 = self._thermal_population(temperature,
                                                subtract=re,
                                                relaxation_hamiltonian=ham,
                                                start=start)

            elif relaxation_theory_limit == "weak_coupling":

                if not relaxation_hamiltonian:
                    Ham = self.get_Hamiltonian()
                else:
                    Ham = relaxation_hamiltonian

                # we get this in EXCITON BASIS
                with qr.eigenbasis_of(Ham):
                    H = Ham.data

                start = self.Nb[0] # this is where excited state starts

                # we subtract lowest energy to ease the calcultion,
                # but we do not remove reorganization enegies
                subt = numpy.zeros(H.shape[0])
                subtfil = numpy.amin(numpy.array([H[ii,ii] \
                                    for ii in range(start, H.shape[0])]))
                subt.fill(subtfil)

                rho0 = self._thermal_population(temperature,\
                            subtract = subt,
                            relaxation_hamiltonian=H,
                            start=start)

            else:
                raise Exception("Unknown relaxation_theory_limit")

            self.rho0 = rho0
            return DensityMatrix(data=self.rho0)

        else:
            raise Exception("Unknown condition type")
        #
        # TESTED


    def get_temperature(self):
        """Returns temperature associated with this aggregate


        The temperature originates from the system-bath interaction

        """

        # aggregate must be built before we call this method
        if not self._built:
            raise Exception()

        return self.sbi.CC.get_temperature()
        #
        # TESTED


    def get_electronic_groundstate(self):
        """Indices of states in electronic ground state


        Returns indices of all states in the electronic
        ground state of the system.

        """

        Ng = self.Nb[0]
        lst = [k for k in range(Ng)]

        return tuple(lst)


    def get_excitonic_band(self, band=1):
        """Indices of states in a given excitonic band.


        Returns indices of all states in the excitonic band
        with number of excitons equal to `band`

        Parameters
        ----------

        band : int
            Specifies which band should be returned.

        """
        Nbefore = 0
        for ii in range(band):
            Nbefore += self.Nb[ii]
        Nin = self.Nb[band]
        lst = [k for k in range(Nbefore, Nbefore+Nin)]

        return tuple(lst)


    def get_transition(self, Nf, Ni):
        """Returns relevant info about the energetic transition

        Parameters
        ----------

        Nf : {int, ElectronicState, VibronicState}
            Final state of the transition

        Ni : {int, ElectronicState VibronicState}
            Initial state of the transition

        """
        if (isinstance(Nf, ElectronicState)
            and isinstance(Ni, ElectronicState)):

            if self.Ntot == self.Nel:
                iNf = Nf.index
                iNi = Ni.index
            else:
                raise Exception("The Hamiltonian is not pure electronic")

        elif (isinstance(Nf, VibronicState)
            and isinstance(Ni, VibronicState)):
            vsig = Nf.get_vibsignature()
            esig = Nf.get_ElectronicState().get_signature()
            iNf = self.vibsigs.index((esig, vsig))

            #print(esig, vsig, iNf)
            vsig = Ni.get_vibsignature()
            esig = Ni.get_ElectronicState().get_signature()
            iNi = self.vibsigs.index((esig, vsig))
            #print(esig, vsig, iNi)

        else:
            iNf = Nf
            iNi = Ni

        #
        # if Nf and Ni are not of the same type, it will lead to Exception
        #
        energy = self.convert_energy_2_current_u(self.HH[iNf,iNf]
                                                -self.HH[iNi,iNi])
        trdipm = self.DD[iNf,iNi,:]

        return (energy, trdipm)


    def has_SystemBathInteraction(self):
        """Returns True if the Aggregate is embedded in a defined environment

        """

        # aggregate must be built before we call this method
        if not self._built:
            raise Exception()

        if (self.sbi is not None) and self._has_system_bath_interaction:
            return True

        return False

    def get_SystemBathInteraction(self):
        """Returns the aggregate SystemBathInteraction object

        """
        if self._built:
            return self.sbi
        else:
            raise Exception("Aggregate object not built")


    def set_SystemBathInteraction(self, sbi):
        """Sets the SystemBathInteraction operator for this aggregate

        """
        # FIXME: check its compatibility
        self.sbi = sbi
        self.sbi.set_system(self)



    def get_Hamiltonian(self):
        """Returns the aggregate Hamiltonian

        """
        if self._built:
            return self.HamOp #Hamiltonian(data=self.HH)
        else:
            raise Exception("Aggregate object not built")

    def get_electronic_Hamiltonian(self, full=False):
        """Returns the aggregate electronic Hamiltonian

        In case this is a purely electronic aggregate, the output
        is identical to get_Hamiltonian()

        """

        HH = numpy.zeros((self.Nel, self.Nel), dtype=qr.REAL)
        for (a, sta) in self.elstates(mult=self.mult):
            HH[a,a] = sta.energy()
            for (b, stb) in self.elstates(mult=self.mult):
                if a != b:
                    HH[a,b] = self.coupling(sta, stb, full=full)
        HHel = Hamiltonian(data=HH)

        return HHel


    def get_TransitionDipoleMoment(self):
        """Returns the aggregate transition dipole moment operator

        """
        if self._built:
            return self.TrDMOp # TransitionDipoleMoment(data=self.DD)
        else:
            raise Exception("Aggregate object not built")
