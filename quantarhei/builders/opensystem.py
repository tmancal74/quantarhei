# -*- coding: utf-8 -*-
import numpy

from ..core.managers import Manager
from ..utils import array_property
from ..qm.hilbertspace.hamiltonian import Hamiltonian
from ..qm.liouvillespace.heom import KTHierarchy
from ..qm.liouvillespace.heom import KTHierarchyPropagator
from ..core.managers import eigenbasis_of
from ..core.managers import energy_units
from ..qm import SelfAdjointOperator
from ..qm import ReducedDensityMatrix
from ..core.dfunction import DFunction
from ..core.triangle import triangle
from ..core.units import kB_intK
from .. import REAL
from .. import COMPLEX

from ..core.units import eps0_int, c_int

class OpenSystem():
    """The class representing a general open quantum system
    
    It provides routines to store and extract the most interesting
    characteristics of an open qwuantum system, such as the Hamiltonian,
    SystemBathInteraction object, relaxation tensor and the reduced density
    matrix evolution.
    
    
    """

    # energies of electronic states
    elenergies = array_property('elenergies')
    # transition dipole moments
    dmoments = array_property('dmoments')
    
    def __init__(self, elen):
        
        #
        # Analyze energies
        #
        count = 0
        nband = 0
        which_band = []
        i_en = []
        for sub in elen:
            try:
                for ssub in sub:
                    count += 1
                    i_en.append(ssub)
                    which_band.append(nband)
                nband += 1
            except:
                count += 1
                i_en.append(sub)
                which_band.append(nband)
                nband += 1

        self.which_band = which_band

        #
        #  We have multiplicity of bands
        #
        self.mult = nband - 1

        self.Nb = numpy.zeros(nband, dtype=numpy.int32)
        for bb in range(count):
            self.Nb[self.which_band[bb]] += 1

        self.elenergies = i_en
        self.elenergies = Manager().convert_energy_2_internal_u(self.elenergies)

        self.Nel = count
        self.nel = self.Nel
        #
        #  The system has to be built before it can be used
        #
        self._built = False

        self._diagonalized = False
        
        self._has_nat_lifetime = [False]*self.nel
        self._nat_lifetime = numpy.zeros(self.nel)
        self._nat_linewidth = numpy.zeros(self.nel) 

        #
        #  Energy gap correlation functions
        #
        self._is_mapped_on_egcf_matrix = False
        self._has_egcf_matrix = False #
        self.egcf_matrix = None

        self.triangle = triangle(self.Nel)
        
        self._egcf_initialized = True
        self._has_egcf = self.triangle.get_list(init=False)
        self.egcf = self.triangle.get_empty_list()
        

        #
        #   System-bath interaction
        #
        self._has_system_bath_interaction = False 


        #
        #  Relaxation tensor information
        #
        self.RelaxationTensor = None
        self.RelaxationHamiltonian = None
        self._has_relaxation_tensor = False
        self._relaxation_theory = "standard_Redfield"
        self.has_Iterm = False
        

        #
        #  Electronic and other states
        #
        #self.Nel = 0
        self.Ntot = self.Nel
        

        # 
        #  Transition dipole moment 
        #
        self.dmoments = numpy.zeros((self.Nel,self.Nel,3))
        
        self.HH = None

        #
        # RWA
        #
        self.el_rwa_indices = None
        self.has_rwa = False

        #
        #  Spectroscopic information
        #
        self.WPM = None # weighted participation matrix (for spectroscopy)
        self._has_wpm = False
        

    def diagonalize(self):
        """Prepage the diagonal form of the Hamiltonian
        
        In the present case, it is already diagonal

        """
        if self._diagonalized:
            return

        self.HD = numpy.diag(self.HH)

        # diagonal transformation matrix
        self.SS = numpy.diag(numpy.ones_like(self.HD))

        self._diagonalized = True


    def get_Hamiltonian(self):
        """ Returns the Hamiltonian operator

        """
        if self._built:

            return self.HamOp

        else:

            raise Exception("System has to be built first")


    def get_TransitionDipoleMoment(self):
        """ Returns the asystem's transition dipole moment operator

        """
        if self._built:
            return self.TrDMOp # TransitionDipoleMoment(data=self.DD)
        else:
            raise Exception("Aggregate object not built")


    def set_transition_environment(self, transition, egcf):
        """Sets a correlation function for a transition on this monomer
        
        Parameters
        ----------
        transition : tuple
            A tuple describing a transition in the molecule, e.g. (0,1) is
            a transition from the ground state to the first excited state.

        egcf : CorrelationFunction
            CorrelationFunction object 
            
        (Should go to Open System)
        
        """

        if self._is_mapped_on_egcf_matrix:
            raise Exception("This system is mapped" 
                            + " on a CorrelationFunctionMatrix")

        if not (self._has_egcf[self.triangle.locate(transition[0],
                                                    transition[1])]):
                                                        
            self.egcf[self.triangle.locate(transition[0],
                                           transition[1])] = egcf
                                           
            self._has_egcf[self.triangle.locate(transition[0],
                                                transition[1])] = True
            self._has_system_bath_coupling = True                                          
                                                
        else:
            raise Exception("Correlation function already speficied" +
                            " for this monomer")


    def unset_transition_environment(self, transition):
        """Unsets correlation function from a transition on this monomer
        
        This is needed if the environment is to be replaced
        
        Parameters
        ----------
        transition : tuple
            A tuple describing a transition in the molecule, e.g. (0,1) is
            a transition from the ground state to the first excited state.
 
        (Should go to OpenSystem)
            
        """
        
        if self._is_mapped_on_egcf_matrix:
            raise Exception("This monomer is mapped" 
                            + " on a CorrelationFunctionMatrix") 
             
        if self._has_egcf[self.triangle.locate(transition[0], transition[1])]:

            self.egcf[self.triangle.locate(transition[0],
                                           transition[1])] = None
                                           
            self._has_egcf[self.triangle.locate(transition[0],
                                                transition[1])] = False
        
        self._has_system_bath_coupling = False


    def get_transition_environment(self, transition):
        """Returns energy gap correlation function of a monomer
        
        Parameters
        ----------
        transition : tuple
            A tuple describing a transition in the molecule, e.g. (0,1) is
            a transition from the ground state to the first excited state.
 
        (Should go to OpenSystem)    
        
        """

        if self._has_egcf[self.triangle.locate(transition[0] ,transition[1])]:
            return self.egcf[self.triangle.locate(transition[0],
                                                  transition[1])]
 
        if self._is_mapped_on_egcf_matrix:
            n = self.egcf_mapping[0]
            iof = self.egcf_matrix.get_index_by_where((n,n))
            
            if iof >= 0:
                return self.egcf_matrix.cfunc[iof]

        raise Exception("No environment set for the transition")   


    def set_egcf_mapping(self, transition, correlation_matrix, position):
        """ Sets a correlation function mapping for a selected transition.
        
        The monomer can either have a correlation function assigned to it,
        or it can be a part of a correlation matrix. Here the mapping to the
        correlation matrix is specified. 
        
        Parameters
        ----------
        transition : tuple
            A tuple describing a transition in the molecule, e.g. (0,1) is
            a transition from the ground state to the first excited state.
            
        correlation_matrix : CorrelationFunctionMatrix
            An instance of CorrelationFunctionMatrix
            
        position : int
            Position in the CorrelationFunctionMatrix corresponding
            to the monomer. 
            
            
        (Should go to OpenSystems)
        
        """
        
        if not (self._has_egcf[self.triangle.locate(transition[0],
                                                    transition[1])]):
                                                        
            if not (self._is_mapped_on_egcf_matrix):
                
                self.egcf_matrix = correlation_matrix
                self.egcf_transitions = []
                self.egcf_mapping = []

                
            self.egcf_transitions.append(transition)
            self.egcf_mapping.append(position)
            self._has_egcf[self.triangle.locate(transition[0],
                                                transition[1])] = True

            self._is_mapped_on_egcf_matrix = True
            self._has_system_bath_coupling = True

        else:
            
            raise Exception("Monomer has a correlation function already")            
    

    def get_RWA_suggestion(self):
        """Returns a suggestion of the RWA frequency
        
        """
        #return self.convert_energy_2_current_u(self.HH[1,1])
    
        indxs = self.get_band(band=1)
        e_av = 0.0
        for ii in indxs:
            e_av += self.HH[ii,ii]
        e_av = e_av/len(indxs)

        return self.convert_energy_2_current_u(e_av)


    def get_band(self, band=1):
        """ Returns a tuple of state indices for a given band of states
        
        """

        Nbefore = 0
        for ii in range(band):
            Nbefore += self.Nb[ii]
        Nin = self.Nb[band]
        lst = [k for k in range(Nbefore, Nbefore+Nin)]

        return tuple(lst)


    def set_dipole(self, N, M, vec=None):
        """Sets transition dipole moment for an electronic transition
        
        
        There are two ways how to use this function:
            
            1) recommended
                
                set_dipole((0,1),[1.0, 0.0, 0.0])
                
                here N represents a transition by a tuple, M is the dipole
                
            2) deprecated (used in earlier versions of quantarhei)
                
                set_dipole(0,1,[1.0, 0.0, 0.0])
                
                here the transition is characterized by two integers
                and the last arguments is the vector
                
        
        
        Examples
        --------
 
        >>> m = OpenSystem([0.0, 1.0])
        >>> m.set_dipole((0,1),[1.0, 0.0, 0.0])
        >>> m.get_dipole((0,1))
        array([ 1.,  0.,  0.])

        
        """
        if vec is None:
            n = N[0]
            m = N[1]
            vc = M
        else:
            n = N
            m = M
            vc = vec
            
        if n == m:
            raise Exception("M must not be equal to N")
        try:
            self.dmoments[n, m, :] = vc
            self.dmoments[m, n, :] = numpy.conj(vc)
        except:
            raise Exception()


    def get_dipole(self, N, M=None):
        """Returns the dipole vector for a given electronic transition
        
        There are two ways how to use this function:
            
            1) recommended
                
                get_dipole((0,1),[1.0, 0.0, 0.0])
                
                here N represents a transition by a tuple, M is the dipole
                
            2) deprecated (used in earlier versions of quantarhei)
                
                get_dipole(0,1,[1.0, 0.0, 0.0])
                
                here the transition is characterized by two integers
                and the last arguments is the vector        
                
        
        Examples
        --------
        
        
        >>> m = OpenSystem([0.0, 1.0])
        >>> m.set_dipole((0,1),[1.0, 0.0, 0.0])
        >>> m.get_dipole((0,1))
        array([ 1.,  0.,  0.])
        
        
        """
        if M is None:
            n = N[0]
            m = N[1]
        else:
            n = N
            m = M
            
        try:
            return self.dmoments[n, m, :]
        except:
            raise Exception()




    # FIXME: There are two functions with similar results 
    def map_egcf_to_states(self, mpx):
        """Returns a mapping between states and correlation functions"""

        ss = self.SS
        Ng = self.Nb[0]
        Ne1 = self.Nb[1] + Ng

        Ne1 = self.Nel

        WPM = numpy.einsum("na,nb,ni->abi",ss[Ng:Ne1,Ng:Ne1]**2,
                                            ss[Ng:Ne1,Ng:Ne1]**2,mpx)

        return WPM


    def get_lineshape_functions(self):
        """Returns lineshape functions defined for this system
        
        """
        sbi = self.get_SystemBathInteraction()
        return sbi.get_goft_storage()


    def map_lineshape_to_states(self, mpx):
        """Maps the participation matrix on g(t) functions storage
        
        This should work for a simple mapping matrix and first excited band.
        Specialized mapping such as those for 2 exciton states is implemented
        in classes that inherite from here.
        
        """

        ss = self.SS
        Ng = self.Nb[0]
        Ne1 = self.Nb[1] + Ng
        
        if self.mult > 1:
            raise Exception("Participation matrix not implemented for"+
                            "multiplicity higher than mult=1.")

        WPM = numpy.einsum("na,nb,ni->abi",ss[Ng:Ne1,Ng:Ne1]**2,
                                           ss[Ng:Ne1,Ng:Ne1]**2,mpx) 
            
        return WPM
    

    def get_weighted_participation(self):
        """Returns a participation matrix weighted by the mapping onto g(t) storage
        
        
        """
        self.diagonalize()
        
        if self._has_wpm:
            return self.WPM
        
        # mapping of the function storage
        gg = self.get_lineshape_functions()
        mpx = gg.get_mapping_matrix()

        # here we get the mapping of energy gap correlation functions to states
        self.WPM = self.map_egcf_to_states(mpx)
        self._has_wpm = True               
        
        return self.WPM


    def get_eigenstate_energies(self):
        """Returns the energies of the system's eigenstates
        
        """
        self.diagonalize()
        
        return self.HD


    def get_F4d(self, which="bbaa"):
        """Get vector of transition dipole moments for 4 wave-mixing 
        
        F4 is one of the vectors/matrices needed for orientational average
        of the third order non-linear signal.
        
        Parameters
        ----------
        
        which : string
            A string characterizing the type of dipole moment product
            
        """
        
        def DD(aa,bb):
            return numpy.dot(self.DD[aa[1],aa[0]],self.DD[bb[1],bb[0]])
        
        def _setF4(F4, x1, x2, x3, x4):
            F4[0] = DD(x4,x3)*DD(x2,x1)
            F4[1] = DD(x4,x2)*DD(x3,x1)
            F4[2] = DD(x4,x1)*DD(x3,x2) 
        
        gg = 0
        Ng = self.Nb[0]
        band1 = self.get_band(1)
        N1b = self.Nb[1] 
        if self.mult > 1:
            band2 = self.get_band(2)
            N2b = self.Nb[2]
            
        if which == "abba":
            F4 = numpy.zeros((N1b,N1b,3), dtype=REAL)
            for aa in band1:
                a = aa - Ng
                for bb in band1:
                    b = bb - Ng
                    _setF4(F4[a,b,:],(aa,gg),(bb,gg),(bb,gg),(aa,gg))
                    
            
        elif which == "baba":
            F4 = numpy.zeros((N1b,N1b,3), dtype=REAL)
            for aa in band1:
                a = aa - Ng
                for bb in band1:
                    b = bb - Ng
                    _setF4(F4[a,b],(aa,gg),(bb,gg),(aa,gg),(bb,gg))
                    
        elif which == "bbaa":
            F4 = numpy.zeros((N1b,N1b,3), dtype=REAL)
            for aa in band1:
                a = aa - Ng
                for bb in band1:
                    b = bb - Ng
                    _setF4(F4[a,b,:],(aa,gg),(aa,gg),(bb,gg),(bb,gg))
                    
        elif which == "fbfaba":
            F4 = numpy.zeros((N2b,N1b,N1b,3), dtype=REAL)
            for aa in band1:
                a = aa - Ng
                for bb in band1:
                    b = bb - 1
                    for ff in band2:
                        f = ff - Ng - N1b
                        _setF4(F4[f,a,b,:],(aa,gg),(bb,gg),(ff,aa),(ff,bb))
                    
        elif which == "fafbba":
            F4 = numpy.zeros((N2b,N1b,N1b,3), dtype=REAL)    
            for aa in band1:
                a = aa - Ng
                for bb in band1:
                    b = bb - Ng
                    for ff in band2:
                        f = ff - Ng - N1b
                        _setF4(F4[f,a,b,:], (aa,gg), (bb,gg), (ff,bb), (ff,aa)) 
                        
        elif which == "fbfbaa":
            F4 = numpy.zeros((N2b,N1b,N1b,3), dtype=REAL)    
            for aa in band1:
                a = aa - Ng
                for bb in band1:
                    b = bb - Ng
                    for ff in band2:
                        f = ff - Ng - N1b
                        _setF4(F4[f,a,b,:], (aa,gg), (aa,gg), (ff,bb), (ff,bb)) 
                    
        return F4
        

    def get_SystemBathInteraction(self):
        """Returns the aggregate SystemBathInteraction object

        """
        if self._built:
            return self.sbi
        else:
            raise Exception("The object not built.")       
    
    

    def get_RelaxationTensor(self, timeaxis,
                       relaxation_theory=None,
                       as_convolution_kernel=False,
                       time_dependent=False,
                       secular_relaxation=False,
                       relaxation_cutoff_time=None,
                       coupling_cutoff=None,
                       recalculate=True,
                       as_operators=False):
        """Returns a relaxation tensor corresponding to the system


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
        from ..qm import ModRedfieldRelaxationTensor
        from ..qm import TDModRedfieldRelaxationTensor
        from ..qm import FoersterRelaxationTensor
        from ..qm import TDFoersterRelaxationTensor
        from ..qm import NEFoersterRelaxationTensor
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
        theories["noneq_Foerster"] = ["noneq_Foerster", "neF"]
        
        #
        # Future
        #
        theories["modified_Redfield"] = ["modifield_Redfield", "mR"]
        theories["noneq_modified_Redfield"] = ["noneq_modified_Redfield",
                                               "nemR"]
        theories["generalized_Foerster"] = ["generalized_Foerster", "gF",
                                            "multichromophoric_Foerster"]
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

        elif relaxation_theory in theories["modified_Redfield"]:

            if time_dependent:

                # Time dependent standard Refield

                ham.protect_basis()
                with eigenbasis_of(ham):
                    relaxT = TDModRedfieldRelaxationTensor(ham, sbi,
                                        cutoff_time=relaxation_cutoff_time,
                                        as_operators=as_operators)
                    if secular_relaxation:
                        relaxT.secularize()
                ham.unprotect_basis()

            else:

                # Time independent standard Refield


                ham.protect_basis()

                with eigenbasis_of(ham):
                    relaxT = ModRedfieldRelaxationTensor(ham, sbi,
                                                    as_operators=as_operators)

                    if secular_relaxation:
                        relaxT.secularize()

                ham.unprotect_basis()


            self.RelaxationTensor = relaxT
            self.RelaxationHamiltonian = ham
            self._has_relaxation_tensor = True
            self._relaxation_theory = "modified_Redfield"

            return relaxT, ham


        elif relaxation_theory in theories["standard_Foerster"]:

            if time_dependent:

                # Time dependent standard Foerster
                relaxT = TDFoersterRelaxationTensor(ham, sbi)
                dat = numpy.zeros((ham.dim,ham.dim),dtype=REAL)
                for i in range(ham.dim):
                    dat[i,i] = ham._data[i,i]
                ham_0 = Hamiltonian(data=dat)
                ham_0.set_rwa(ham.rwa_indices)

            else:

                # Time independent standard Foerster

                #
                # This is done strictly in site basis
                #

                relaxT = FoersterRelaxationTensor(ham, sbi, 
                                                  pure_dephasing=True)
                dat = numpy.zeros((ham.dim,ham.dim),dtype=REAL)
                for i in range(ham.dim):
                    dat[i,i] = ham._data[i,i]
                ham_0 = Hamiltonian(data=dat)
                ham_0.set_rwa(ham.rwa_indices)

            # The Hamiltonian for propagation is the one without
            # resonance coupling
            self.RelaxationTensor = relaxT
            self.RelaxationHamiltonian = ham_0
            self._has_relaxation_tensor = True
            self._relaxation_theory = "standard_Foerster"

            return relaxT, ham_0

        elif relaxation_theory in theories["noneq_Foerster"]:

            if time_dependent:


                # Time dependent standard Foerster
                relaxT = NEFoersterRelaxationTensor(ham, sbi, 
                                            as_kernel=as_convolution_kernel)
                
                
                dat = numpy.zeros((ham.dim,ham.dim),dtype=REAL)
                
                #for i in range(ham.dim):
                #    dat[i,i] = ham._data[i,i]
                
                ham_0 = Hamiltonian(data=dat)
                
                # this type of theory has an inhomogeneous term
                self.has_Iterm = True

            else:

                # Fixme Check that it makes sense to have here
                # the time-independent theory !!!
                # Time independent standard Foerster

                #
                # This is done strictly in site basis
                #

                relaxT = FoersterRelaxationTensor(ham, sbi)
                dat = numpy.zeros((ham.dim,ham.dim),dtype=REAL)
                for i in range(ham.dim):
                    dat[i,i] = ham._data[i,i]
                ham_0 = Hamiltonian(data=dat)

            # The Hamiltonian for propagation is the one without
            # resonance coupling
            self.RelaxationTensor = relaxT
            self.RelaxationHamiltonian = ham_0
            self._has_relaxation_tensor = True
            self._relaxation_theory = "noneq_Foerster"

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

                    #print("Last line of the context", 
                    #      Manager().get_current_basis())
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
                        vKK = numpy.zeros((eKK.shape[0], ham.dim, ham.dim),
                                          dtype=REAL)

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
                       time_nonlocal=False,
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
                       as_convolution_kernel=time_nonlocal,
                       time_dependent=time_dependent,
                       secular_relaxation=secular_relaxation,
                       relaxation_cutoff_time=relaxation_cutoff_time,
                       coupling_cutoff=coupling_cutoff,
                       recalculate=recalculate,
                       as_operators=as_operators)
        

        if time_nonlocal:
            # here we create time non-local propagator
            prop = None
            raise Exception()
        else:
            # FIXME: is the eigenbases needed???
            #with eigenbasis_of(ham):
            prop = ReducedDensityMatrixPropagator(timeaxis, ham, relaxT)

        return prop


    #FIXME: There must be a general theory here
    def get_RedfieldRateMatrix(self):
        """Returns Redfield rate matrix
        
        """

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
        """Returns Förster rate matrix for the open system
        
        
        """

        from ..qm import FoersterRateMatrix
        
        if self._built:        
            ham = self.get_Hamiltonian()
            sbi = self.get_SystemBathInteraction()
        else:
            raise Exception()

        return FoersterRateMatrix(ham, sbi)

    

    def get_KTHierarchy(self, depth=2):
        """Returns the Kubo-Tanimura hierarchy of an open system
        
        """
        
        HH = self.get_Hamiltonian()
        HH.set_rwa([0,1])
        sbi = self.get_SystemBathInteraction()
        return KTHierarchy(HH, sbi, depth=depth)
    
    
    def get_KTHierarchyPropagator(self, depth=2):
        """Returns a propagator based on the Kubo-Tanimura hierarchy
        
        """
        
        kth = self.get_KTHierarchy(depth)
        ta = kth.sbi.TimeAxis
        
        return KTHierarchyPropagator(ta, kth)
    
    
    def get_excited_density_matrix(self, condition="delta", polarization=None):
        """Returns the density matrix corresponding to excitation condition
        
        
        """
        
        dip = self.get_TransitionDipoleMoment()
        if polarization is None:
            dip = dip.get_dipole_length_operator()
        else:
            # FIXME: This method is not implemented yet
            dip = dip.get_dipole_projection(polarization)
            
        with energy_units("int"):
            rho0 = self.get_thermal_ReducedDensityMatrix()
        
        
        if isinstance(condition, str):
            cond = condition
            
        else:
            cond = condition[0]

        
        if cond == "delta":
            
            rdi = numpy.dot(dip.data,numpy.dot(rho0.data,dip.data))
            rhoi = ReducedDensityMatrix(data=rdi)
            
            return rhoi

        elif cond == "pulse_spectrum":
            
            spectrum = condition[1]
            
            HH = self.get_Hamiltonian()
            
            if isinstance(spectrum, DFunction):
                
                dat = numpy.zeros((HH.dim,HH.dim), dtype=REAL)
                with eigenbasis_of(HH):
                    
                    for ii in range(HH.dim):
                        for jj in range(ii):
                            # frequency will always be >= 0.0
                            freque = HH.data[ii,ii] - HH.data[jj,jj]
                            if ((spectrum.axis.max > freque) 
                                and (spectrum.axis.min < freque)):
                                weight = numpy.sqrt(spectrum.at(freque))
                            else:
                                weight = 0.0
                            
                            dat[ii,jj] = weight*dip.data[ii,jj]
                            dat[jj,ii] = dat[ii,jj]
                            
                    dip = SelfAdjointOperator(data=dat)        
                    rdi = numpy.dot(dip.data,numpy.dot(rho0.data,dip.data))
                    rhoi = ReducedDensityMatrix(data=rdi)
                            
                return rhoi
            
            else:
                
                raise Exception("Spectrum must be specified through"+
                                " a DFunction object")

        else:
            print("Excitation condition:", condition)
            raise Exception("Excition condition not implemented.")
 
    
 
    def get_thermal_ReducedDensityMatrix(self):
        """Returns equilibrium density matrix for a give temparature
        
        """
        
        H = self.get_Hamiltonian() 
        T = self.get_temperature()
        dat = numpy.zeros(H._data.shape,dtype=COMPLEX)
        
        with eigenbasis_of(H):
            
            if numpy.abs(T) < 1.0e-10:
                dat[0,0] = 1.0
            
            else:
            
                dsum = 0.0
                
                for n in range(H._data.shape[0]):
                    dat[n,n] = numpy.exp(-H.data[n,n]/(kB_intK*T))
                    dsum += dat[n,n]

                dat *= 1.0/dsum
            
            rdm = ReducedDensityMatrix(data=dat)
                
        return rdm 
    

    def integrate_deposited_energy(self, rho_t, field, ome=None):
        """ Integrates the energy deposited into the system by light

        By default, rho_t and field contain the optical frequency. If ome and dt are defined,
        it is assumed that rho_t was calculated using a single frequency (global) RWA. In that
        case the field is assumed to be the field's envelop.

        Parameters
        ----------

        rho_t : DensityMatrixEvolution
            The dynamics of the system in terms of its time-dependent density matrix

        field: real array
            Positive frequency (e^{-i\omega t}) component of the electric field or its envelop
            if RWA is assumed.

        ome: real
            Field RWA frequency

        dt: real
            Time step on which rho_t and field are defined
            
        """

        # transition dipole moment 
        DD = self.get_TransitionDipoleMoment()


        # X component of the transition dipole moment
        # FIXME: here and in the propagator, we have to find some way how to use the whole of DD and
        #        the field polaritation
        mu = DD.data[:,:,0]

        # upper and lower triangle of the dipole moment
        #mu_u = numpy.triu(mu, k=1)
        mu_l = numpy.tril(mu, k=-1)
        
        # complex components of the field
        Ep = field
        #Em = numpy.conj(field)

        # quasi-rotating wave approximation (resonant parts only)
        #trdip_u = numpy.einsum("ij,tji->t", mu_u, rho_t.data)
        trdip_l = numpy.einsum("ij,tji->t", mu_l, rho_t.data)

        # numerical integration with the derivative of the field
        #integ_u = integral_g_dfdt(trdip_u, Em)
        integ_l = integral_g_dfdt(trdip_l, Ep)
        
        if ome is not None:
            dt = rho_t.TimeAxis.step
            ome_in = Manager().convert_energy_2_internal_u(ome)
            Epom = 1j*ome_in*Ep
            
            integ_rwa = integral_g_f(trdip_l, Epom, dt)
        else:
            integ_rwa = 0.0

        val = -2.0*numpy.real(integ_l) + 2.0*numpy.real(integ_rwa) # + integ_u)
        
        return Manager().convert_energy_2_current_u(val)

    def set_electronic_natural_lifetime(self, N, epsilon_r=1.0):
        
        rate = 0.0
        eps = eps0_int*epsilon_r
        
        try: 
            # loop over all states with lower energy
            for n in range(N):
                dip2 = numpy.dot(self.dmoments[n,N,:],
                                 self.dmoments[n,N,:])
                ome = self.elenergies[N] - self.elenergies[n]                        
        
                rate += dip2*(ome**3)/(3.0*numpy.pi*eps*(c_int**3))
                
        except:
            raise Exception("Calculation of rate failed")
            
        if rate > 0.0:
            
            lftm = 1.0/rate
            self._has_nat_lifetime[N] = True
            self._nat_lifetime[N] = lftm
            # FIXME: calculate the linewidth
            self._nat_linewidth[N] = 1.0/lftm
            self._saved_epsilon_r = epsilon_r
            
        else:
            #raise Exception("Cannot calculate natural lifetime")
            self._has_nat_lifetime[N] = True
            self._nat_lifetime[N] = numpy.inf   
            
            
    def get_electronic_natural_lifetime(self, N, epsilon_r=1.0):
        """Returns natural lifetime of a given electronic state
        
        
        Examplex
        --------
        
        Molecule where transition dipole was not set does not calculate lifetime
        
        >>> from quantarhei import Molecule
        >>> m = Molecule([0.0, 1.0])
        >>> m.get_electronic_natural_lifetime(1)
        Traceback (most recent call last):
        ...
        AttributeError: 'Molecule' object has no attribute '_saved_epsilon_r'
        
        
        Setting transition dipole moment allows calculation of the lifetime
        
        >>> m = Molecule([0.0, 1.0])
        >>> m.set_dipole((0,1), [5.0, 0.0, 0.0])
        >>> print(m.get_electronic_natural_lifetime(1))
        852431568.119

        The value is recalculated only when relative permitivity is changed

        >>> m = Molecule([0.0, 1.0])
        >>> m.set_dipole((0,1), [5.0, 0.0, 0.0])
        >>> tl1 = m.get_electronic_natural_lifetime(1)
        >>> m.set_dipole((0,1), [2.0, 0.0, 0.0]) 
        >>> tl2 = m.get_electronic_natural_lifetime(1)
        >>> tl3 = m.get_electronic_natural_lifetime(1, epsilon_r=2.0)
        >>> (tl2 == tl1) and (tl2 != tl3)
        True
        
        
        """     
        
        
        if not self._has_nat_lifetime[N]:
            self.set_electronic_natural_lifetime(N,epsilon_r=epsilon_r)
            

        if self._saved_epsilon_r != epsilon_r:
            self.set_electronic_natural_lifetime(N,epsilon_r=epsilon_r)

        return self._nat_lifetime[N]

#
# Auxiliary routines
#
def integral_g_dfdt(g,f):
    """
    Numerically compute ∫ g(t) * df/dt(t) dt using central differences.
    Assumes df/dt = 0 at boundaries.
    
    Parameters:
        g (ndarray): Function values to multiply with df/dt.
        f (ndarray): Function values at discrete time points.
        # dt (float): Time step between points. (Time step cancels out)
        
    Returns:
        float: Approximation of the integral.
    """
    dfdt = numpy.empty_like(f)
    dfdt[0] = (f[1] - f[0])
    dfdt[-1] = (f[-1] - f[-2])
    dfdt[1:-1] = (f[2:] - f[:-2])/2.0
    
    return numpy.sum(g*dfdt)


def integral_g_f(g, f, dt):
    """
    Numerically compute ∫ g(t) * f(t) dt using the trapezoidal rule.
    Assumes uniform time step (dt cancels out, as in the df/dt version).

    Parameters:
        g (ndarray): Function values.
        f (ndarray): Function values.

    Returns:
        float: Approximate integral (up to a constant dt).
    """
    integrand = g*f
    result = 0.5*(integrand[0] + integrand[-1]) + numpy.sum(integrand[1:-1])
    return result*dt 
