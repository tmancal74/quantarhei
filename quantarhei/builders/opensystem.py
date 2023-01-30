# -*- coding: utf-8 -*-
import numpy


from ..qm.hilbertspace.hamiltonian import Hamiltonian


class OpenSystem:
    """The class representing a general open quantum system
    
    It provides routines to store and extract the most interesting
    characteristics of an open qwuantum system, such as the Hamiltonian,
    SystemBathInteraction object, relaxation tensor and the reduced density
    matrix evolution.
    
    
    """
    


    def __init__(self):
        
        self._built = False
        
        self.RelaxationTensor = None
        self.RelaxationHamiltonian = None
        self._has_relaxation_tensor = False
        self._relaxation_theory = "standard_Redfield"
        self.has_Iterm = False
        
        self.Nel = 0
        self.Ntot = 0
        


    def get_Hamiltonian(self):
        """Returns the system Hamiltonian 
        
        
        This is a method to be implemented by the particular system
        
        """
        
        return None


    def get_SystemBathInteraction(self):
        """Returns the system-bath interaction definition
        
        
        This is a method to be implemented by the particular system
        
        """
        
        return None        
    
    

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

        elif relaxation_theory in theories["noneq_Foerster"]:

            if time_dependent:

                # Time dependent standard Foerster
                relaxT = NEFoersterRelaxationTensor(ham, sbi, 
                                            as_kernel=as_convolution_kernel)
                dat = numpy.zeros((ham.dim,ham.dim),dtype=numpy.float64)
                
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
                dat = numpy.zeros((ham.dim,ham.dim),dtype=numpy.float64)
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

    