# -*- coding: utf-8 -*-
"""Multi-level molecule (monomer)

    The molecule is defined by the vector of energies of its states
    and by the transition dipole moments between allowed transitions.

    >>> m = Molecule([0.0, 1.0], "Alphonso")
    >>> print("Hi, my name is", m.get_name()+"!")
    Hi, my name is Alphonso!


    Information about the molecule can be obtained simply by printing it
    
    >>> print(m)
    <BLANKLINE>
    quantarhei.Molecule object
    ==========================
       name = Alphonso  
       position = None
       number of electronic states = 2
       # State properties
       State nr: 0 (ground state)
          electronic energy = 0.0 1/fs
          number of vibrational modes = 0
    <BLANKLINE>
       State nr: 1
          electronic energy = 1.0 1/fs
          transition 0 -> 1 
          transition dipole moment = [0.0, 0.0, 0.0]
          number of vibrational modes = 0
    <BLANKLINE>

    
    >>> import quantarhei as qr
    >>> mol = qr.TestMolecule("two-levels-1-mode")
    >>> print(mol)
    <BLANKLINE>
    quantarhei.Molecule object
    ==========================
       name = two-levels-1-mode  
       position = None
       number of electronic states = 2
       # State properties
       State nr: 0 (ground state)
          electronic energy = 0.0 1/fs
          number of vibrational modes = 1
          # Mode properties
          mode no. = 0 
             frequency = 1.0 1/fs
             shift = 0.0
             nmax = 2
       State nr: 1
          electronic energy = 1.0 1/fs
          transition 0 -> 1 
          transition dipole moment = [0.0, 0.0, 0.0]
          number of vibrational modes = 1
          # Mode properties
          mode no. = 0 
             frequency = 1.0 1/fs
             shift = 0.0
             nmax = 2
    

    Class Details
    -------------
    
    

"""


import numpy

from ..utils import array_property
from ..utils import Integer
#from ..utils import check_numpy_array

from ..core.managers import UnitsManaged, Manager
from ..core.managers import eigenbasis_of

from . import Mode

from ..core.triangle import triangle
from ..core.unique import unique_list
from ..core.unique import unique_array

from ..core.units import kB_intK, eps0_int, c_int

from ..qm import Hamiltonian
from ..qm import TransitionDipoleMoment
from ..qm.oscillators.ho import operator_factory

from ..qm import SystemBathInteraction
from ..qm.corfunctions.cfmatrix import CorrelationFunctionMatrix

from ..qm import ReducedDensityMatrix

from ..core.saveable import Saveable

import quantarhei as qr

 
    
class Molecule(UnitsManaged, Saveable):
    """Multi-level molecule (monomer)


    Parameters
    ----------
    
    name : str
        Monomer descriptor; a string identifying the monomer

    elenergies : list of real numbers
        List of electronic energies, one per state. It includes ground state
        energy. It wise to chose the ground state energy as zero. 

    
    """
    
    # position of the monomer 
    position = array_property('position',shape=(3,))
    # energies of electronic states
    elenergies = array_property('elenergies')
    # transition dipole moments
    dmoments = array_property('dmoments')    
    
    # number of electronic statesarray_property
    nel      = Integer('nel')
    # number of allowed transitions
    nal      = Integer('nal')
    # number of vibrational modes
    nmod     = Integer('nmod')
    
    
    def __init__(self,  elenergies=[0.0,1.0], name=None): #,dmoments):
    
        #self.manager = Manager()
        
        # monomer name
        if name is None:
            # FIXME: generate unique name
            self.name = "xxx"
        else:
            self.name = name  #
            
        #
        # set energies
        # 

        # convert to internal_units
        #self.elenergies = self.manager.convert_energy_2_internal_u(elenergies)
        self.elenergies = elenergies
        self.elenergies = self.convert_energy_2_internal_u(self.elenergies) #
        self.nel = len(elenergies)    #
        
        
        
        # FIXME: check the order of energies (increasing order has 
        # to be enforced) no vibrational modes is a defauls
        self.nmod = 0
        self.modes = []  
        
        # allowed transitions are now only between the ground state and
        # the rest of the excited states are dark
        self.allowed_transitions = []
         
        # transition dipole moments
        self.dmoments = numpy.zeros((self.nel,self.nel,3)) 
        
        # matrix of the transition widths
        self.widths = None
        
        # matrix of dephasing rates
        self.dephs = None
        
        self._has_nat_lifetime = [False]*self.nel
        self._nat_lifetime = numpy.zeros(self.nel)
        self._nat_linewidth = numpy.zeros(self.nel)        
        
        # FIXME: some properties using triangle should be initialized 
        # only when used. triangle should stay here, other things
        # should be moved to the setters (see set_adiabatic_coupling)
        self.triangle = triangle(self.nel)
        
        self._egcf_initialized = True
        self._has_egcf = self.triangle.get_list(init=False)
        self.egcf = self.triangle.get_empty_list()
        self._has_egcf_matrix = False
        
        self._adiabatic_initialized = False
        
        # 
        # We can have an environment at each mode and each electronic state
        #
        self._mode_env_initialized = False
        
        
        self._is_mapped_on_egcf_matrix = False
        
        self._has_system_bath_coupling = False
        
        # data attribute can hold PDB coordinates or something else
        self.model = None
        self.data = None
        self._data_type = None

        
    def get_name(self):
        """Returns the name of the molecule
        
        Examples
        --------
        
        >>> m = Molecule([0.0, 1.0], name="Jane")
        >>> m.get_name()
        'Jane'
        
        """
        return self.name

        
    def set_name(self, name):
        """Sets the name of the Molecule object
        
        Examples
        --------
        
        >>> m = Molecule([0.0, 1.0])
        >>> m.set_name("Jane")
        >>> print(m.get_name())
        Jane
        
        """
        self.name = name
        
        
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
        """
        
        if not (self._has_egcf[self.triangle.locate(transition[0],
                                                    transition[1])]):
                                                        
            if not (self._is_mapped_on_egcf_matrix):
                
                self.egcf_matrix = correlation_matrix
                self.egcf_transitions = []
                self.egcf_transitions.append(transition)
                self.egcf_mapping = []
                self.egcf_mapping.append(position)
                self._is_mapped_on_egcf_matrix = True
                self._has_system_bath_coupling = True
                
            else:
                
                if self.egcf_matrix is correlation_matrix:
                    
                    if self.egcf_transitions.count(transition) > 0:
                        raise Exception("Correlation function already \
                         assigned to the transition")
                    else:
                        self.egcf_transitions.append(transition)
                        self.egcf_mapping.append(position)                        

        else:
            
            raise Exception("Monomer has a correlation function already")            
        

    def set_transition_environment(self, transition, egcf):
        """Sets a correlation function for a transition on this monomer
        
        Parameters
        ----------
        transition : tuple
            A tuple describing a transition in the molecule, e.g. (0,1) is
            a transition from the ground state to the first excited state.

        egcf : CorrelationFunction
            CorrelationFunction object 
            
            
        """
        if self._is_mapped_on_egcf_matrix:
            raise Exception("This monomer is mapped \
            on a CorrelationFunctionMatrix")

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
        
        if self._is_mapped_on_egcf_matrix:
            raise Exception("This monomer is mapped \
            on a CorrelationFunctionMatrix")  
             
        if self._has_egcf[self.triangle.locate(transition[0], transition[1])]:

            self.egcf[self.triangle.locate(transition[0],
                                           transition[1])] = None
                                           
            self._has_egcf[self.triangle.locate(transition[0],
                                                transition[1])] = False
        
    #@deprecated
    def set_egcf(self, transition, egcf):
        self.set_transition_environment(transition, egcf)
        
            
    def get_transition_environment(self, transition):
        """Returns energy gap correlation function of a monomer
        
        
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

   
    #@deprecated
    def get_egcf(self, transition): 
        
        if self._is_mapped_on_egcf_matrix:
            n = self.egcf_mapping[0]
            return self.egcf_matrix.get_coft(n,n)
            
        return self.get_transition_environment(transition).data         
            
    
    def add_Mode(self, mod):
        """Adds a vibrational mode to the monomer
        
        Parameters
        ----------
        
        mod : quantarhei.Mode
            Intramolecular vibrational mode
            
            
        Examples
        --------
        
        >>> mode = qr.Mode()
        >>> mol = Molecule([0.0, 2.0])
        >>> print(mol.get_number_of_modes())
        0
        >>> mol.add_Mode(mode)
        >>> print(mol.get_number_of_modes())
        1
        
        """
        if isinstance(mod, Mode):
            mod.set_Molecule(self)
            self.modes.append(mod)
            self.nmod += 1
        else:
            raise TypeError()
            
    #
    # Some getters and setters 
    #
    def get_Mode(self, N):
        """Returns the Nth mode of the Molecule object
        
        
        Parameters
        ----------
        
        N : int
            Index of the mode to be returned
            
            
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mod = mol.get_Mode(1)
        Traceback (most recent call last):
            ...
        IndexError: list index out of range
        
        >>> mod = mol.get_Mode(0)
        >>> mod.get_energy(1)
        1.0
        
        """
        
        return self.modes[N]

   
    def get_number_of_modes(self):
        """Retruns the number of modes in this molecule
        
        
        Examples
        --------
        

        >>> m = Molecule([0.0, 1.0])
        >>> m.get_number_of_modes()
        0

        >>> import quantarhei as qr
        >>> m = qr.TestMolecule("two-levels-1-mode")
        >>> m.get_number_of_modes()
        1
        
        """
        return len(self.modes)

        
    def get_dipole(self, N, M):
        try:
            return self.dmoments[N, M, :]
        except:
            raise Exception()

            
    def set_dipole(self, N, M, vec):
        if N == M:
            raise Exception("M must not be equal to N")
        try:
            self.dmoments[N, M, :] = vec
            self.dmoments[M, N, :] = numpy.conj(vec)
        except:
            raise Exception()  

            
    def set_transition_width(self, transition, width):
        """Sets the width of a given transition
        
        
        Parameters
        ----------
        
        transition : {tuple, list}
            Quantum numbers of the states between which the transition occurs
            
        width : float
            The full width at half maximum (FWHM) of a Gaussian lineshape,
            or half width at half maximum (HWHM) of a Lorentzian lineshape
            
            
        """
        cwidth = qr.Manager().convert_energy_2_internal_u(width)
        if self.widths is None:
            N = self.elenergies.shape[0]
            self.widths = numpy.zeros((N, N), dtype=qr.REAL)
        self.widths[transition[0], transition[1]] = cwidth
        self.widths[transition[1], transition[0]] = cwidth


    def get_transition_width(self, transition):
        """Returns the transition width
        
        
        Returns the full width at half maximum (FWHM) of a Gaussian lineshape,
        or half width at half maximum (HWHM) of a Lorentzian lineshape
        
        Parameters
        ----------
        
        transition : {tuple, list}
            Quantum numbers of the states between which the transition occurs
            
        """
        
        if self.widths is None:
            return 0.0
        
        return self.widths[transition[0], transition[1]]

    
    def set_transition_dephasing(self, transition, deph):
        """Sets the dephasing rate of a given transition
        
        
        Parameters
        ----------
        
        transition : {tuple, list}
            Quantum numbers of the states between which the transition occurs
            
        deph : float
            Dephasing rate of the transition
            
            
        """        
        if self.dephs is None:
            N = self.elenergies.shape[0]
            self.dephs = numpy.zeros((N, N), dtype=qr.REAL)
        self.dephs[transition[0], transition[1]] = deph
        self.dephs[transition[1], transition[0]] = deph


    def get_transition_dephasing(self, transition):
        """Returns the dephasing rate of a given transition
        
        
        Parameters
        ----------
        
        transition : {tuple, list}
            Quantum numbers of the states between which the transition occurs
            
        """
        
        if self.dephs is None:
            return 0.0
        
        return self.dephs[transition[0], transition[1]]


    def get_energy(self, N):
        """Returns energy of the Nth state of the molecule
        
        
        Parameters
        ----------
        
        N : int
            Index of the state
            
            
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mol.get_energy(1)
        1.0
        
        This methods reacts to the `energy_units` context manager
        
        >>> with qr.energy_units("1/cm"):
        ...     print("{0:.8f}".format(mol.get_energy(1)))
        5308.83745888
        
        """
        try:
            return self.convert_energy_2_current_u(self.elenergies[N])
        except:
            raise Exception()

      
    def set_energy(self, N, en):
        """Sets the energy of the Nth state of the molecule
        
        
        Parameters
        ----------
        
        N : int
            Index of the state
            
        en : float
            Energy to be assigned to the Nth state.
            
            
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mol.set_energy(1, 1.5)
        >>> mol.get_energy(1)
        1.5

        This method reacts to the `energy_units` context manager

        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> with qr.energy_units("1/cm"):
        ...     mol.set_energy(1, 5308.8374588761453)
        >>> mol.get_energy(1)
        1.0

        """
        self.elenergies[N] = self.convert_energy_2_internal_u(en)
        
            
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
        
        """            
        if not self._has_nat_lifetime[N]:
            self.set_electronic_natural_lifetime(N,epsilon_r=epsilon_r)
            

        if self._saved_epsilon_r != epsilon_r:
            self.set_electronic_natural_lifetime(N,epsilon_r=epsilon_r)

        return self._nat_lifetime[N]
    
        
    def get_thermal_ReducedDensityMatrix(self):
        
        H = self.get_Hamiltonian() 
        T = self.get_temperature()
        dat = numpy.zeros(H._data.shape,dtype=numpy.complex)
        
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
        
        
    def get_temperature(self):
        """Returns temperature of the molecule
        
        Checks if the setting of environments is consistent and than
        takes the temperature from one of the energy gaps. If no 
        environment (correlation function) is assigned to this 
        molecule, we assume zero temperature.
        
        
        """
        if self.check_temperature_consistent():
        
            try:
                egcf =  self.get_transition_environment([0,1])
            except:
                egcf = None
        
            if egcf is None:
                return 0.0
            else:
                return egcf.get_temperature()
                
        else:
            
            raise Exception("Molecular environment has"+
            "an inconsisten temperature")
        
        
    def check_temperature_consistent(self):
        """Checks that the temperature is the same for all components
        
        """
        #FIXME: implement the check
        return True

     
    def set_adiabatic_coupling(self,state1,state2,coupl):
        """Sets adiabatic coupling between two states
        
        
        """
        if not self._adiabatic_initialized:
            self._has_adiabatic = self.triangle.get_list(init=False)
            self.adiabatic_coupling = self.triangle.get_empty_list()
            self._adiabatic_initialized = True
            
        cp = self.convert_energy_2_internal_u(coupl)
        
        self.adiabatic_coupling[self.triangle.locate(state1,state2)] = cp
        self._has_adiabatic[self.triangle.locate(state1,state2)] = True


    def get_adiabatic_coupling(self,state1,state2):
        """Returns adiabatic coupling between two states
        
        
        """
        return self.adiabatic_coupling[self.triangle.locate(state1,state2)]
        
         
    def get_electronic_natural_linewidth(self,N):
        """Returns natural linewidth of a given electronic state
        
        """
            
        if not self._has_nat_lifetime[N]:
            self.get_electronic_natural_lifetime(N)
                
        return self._nat_lifetime[N]
           
           
    def get_Hamiltonian(self):
        """Returns the Hamiltonian of the Molecule object
        
        
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> H = mol.get_Hamiltonian()
        >>> print(H.dim)
        4
        
        >>> print(H.data)
        [[ 0.  0.  0.  0.]
         [ 0.  1.  0.  0.]
         [ 0.  0.  1.  0.]
         [ 0.  0.  0.  2.]]
        
        """
        # list of vibrational Hamiltonians
        lham = [None]*self.nel  
        # list of Hamiltonian dimensions
        ldim = [None]*self.nel
    
        # loop over electronic states
        for i in range(self.nel):
        
            if self.nmod > 0:
                # loop over modes
                for j in range(self.nmod):
                    # FIXME: enable more than one mode
                    if j > 0: # limits the number of modes to 1
                        raise Exception("Not yet implemented") 
                    
                    # number of vibrational states in this electronic state
                    Nvib = self.modes[j].get_nmax(i)
                
                    # shift of the PES
                    dd = self.modes[j].get_shift(i)
                    en = self.modes[j].get_energy(i)
 
                    # create the Hamiltonian
                    of = operator_factory(N=Nvib)
                    aa = of.anihilation_operator()
                    ad = of.creation_operator()
                    ones = of.unity_operator()    
                
                    hh = en*(numpy.dot(ad,aa) - (dd/numpy.sqrt(2.0))*(ad+aa)
                    + dd*dd*ones/2.0)
                
                lham[i] = hh + self.elenergies[i]*ones
                ldim[i] = hh.shape[0]
                
            else:
                hh = numpy.zeros((1,1),dtype=numpy.float)
                hh[0,0] = self.elenergies[i]
                lham[i] = hh
                ldim[i] = 1
            
        # dimension of the complete Hamiltonian
        totdim = numpy.sum(ldim)
        
        # this will be the Hamiltonian data
        ham = numpy.zeros((totdim,totdim),dtype=numpy.float)
        
        #
        # electronically diagonal part
        #
        lbound = 0
        ub = numpy.zeros(self.nel)
        lb = numpy.zeros(self.nel)
        # loop over electronic states
        for i in range(self.nel):
            ubound = lbound + ldim[i]
            ham[lbound:ubound,lbound:ubound] = lham[i]
            ub[i] = ubound
            lb[i] = lbound
            lbound = ubound
            
        lb,ub = self._sub_matrix_bounds(ldim)    
        for i in range(self.nel):
            ham[lb[i]:ub[i],lb[i]:ub[i]] = lham[i]    
        
        #
        # adiabatic coupling
        #
        if self._adiabatic_initialized:
            for i in range(self.nel):
                for j in range(i+1,self.nel):
                    if self._has_adiabatic[self.triangle.locate(i,j)]:
                        J = self.get_adiabatic_coupling(i,j)
                        hj = numpy.zeros((ldim[i],ldim[j]),dtype=numpy.float)
                        # FIXME: this works only if the frequencies
                        # of the oscillators are the same
                        for k in range(min([ldim[i],ldim[j]])):
                            hj[k,k] = J
                        ham[lb[i]:ub[i],lb[j]:ub[j]] = hj
                        ham[lb[j]:ub[j],lb[i]:ub[i]] = hj.T
                        
            
        return Hamiltonian(data=ham)
        
    def _sub_matrix_bounds(self,ldim):
        lbound = 0
        ub = numpy.zeros(self.nel,dtype=numpy.int)
        lb = numpy.zeros(self.nel,dtype=numpy.int)
        # loop over electronic states
        for i in range(self.nel):
            ubound = lbound + ldim[i]
            #ham[lbound:ubound,lbound:ubound] = lham[i]
            ub[i] = ubound  
            lb[i] = lbound
            lbound = ubound   
        return lb,ub
        
    def _ham_dimension(self):

        # list of Hamiltonian dimensions
        ldim = [None]*self.nel

        # loop over electronic states
        for i in range(self.nel):
        
            Nvib = 1
            
            # loop over modes
            for j in range(self.nmod):
                # FIXME: enable more than one mode
                if j > 0: # limits the number of modes to 1
                    raise Exception("Not yet implemented") 
                    
                # number of vibrational states in this electronic state
                Nvib = self.modes[j].get_nmax(i)
                
            ldim[i] = Nvib
            
        # dimension of the complete Hamiltonian
        totdim = numpy.sum(ldim)

        return totdim, ldim        
               
        
    def get_TransitionDipoleMoment(self):


        totdim,ldim = self._ham_dimension()

        # This will be the operator data
        dip = numpy.zeros((totdim,totdim,3),dtype=numpy.float) 
        
        lb,ub = self._sub_matrix_bounds(ldim)
        
        # FIXME: sofar only the lowest state of all is the start of 
        # optical transitions
        for k in range(1,self.nel):
            # FIXME: all just for one mode
            
            # get dipome moment vector
            dp = self.dmoments[0,k,:]
 
            dd = numpy.zeros((ldim[0],ldim[k]),dtype=numpy.float)        
        
            # now we run through the vibrational states of the groundstate
            # and of the excited state. Vibrational states we use are
            # unshifted and only n->n transitions are allowed. We therefore
            # run upto which ever index is shorter
            Nvib = min([ldim[0],ldim[k]]) 
            for l in range(3):
                for a in range(Nvib):
                    dd[a,a] = dp[l]
                
                dip[lb[0]:ub[0],lb[k]:ub[k],l] = dd
                dip[lb[k]:ub[k],lb[0]:ub[0],l] = dd.T
            
        return TransitionDipoleMoment(data=dip) 
    
    
    def get_SystemBathInteraction(self,timeAxis):
        """Returns a SystemBathInteraction object of the molecule 


        """ 
        nob = 0
        cf = unique_list()
        
        #
        # Look for pure dephasing environmental interactions on transitions
        #
        d = {}
        where = {}
        for i in range(self.nel):
            if i > 0:
                break # transitions not from ground state are not counted
            for j in range(i+1,self.nel):
                eg = self.egcf[self.triangle.locate(i,j)]
                if eg is not None:
                    # we save where the correlation function comes from
                    # we will have a list of excited states
                    if eg in where.keys():
                        ll = where[eg]
                        ll.append((nob,nob))  # this is a diagonal element
                                              # of the correlation function
                                              # matrix. nob is counting
                                              # the baths
                    else:
                        where[eg] = [(nob,nob)]
                    # for each bath, save the state of the traprint(mols_c2)
                    # nsition g -> j
                    d[nob] = j
                        
                    nob += 1
                    # save eg to unique list
                    cf.append(eg)
            
        # number of transition bath
        ntr = nob

        
        #
        # Look for linear environmental interactions with vibrational modes
        #           
        for i in range(self.nmod):
            for j in range(self.nel):
                eg = self.get_mode_environment(i,j)
                if eg is not None:
                    # as above
                    if eg in where.keys():
                        ll = where[eg]
                        ll.append((nob,nob))  
                    else:
                        where[eg] = [(nob,nob)]
                    # for each bath, save the combination of mod and elstate
                    d[nob] = (i,j)
                   
                    nob +=1
                    cf.append(eg)
        
        # number of mode environments
        nmd = nob - ntr

        # number of different instances of correlation functions
        nof = cf.get_number_of_unique_elements()
                
        # number of different baths nob = number of transition environments +
        # number of mode environments
        cfm = CorrelationFunctionMatrix(timeAxis,nob,nof)
        for i in range(nof):
            el = cf.get_element(i)
            wr = where[el]
#            print(wr)
#            for w in wr:
#                char = d[w[0]]
#                print("fce ",el, " is a bath no.", w[0]," of ",char)
            cfm.set_correlation_function(el,wr,i+1)

        #
        # System operators corresponding to the correlation functions
        # 

        sys_operators = []
        totdim,ldim = self._ham_dimension()
        # 
        # First, transition fluctuations.
        # We need to find projector on a given excited electronic state
        for n in range(ntr):
            KK = numpy.zeros((totdim,totdim),dtype=numpy.float)

            state = d[n]

            states_before = 0
            for k in range(state):
                states_before += ldim[k]
            states_inc = states_before +ldim[state]
            # fill 1 on diagonal corresponding to an electronic state "state"
            KK[states_before:states_inc,
               states_before:states_inc] = numpy.diag(
                               numpy.ones(ldim[state],dtype=numpy.float))
            
            sys_operators.append(KK)
            
        #
        # Linear coupling with oscillators corresponding to a given
        # electronic state 
        # FIXME: works only for one oscillator
        for n in range(ntr,nob):
            KK = numpy.zeros((totdim,totdim),dtype=numpy.float)
            
            state = d[n][1]
            mod   = d[n][0]
            
            # prepare coordinate q = (a^+ + a)/sqrt(2)
            shift = self.modes[mod].get_shift(state)
            of = operator_factory(N=ldim[state])
            aa = of.anihilation_operator()
            ad = of.creation_operator()
            ones = of.unity_operator()
            
            q = (aa + ad)/numpy.sqrt(2.0) - shift*ones
            
            states_before = 0
            for k in range(state):
                states_before += ldim[k]
            states_inc = states_before +ldim[state]
            # fill with shifted coordinate
            KK[states_before:states_inc,
               states_before:states_inc] = q
                           
            sys_operators.append(KK)
        
        return SystemBathInteraction(sys_operators,cfm)
        
    
    def set_mode_environment(self, mode=0, elstate=0, corfunc=None):
        """Sets mode environment 
        
        
        Sets the environment (bath correlation function) interacting with a
        a given mode in a given electronic state.
        
        Parameters
        ----------
        mode : int
            index of the mode
            
        elstate : int
            index of the electronic state
            
        corfunc: quantarhei.CorrelationFunction
            CorrelationFunction object
            
        
        """
        
        if corfunc is None:
            return
            
        if not self._mode_env_initialized:
            self._has_mode_env = numpy.zeros((self.nmod,self.nel),
                                             dtype=numpy.bool)
            self._mode_env = unique_array(self.nmod,self.nel)
            self._mode_env_initialized = True
            
        self._mode_env.set_element(mode,elstate,corfunc)
        self._has_mode_env[mode,elstate] = True
        
        
    def get_mode_environment(self,mode,elstate):
        """Returns mode environment 
        
        
        Returns the environment (bath correlation function) interacting with a
        a given mode in a given electronic state.
        
        Parameters
        ----------
        mode : int
            index of the mode
            
        elstate : int
            index of the electronic state
                      
        
        """        
        if self._mode_env_initialized and self._has_mode_env[mode,elstate]:
            return self._mode_env.get_element(mode,elstate)
        else:
            return None



    def __str__(self):
        """String representation of the Molecule object
        
        """
        
        out  = "\nquantarhei.Molecule object" 
        out += "\n=========================="
        out += "\n   name = %s  \n" % self.name
        try:
            out += "   position = [%r, %r, %r] \n" % (self.position[0],
                                                  self.position[1],
                                                  self.position[2])
        except:
            out += "   position = None\n"
            

            
        out += "   number of electronic states = %i" % self.nel
        out += "\n   # State properties"
        #out += "\n   -----------------"

        eunits = self.unit_repr(utype="energy")

        for n in range(self.nel):
            if n == 0:
                out += "\n   State nr: %i (ground state)" % n
            else:
                out += "\n   State nr: %i" % n
            # state energy
            ene = self.convert_energy_2_current_u(self.elenergies[n])
            out += "\n      electronic energy = %r %s" % (ene,eunits)
            # transition dipole moments
            for j in range(n):
                out += "\n      transition %i -> %i " % (j, n)
                out += "\n      transition dipole moment = [%r, %r, %r]" % (
               self.dmoments[n,j][0], 
               self.dmoments[n,j][1], self.dmoments[n,j][2])
            out += "\n      number of vibrational modes = %i" % self.nmod
            out += "\n"
            if self.nmod > 0:
                out += "      # Mode properties"
                #out += "\n      ----------------"
                for m1 in range(self.nmod):
                    out += "\n      mode no. = %i " % m1
                    out += ("\n         frequency = %r %s" %
                         (self.modes[m1].get_energy(n,no_conversion=False),
                          self.unit_repr(utype="energy")))
                    out += ("\n         shift = %r" %
                           self.modes[m1].get_shift(n))
                    out += ("\n         nmax = %i" %
                           self.modes[m1].get_nmax(n))
            
        return out

      
    def liouville_pathways_1(self, eUt=None, ham=None, dtol=0.01, ptol=1.0e-3,
                             etol=1.0e-6, verbose=0, lab=None):
        """ Generator of the first order Liouville pathways 
        
        
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
        dip_tol = numpy.sqrt(self.D2_max)*dtol
        evf_tol = etol
                        
        if eUt is None:
            
            # secular absorption spectrum calculation
            eUt2_dat = None
            sec = True
        
        else:
            
            raise Exception("Not implemented yet")

        lst = []        

        if sec:
            generate_1orderP_sec(self, lst,
                                 pop_tol, dip_tol, verbose)
        else:
            raise Exception("Not implemented yet")                                
        
        if lab is not None:
            for l in lst:
                l.orientational_averaging(lab)
         
        return lst   

def generate_1orderP_sec(self, lst,
                         pop_tol, dip_tol, verbose):
    
    ngs = (0) # self.get_electronic_groundstate()
    nes = (1) #self.get_excitonic_band(band=1)
    
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
        if True: #self.rho0[i1g,i1g] > pop_tol:
    
            D2 = numpy.dot(self.dmoments[0,1,:], self.dmoments[0,1,:])
            for i2e in nes:
                
                if D2 > dip_tol: #self.D2[i2e,i1g] > dip_tol:
                
                                
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

                        lp = \
                        diag.liouville_pathway("NR",
                                               i1g,
                                    aggregate=self,
                                    order=1,pname="P1",
                                    popt_band=1,
                                    relax_order=1)
                        
                        # first transition lineshape
                        width1 = \
                            self.get_transition_width((i2e, i1g))
                        deph1 = \
                            self.get_transition_dephasing((i2e, 
                                               i1g))
                        
                        #      |g_i1> <g_i1|                                                           
                        lp.add_transition((i2e,i1g),+1,
                                          interval=1, 
                                          width=width1, 
                                          deph=deph1)
                        #      |e_i2> <g_i1|
                        lp.add_transition((i1g,i2e),+1,
                                          interval=1, 
                                          width=width1, 
                                          deph=deph1)
                        #      |g_i1> <g_i1|

                    except:
                        
                        break
                    
                    lp.build()
                    lst.append(lp)
                    k += 1

def PiMolecule(Molecule):
    
    def __init__(self, name=None, elenergies=[0.0,1.0], data=None):
        super().__init__(name=None, elenergies=[0.0,1.0], data=None)
        self.data = data
        