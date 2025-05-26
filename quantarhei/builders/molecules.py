# -*- coding: utf-8 -*-
"""Multi-level molecule (monomer)

    The molecule is defined by the vector of energies of its states
    and by the transition dipole moments between allowed transitions.

    >>> m = Molecule([0.0, 1.0])
    >>> print(m.Nel)
    2


    Information about the molecule can be obtained simply by printing it
    
    >>> print(m)
    <BLANKLINE>
    quantarhei.Molecule object
    ==========================
       name =   
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

from ..core.managers import UnitsManaged, Manager
from ..core.managers import eigenbasis_of
from ..core.managers import energy_units

from . import Mode

from ..core.triangle import triangle
from ..core.unique import unique_list
from ..core.unique import unique_array

from ..core.units import kB_intK, eps0_int, c_int
from ..core.units import conversion_facs_length, conversion_facs_edipole

from ..qm import Hamiltonian
from ..qm import TransitionDipoleMoment

from ..qm.oscillators.ho import operator_factory

from ..qm import SystemBathInteraction
from ..qm.corfunctions.cfmatrix import CorrelationFunctionMatrix

from ..core.saveable import Saveable
from .opensystem import OpenSystem

from .. import REAL

 
    
class Molecule(UnitsManaged, Saveable, OpenSystem):
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

    # These will be inherited from OpenSystem
    
    # energies of electronic states
    #elenergies = array_property('elenergies')
    # transition dipole moments
    #dmoments = array_property('dmoments')    
    
    # number of electronic statesarray_property
    nel      = Integer('nel')
    # number of allowed transitions
    nal      = Integer('nal')
    # number of vibrational modes
    nmod     = Integer('nmod')
    
    
    def __init__(self,  elenergies=[0.0,1.0], name=None): #,dmoments):
    
        OpenSystem.__init__(self, elenergies)
        #self.manager = Manager()
        
        # monomer name
        if name is None:
            # FIXME: generate unique name
            self.name = ""
        else:
            self.name = name  #
            
        #
        # set energies
        # 

        # convert to internal_units
        #self.elenergies = self.manager.convert_energy_2_internal_u(elenergies)
        energy = []
        self.Nb = numpy.zeros(len(elenergies),dtype=numpy.int64)
        self.which_band = []
        for n in range(self.Nb.size):
            if isinstance(elenergies[n],list):
                self.Nb[n] = len(elenergies[n])
                for ii in range(self.Nb[n]):
                    self.which_band.append(n)
                    energy.append(elenergies[n][ii])
            else:
                self.Nb[n] = 1
                self.which_band.append(n)
                energy.append(elenergies[n])
        self.which_band = numpy.array(self.which_band,dtype=numpy.int64)
        
        self.elenergies = energy
        self.elenergies = self.convert_energy_2_internal_u(self.elenergies) #
        self.nel = len(self.elenergies)    #
        
        
        # FIXME: check the order of energies (increasing order has 
        # to be enforced) no vibrational modes is a default
        self.nmod = 0
        self.modes = []  
        
        # allowed transitions are now only between the ground state and
        # the rest of the excited states are dark
        self.allowed_transitions = []
         
        # transition dipole moments
        self.dmoments = numpy.zeros((self.nel,self.nel,3))
        
        # transition velocity dipole moments
        self.dvmoments = numpy.zeros((self.nel,self.nel,3),dtype=numpy.complex128)
        
        # transition magnetic dipole moments
        self.mmoments = numpy.zeros((self.nel,self.nel,3),dtype=numpy.complex128)
        
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
        self._has_transition_velocity = False
        
        #
        self._adiabatic_initialized = False
 
        # we can specify diabatic coupling matrix
        self._diabatic_initialized = False

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
        
        # how many optical bands the molecule has
        # the band is defined as a set of states separated from 
        # another band by optical frequency.
        # Currently we always have just ONE band, although we might have
        # more states in it
        self.mult = 1
        
        # here we will store the Hamiltonian, once it is constructed
        self.HH = None
        
        self.el_rwa_indices = None
        self.has_rwa = False
        
        self.build()
        

    def build(self):
        """Building routine for the molecule
        
        
        Unlike with the Aggregate, it is not compulsory to call build()
        before we start using the Molecule
        
        """
        
        self.Nel = self.nel
        self._built = True


    #
    # I am systematically removing any "naming" in Quantarhei
    #
        
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
        
    def set_electronic_rwa(self, rwa_indices):
        """Sets which electronic states belong to different blocks
        
        Setting the indices of blocks should allow the construction
        of Rotating Wave Approximation for the Hamiltonian. This in turn
        is used for the calculations of optical spectra.
        
        
        Examples
        --------
        
        >>> with energy_units("1/cm"):
        ...     mol1 = Molecule([0.0, 10000.0])
        >>> mol1.set_electronic_rwa([0, 1])
        >>> H = mol1.get_Hamiltonian()
        >>> H.has_rwa
        True
        
        
        >>> with energy_units("1/cm"):
        ...     mol2 = Molecule([0.0, 10000.0, 10500.0, 20000.0])

        >>> with energy_units("1/cm"):
        ...     mod2 = Mode(frequency=200.0)
        ...     mol2.add_Mode(mod2) 
        >>> mod2.set_nmax(0, 5)
        >>> mod2.set_nmax(1, 4)
        >>> mod2.set_nmax(2, 5)
        >>> mod2.set_nmax(3, 3)
                   
        >>> mol2.set_electronic_rwa([0, 1, 3])
        >>> H = mol2.get_Hamiltonian()
        >>> H.has_rwa
        True
        
        >>> print(H.rwa_indices)
        [ 0  5 14]
        
        """
        
        self.el_rwa_indices = rwa_indices
        if rwa_indices is not None:
            self.has_rwa = True
        else:
            self.has_rwa = False
        
        if self.HH is not None:
            vib_rwa_indices = self._calculate_rwa_indices(rwa_indices)
            self.HH.set_rwa(vib_rwa_indices)
            
        
        
    def _calculate_rwa_indices(self, el_rwa_indices):
        """Extends the rwa_indices by vibrational states if necessary
        
        """
        if self.nmod == 0:
            # if there are no modes, than we can use 
            # electronic rwa_indices directly
            rwa_indices = el_rwa_indices

        else:
            # if modes are present, we get the number of vib states
            # in each block
            
            rwa_indices=[None]*len(el_rwa_indices)
            
            ib = 0 # block index
            rwa_indices[ib] = 0
            
            ib += 1  
            next_el_index = el_rwa_indices[ib] # next is the start of the 
            count = 0
            
            # loop over electronic states
            for kk in range(self.nel):
                
                inst_c = 1 # at leat one state per electronic state

                # if this is the start of the new block
                if kk == next_el_index:
                    
                    rwa_indices[ib] = rwa_indices[ib-1] + count
                    
                    ib += 1 # next block
                    if ib < len(el_rwa_indices):
                        next_el_index = el_rwa_indices[ib]
                    else:
                        next_el_index = self.nel + 1
                    count = 0 # reset the count of states

                # number of states in an electronic state
                for mk in range(self.nmod):  # over all modes
                    # number of vib. states
                    inst_c = inst_c*self.modes[mk].get_nmax(kk)
                    
                count += inst_c
            
        return rwa_indices
        
        
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
            
            
            
        Examples
        --------

        A set of three monomers
 
        >>> en1 = [0.0,12100, 13000] #*cm2int]
        >>> en2 = [0.0,12000] #*cm2int]
        >>> with energy_units("1/cm"):
        ...     m1 = Molecule(en1)
        ...     m2 = Molecule(en1)
        ...     m3 = Molecule(en2)

        Bath correlation functions to describe molecular environment
        
        >>> from .. import TimeAxis
        >>> from .. import CorrelationFunction        
        >>> time = TimeAxis(0.0, 2000, 1.0) # in fs
        >>> temperature = 300.0 # in Kelvins
        >>> cfce_params1 = dict(ftype="OverdampedBrownian",
        ...           reorg=30.0,
        ...           cortime=60.0,
        ...           T=temperature)
        >>> cfce_params2 = dict(ftype="OverdampedBrownian",
        ...           reorg=30.0,
        ...           cortime=60.0,
        ...           T=temperature)
        >>> with energy_units("1/cm"):
        ...    cf1 = CorrelationFunction(time,cfce_params1)
        ...    cf2 = CorrelationFunction(time,cfce_params2)
        
        Environment of the molecules is collected to a matrix. A smaller
        number of correlation functions can be assigned to a large
        number of "sites".
        
        >>> cm = CorrelationFunctionMatrix(time,3)
        >>> ic1 = cm.set_correlation_function(cf1,[(0,0),(2,2)])
        >>> ic2 = cm.set_correlation_function(cf2,[(1,1)])

        The sites of the in the matrix are assigned to molecules. 
        
        >>> m1.set_egcf_mapping((0,1), cm, 0)
        >>> m2.set_egcf_mapping((0,1), cm, 1)
        >>> m3.set_egcf_mapping((0,1), cm, 2)
        
        The environment cannot be set twice
        
        >>> m1.set_egcf_mapping((0,1), cm, 2)
        Traceback (most recent call last):
            ...
        Exception: Monomer has a correlation function already
        
        >>> m1.set_egcf_mapping((0,2), cm, 2)
        >>> m1.set_egcf_mapping((0,2), cm, 2)
        Traceback (most recent call last):
            ...
        Exception: Monomer has a correlation function already
        
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
        

    def set_transition_environment(self, transition, egcf):
        """Sets a correlation function for a transition on this monomer
        
        Parameters
        ----------
        transition : tuple
            A tuple describing a transition in the molecule, e.g. (0,1) is
            a transition from the ground state to the first excited state.

        egcf : CorrelationFunction
            CorrelationFunction object 
            
            
        Example
        -------
        
        >>> from ..qm.corfunctions import CorrelationFunction
        >>> from .. import TimeAxis 
        >>> ta = TimeAxis(0.0,1000,1.0)
        >>> params = dict(ftype="OverdampedBrownian",reorg=20,cortime=100,T=300)
        >>> with energy_units("1/cm"):
        ...     cf = CorrelationFunction(ta, params)
        >>> m = Molecule([0.0, 1.0])
        >>> m.set_transition_environment((0,1), cf)
        >>> print(m._has_system_bath_coupling)
        True
        
        When the environment is already set, the next attempt is refused
        
        >>> m.set_transition_environment((0,1), cf)
        Traceback (most recent call last):
        ...
        Exception: Correlation function already speficied for this monomer


        The environment cannot be set when the molecule is mapped on 
        a correlation function matrix

        >>> with energy_units("1/cm"):
        ...     cf1 = CorrelationFunction(ta, params)
        >>> cm = CorrelationFunctionMatrix(ta,1)
        >>> ic1 = cm.set_correlation_function(cf1,[(0,0)])
        >>> m1 = Molecule([0.0, 1.0])
        >>> m1.set_egcf_mapping((0,1), cm, 0)
        >>> m1.set_transition_environment((0,1), cf1)
        Traceback (most recent call last):
            ...
        Exception: This monomer is mapped on a CorrelationFunctionMatrix
        
        
        """
        if self._is_mapped_on_egcf_matrix:
            raise Exception("This monomer is mapped" 
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
 
        Example
        -------
        
        >>> from ..qm.corfunctions import CorrelationFunction
        >>> from .. import TimeAxis 
        >>> ta = TimeAxis(0.0,1000,1.0)
        >>> params = dict(ftype="OverdampedBrownian",reorg=20,cortime=100,T=300)
        >>> with energy_units("1/cm"):
        ...     cf = CorrelationFunction(ta, params)
        >>> m = Molecule([0.0, 1.0])
        >>> m.set_transition_environment((0,1), cf)
        >>> print(m._has_system_bath_coupling)
        True
        
        >>> m.unset_transition_environment((0,1))
        >>> print(m._has_system_bath_coupling)
        False
        
        When the environment is unset, the next attempt to set is succesful
        
        >>> m.set_transition_environment((0,1), cf)
        >>> print(m._has_system_bath_coupling)
        True


        The environment cannot be unset when the molecule is mapped on 
        a correlation function matrix

        >>> with energy_units("1/cm"):
        ...     cf1 = CorrelationFunction(ta, params)
        >>> cm = CorrelationFunctionMatrix(ta,1)
        >>> ic1 = cm.set_correlation_function(cf1,[(0,0)])
        >>> m1 = Molecule([0.0, 1.0])
        >>> m1.set_egcf_mapping((0,1), cm, 0)
        >>> m1.unset_transition_environment((0,1))
        Traceback (most recent call last):
            ...
        Exception: This monomer is mapped on a CorrelationFunctionMatrix

            
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


    #@deprecated
    def set_egcf(self, transition, egcf):
        self.set_transition_environment(transition, egcf)
        
            
    def get_transition_environment(self, transition):
        """Returns energy gap correlation function of a monomer
        
        Parameters
        ----------
        transition : tuple
            A tuple describing a transition in the molecule, e.g. (0,1) is
            a transition from the ground state to the first excited state.
 
    
        Example
        -------
        

        >>> m = Molecule([0.0, 1.0])
        
        Environment of the transition has to be set first
        
        >>> cc = m.get_transition_environment((0,1))
        Traceback (most recent call last):
            ...
        Exception: No environment set for the transition


        Environment is characterized by the bath correlation function
        
        >>> from ..qm.corfunctions import CorrelationFunction
        >>> from .. import TimeAxis 
        >>> ta = TimeAxis(0.0,1000,1.0)
        >>> params = dict(ftype="OverdampedBrownian",reorg=20,cortime=100,T=300)
        >>> with energy_units("1/cm"):
        ...     cf = CorrelationFunction(ta, params)        
        >>> m.set_transition_environment((0,1), cf)
        >>> cc = m.get_transition_environment((0,1))
        >>> cc == cf
        True
        
        When non-existent transition is tried, exception is raised
        
        >>> cc = m.get_transition_environment((0,2))
        Traceback (most recent call last):
            ...        
        Exception: Index out of range
        
        
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
    
        >>> mode = Mode()
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
            
    
    def get_velocity_dipole(self, N, M):
        try:
            if self._has_transition_velocity:
                return self.dvmoments[N, M, :]
            else:
                raise Exception()
        except:
            raise Exception()
            
            
    def get_magnetic_dipole(self, N, M):
        try:
            return self.mmoments[N, M, :]
        except:
            raise Exception()
            
            
    def set_dipole(self, N, M, vec=None):
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
    
    def set_velocity_dipole(self, N, M, vec=None):

        if vec is None:
            n = N[0]
            m = N[1]
            vc = M
        else:
            n = N
            m = M
            vc = vec


        # FIXME: use complex dipole velocity moment (pure imaginary quantity)
        if n == m:
            raise Exception("M must not be equal to N")
        try:
            self.dvmoments[n, m, :] = vc
            self.dvmoments[m, n, :] = numpy.conj(vc)
            self._has_transition_velocity = True
        except:
            raise Exception()  
            
    def set_velocity_dipole_from_dipole(self):
        # FIXME: use complex dipole velocity moment (pure imaginary qunatity)
        #try:
            for N in range(self.dmoments.shape[0]):
                for M in range(N+1,self.dmoments.shape[1]):
                    DE = self.elenergies[M] - self.elenergies[N]
                    self.dvmoments[N, M, :] = -1j*self.dmoments[N,M,:]*DE
                    self.dvmoments[M, N, :] = numpy.conj(self.dvmoments[N, M, :])
            self._has_transition_velocity = True
        #except:
            #raise Exception()  
            
    def set_magnetic_dipole(self, N, M, vec=None):
        # vec is probably in atomic units and our units are [Angstrom*1/fs*Debye]
#        print(vec)

        if vec is None:
            n = N[0]
            m = N[1]
            vc = M
        else:
            n = N
            m = M
            vc = vec


        units = self.manager.get_current_units("energy")
        vec_int = self.manager.convert_energy_2_internal_u(vc) \
                     *conversion_facs_edipole[units] * conversion_facs_length[units]
        
        # magnetic dipole moment stored in internal units
        if n == m:
            raise Exception("M must not be equal to N")
        try:
            self.mmoments[n, m, :] = vec_int
            self.mmoments[m, n, :] = numpy.conj(vec_int)
        except:
            raise Exception()  
        
    def set_magnetic_dipoleR(self, N, M, vec, RR=None):


        if RR is None:
            n = N[0]
            m = N[1]
            vc = M
            R = vec
        else:
            n = N
            m = M
            vc = vec
            R = RR

        # So far vec in atomic units (complex) and R in Angstroms
        try:
            dv = self.dvmoments[n, m, :]
        except:
            DE = self.elenergies[m] - self.elenergies[n]
            dv = -1j*self.dmoments[n,m,:]*DE
        
        vc = numpy.array(vc,dtype=numpy.complex128)
        R = numpy.array(R,dtype="f8")
        
        # dipole velocity (dv) in Int*Debye -> Transform to a.u.
        units = self.manager.get_current_units("energy")
        dv_au =  self.manager.convert_energy_2_current_u(dv)/conversion_facs_edipole[units]
        
        # position is in angstroms -> transform to a.u.
        R_au = R/conversion_facs_length[units]
        
        # look if dipole velocity is defined:
        m_int_au = vc - numpy.cross(R_au,dv_au)
        
        #print(m_int_au)
        self.set_magnetic_dipole(n, m,m_int_au)
        

    # def set_dipole(self, N, M, vec=None):
    #     """Sets transition dipole moment for an electronic transition
        
        
    #     There are two ways how to use this function:
            
    #         1) recommended
                
    #             set_dipole((0,1),[1.0, 0.0, 0.0])
                
    #             here N represents a transition by a tuple, M is the dipole
                
    #         2) deprecated (used in earlier versions of quantarhei)
                
    #             set_dipole(0,1,[1.0, 0.0, 0.0])
                
    #             here the transition is characterized by two integers
    #             and the last arguments is the vector
                
        
        
    #     Examples
    #     --------
 
    #     >>> m = Molecule([0.0, 1.0])
    #     >>> m.set_dipole((0,1),[1.0, 0.0, 0.0])
    #     >>> m.get_dipole((0,1))
    #     array([ 1.,  0.,  0.])

        
    #     """
    #     if vec is None:
    #         n = N[0]
    #         m = N[1]
    #         vc = M
    #     else:
    #         n = N
    #         m = M
    #         vc = vec
            
    #     if n == m:
    #         raise Exception("M must not be equal to N")
    #     try:
    #         self.dmoments[n, m, :] = vc
    #         self.dmoments[m, n, :] = numpy.conj(vc)
    #     except:
    #         raise Exception()
            
        
    # def get_dipole(self, N, M=None):
    #     """Returns the dipole vector for a given electronic transition
        
    #     There are two ways how to use this function:
            
    #         1) recommended
                
    #             get_dipole((0,1),[1.0, 0.0, 0.0])
                
    #             here N represents a transition by a tuple, M is the dipole
                
    #         2) deprecated (used in earlier versions of quantarhei)
                
    #             get_dipole(0,1,[1.0, 0.0, 0.0])
                
    #             here the transition is characterized by two integers
    #             and the last arguments is the vector        
                
        
    #     Examples
    #     --------
        
    #     >>> m = Molecule([0.0, 1.0])
    #     >>> m.set_dipole((0,1),[1.0, 0.0, 0.0])
    #     >>> m.get_dipole((0,1))
    #     array([ 1.,  0.,  0.])
        
        
    #     """
    #     if M is None:
    #         n = N[0]
    #         m = N[1]
    #     else:
    #         n = N
    #         m = M
            
    #     try:
    #         return self.dmoments[n, m, :]
    #     except:
    #         raise Exception()

            
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
        cwidth = Manager().convert_energy_2_internal_u(width)
        if self.widths is None:
            N = self.elenergies.shape[0]
            self.widths = numpy.zeros((N, N), dtype=REAL)
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
            self.dephs = numpy.zeros((N, N), dtype=REAL)
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
        
        
        Examplex
        --------
        
        Molecule where transition dipole was not set does not calculate lifetime
        
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
    
        
    def get_temperature(self):
        """Returns temperature of the molecule
        
        Checks if the setting of environments is consistent and than
        takes the temperature from one of the energy gaps. If no 
        environment (correlation function) is assigned to this 
        molecule, we assume zero temperature.
            
        Examples
        --------
        
        Default temperature is 0 K
        
        
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mol.get_temperature()
        0.0
        
        Molecule gets its temperature from the environment

        >>> from ..qm.corfunctions import CorrelationFunction
        >>> from .. import TimeAxis 
        >>> ta = TimeAxis(0.0,1000,1.0)
        >>> params = dict(ftype="OverdampedBrownian",reorg=20,cortime=100,T=300)
        >>> with energy_units("1/cm"):
        ...     cf = CorrelationFunction(ta, params)
        >>> m = Molecule([0.0, 1.0])
        >>> m.set_transition_environment((0,1), cf)
        >>> print(m.get_temperature())
        300
        
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
            " an inconsisten temperature")
        
        
    def check_temperature_consistent(self):
        """Checks that the temperature is the same for all components
 
    
        Examples
        --------

        Isolated molecule has always consistent temperature
        
        >>> m = Molecule([0.0, 1.0, 1.2])
        >>> print(m.check_temperature_consistent())
        True
        

        >>> from ..qm.corfunctions import CorrelationFunction
        >>> from .. import TimeAxis 
        >>> ta = TimeAxis(0.0,1000,1.0)
        >>> params1 = dict(ftype="OverdampedBrownian",
        ...                reorg=20,cortime=100,T=300)
        >>> params2 = dict(ftype="OverdampedBrownian",
        ...                reorg=20,cortime=100,T=200)
        >>> with energy_units("1/cm"):
        ...     cf1 = CorrelationFunction(ta, params1)
        ...     cf2 = CorrelationFunction(ta, params2)
        >>> m.set_transition_environment((0,1), cf1)
        >>> m.set_transition_environment((0,2), cf2)
        >>> print(m.get_temperature())
        Traceback (most recent call last):
        ...
        Exception: Molecular environment has an inconsisten temperature

                
        """
        
        if self._has_system_bath_coupling:
            
            
            T = -10.0
            for bath in self.egcf:
                
                if bath is not None:
                    if T < 0.0:
                        T = bath.get_temperature()
                    elif T != bath.get_temperature():
                        return False
                
            return True
        
        else: 
            
            return True


   
    def set_diabatic_coupling(self, element, factor, shifts=None):
        """Sets off-diagobal elements of the diabatic potential matrix 
        
        """
        # FIXME: we ignore shifts for now
        Nm = self.get_number_of_modes()
        faclength = len(factor[1])
        
        # energy conversion
        val = self.convert_energy_2_internal_u(factor[0])
        factor[0] = val
        
        if  faclength != Nm:
            raise Exception("Expected "+str(Nm)+
                            " mode, found "+str(faclength)+".")
        
        if not self._diabatic_initialized:
            
            Ne = self.nel
            
            self.diabatic_matrix = []
            for ii in range(Ne):
                self.diabatic_matrix.append([])
                for jj in range(Ne):
                    self.diabatic_matrix[ii].append([])
            
            self._diabatic_initialized = True
        
        self.diabatic_matrix[element[0]][element[1]].append(factor)
        self.diabatic_matrix[element[1]][element[0]].append(factor)



    def _fill_hmatrix(self, HH, en0, coorval):
        """Creates Hamiltonian matrix for given values of the coordinates
        
        """
        
        for ii in range(self.nel):
            for jj in range(self.nel):
                
                HH[ii,jj] = 0.0
                
                if (ii==jj):
                    
                    HH[ii,jj] = self.elenergies[ii]
                    
                    km = 0
                    # loop over modes
                    for md in self.modes:
                        
                        qq = coorval[km]
                            
                        HH[ii,ii] += (md.get_energy(ii)/2.0)* \
                                     (qq-md.get_shift(ii))**2
                                     
                        km += 1
                        
                    en0[ii] = HH[ii,ii]
                    
                else:
                    
                    with energy_units("int"):
                        coup = self.get_diabatic_coupling((ii,jj))
                        lc = len(coup)
                        
                        if lc > 0:
                            # loop over couplings
                            for ll in range(lc):
                                val = coup[ll]

                                if len(val) > 0:
                                    
                                    # loop over modes
                                    for lm in range(self.nmod):
                                        qq = coorval[lm]
                                        # we use linear contribution only
                                        if val[1][lm]==1:
                                            
                                            HH[ii,jj] += val[0]*qq
 

    def get_potential_1D(self, mode, points, other_modes=None):
        """Returns the one dimensional diabatic potentials
        
        """
        
        if other_modes is None:
            # default value for other than the plotted mode
            other_modes = [0.0]*self.nmod
            
        if len(other_modes) != self.nmod:
            raise Exception("Argument 'other_modes' has to have the lenth"+
                            " equal to the number of modes")
            
        coorval = numpy.zeros(self.nmod)
        HH = numpy.zeros((self.nel, self.nel), dtype=REAL)
        pot = numpy.zeros((len(points),self.nel), dtype=REAL)
        pot0 = numpy.zeros((len(points),self.nel), dtype=REAL)
        
        
        # loop over a single coordinate
        kk = 0
        for pt in points:
            
            for ii in range(self.nmod):
                if ii == mode:
                    coorval[ii] = pt
                else:
                    coorval[ii] = other_modes[ii]
            
            self._fill_hmatrix(HH, pot0[kk,:], coorval)

            (en, ss) = numpy.linalg.eigh(HH)
            pot[kk,:] = en
            kk += 1
           
        return (pot, pot0)
                
           
    def get_potential_2D(self, modes, points, other_modes=None):
        """Returns the two dimensional diabatic potentials
        
        """
        if other_modes is None:
            # default value for other than the plotted mode
            other_modes = [0.0]*self.nmod
            
        if len(other_modes) != self.nmod:
            raise Exception("Argument 'other_modes' has to have the lenth"+
                            " equal to the number of modes")
            
        coorval = numpy.zeros(self.nmod)
        HH = numpy.zeros((self.nel, self.nel), dtype=REAL)
        pot = numpy.zeros((len(points[0]),len(points[1]),self.nel),
                          dtype=REAL)
        pot0 = numpy.zeros((len(points[0]),len(points[1]),self.nel),
                          dtype=REAL)
        
        # loop over two coordinates
        k1 = 0
        for pt1 in points[0]:
            k2 = 0
            for pt2 in points[1]:
                
                for ii in range(self.nmod):
                    if ii == modes[0]:
                        coorval[ii] = pt1
                    elif ii == modes[1]:
                        coorval[ii] = pt2
                    else:
                        coorval[ii] = other_modes[ii]   
                        
                self._fill_hmatrix(HH, pot0[k1,k2,:], coorval)
                    
                (en, ss) = numpy.linalg.eigh(HH)
                pot[k1,k2,:] = en
                
                k2 += 1
            k1 += 1                    
                    
        return (pot, pot0)
           


    def plot_potential_1D(self, mode, points, other_modes=None, 
                          nonint=True, states=None,
                          energies=True, show=True,
                          ylims=None):
        """Plots the potentials
        
        """
        
        import matplotlib.pyplot as plt
        
        pot, pot0 = self.get_potential_1D(mode, points, 
                                          other_modes=other_modes)
        
        if states is None:
            sts = []
            for ii in range(self.nel):
                sts.append(ii)
        else:
            sts = states
            
        for ss in sts:
            plt.plot(points, pot[:,ss])
            if nonint:
                plt.plot(points,pot0[:,ss],"--")
           
        if energies:
            HH = self.get_Hamiltonian()
            with eigenbasis_of(HH):
                for ii in range(HH.dim):
                    enr = HH.data[ii,ii]*numpy.ones(len(points), dtype=REAL)
                    plt.plot(points, enr, "-k")
            if nonint:
                for ii in range(HH.dim):
                    enr = HH.data[ii,ii]*numpy.ones(len(points), dtype=REAL)
                    plt.plot(points, enr, "--b")
            
        if ylims is not None:
            plt.ylim(ylims[0],ylims[1])
        if show:
            plt.show()
                

    def plot_stick_spectrum(self, xlims=[0.0,1.0], ylims=None, 
                            show_zero_coupling=False, show=True):
        """Plots the stick spectrum of the molecule 
        
        """
            
        import matplotlib.pyplot as plt
        
        HH = self.get_Hamiltonian()
        dip = self.get_TransitionDipoleMoment()
        
        
        with eigenbasis_of(HH):
            y1 = (dip.data[0,:,0]**2+dip.data[0,:,1]**2+dip.data[0,:,2]**2)/3.0
            x1 = numpy.diag(HH.data)
    
        plt.stem(x1,y1,markerfmt=' ', linefmt='-r')
    
        if show_zero_coupling:
            y0 = (dip.data[0,:,0]**2+dip.data[0,:,1]**2+dip.data[0,:,2]**2)/3.0
            x0 = numpy.diag(HH.data)
    
            plt.stem(x0,y0,markerfmt=' ', linefmt="-b")
            
        plt.xlim(xlims[0],xlims[1])
        if show:
            plt.show()


    def plot_dressed_sticks(self, dfce=None, xlims=[0,1], nsteps=1000,
                            show_zero_coupling=False, show=True):
        """Plots a stick spectrum dessed by a supplied function
        
        
        """

        import matplotlib.pyplot as plt        

        if dfce is None:
            raise Exception("Dressing function must be defined")
        
        HH = self.get_Hamiltonian()
        dip = self.get_TransitionDipoleMoment()
 
        dx = (xlims[1]-xlims[0])/nsteps
        xe = numpy.array([xlims[0] + i*dx for i in range(nsteps)])
        sp1 = numpy.zeros(len(xe), dtype=REAL)
        
        with eigenbasis_of(HH):
            y1 = (dip.data[0,:,0]**2+dip.data[0,:,1]**2+dip.data[0,:,2]**2)/3.0
            x1 = numpy.diag(HH.data)  
            
        for kk in range(HH.dim):
            sp1 += y1[kk]*dfce(xe, x1[kk])

        plt.plot(xe, sp1, "-r")
    
        if show_zero_coupling:
            y0 = (dip.data[0,:,0]**2+dip.data[0,:,1]**2+dip.data[0,:,2]**2)/3.0
            x0 = numpy.diag(HH.data)
            sp0 = numpy.zeros(len(xe), dtype=REAL)
        
            for kk in range(HH.dim):
                sp0 += y0[kk]*dfce(xe, x0[kk])
            
            plt.plot(xe, sp0, "-b")
        
        
        plt.xlim(xlims[0],xlims[1])
        if show:
            plt.show()



    def get_diabatic_coupling(self, element):
        """Returns list of coupling descriptors
             
        """
        if self._diabatic_initialized:
            # FIXME: only first factor used
            factor = self.diabatic_matrix[element[0]][element[1]]
            ven = []
            if len(factor) > 0:
                for fc in factor:
                    fcv = fc[0]
                    val = self.convert_energy_2_current_u(fcv)
                    ven.append([val, fc[1]])
                return  ven
            else:
                return []
        else:
            return []
    
    
    def get_diabatic_shifts(self, order=1):
        
        raise Exception("Shifts not implemented")






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
           
    def _overlap_other(self, tpl1, tpl2, k):
        dif = 0
        for i in range(len(tpl1)):
            if i != k:
                dif += numpy.abs(tpl1[i]-tpl2[i])
                
        if dif == 0:
            return 1.0
        else:
            return 0.0
        

    def _overlap_all(self, tpl1, tpl2):
        dif = 0
        for i in range(len(tpl1)):
            dif += numpy.abs(tpl1[i]-tpl2[i])
                
        if dif == 0:
            return 1.0
        else:
            return 0.0 

    
    
           
    def get_Hamiltonian(self, multi=True, recalculate=False):
        """Returns the Hamiltonian of the Molecule object
        
        
        Examples
        --------
        
        After the molecule is created, its Hamiltonian can be obtained
        
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
        
        
        The Hamiltonian is stored, and next time it is retrieved
        we obtain the same object.
        
        >>> print(H.has_rwa)
        False
        
        We can manipulate with the Hamiltonian in between retrievals
        >>> H.set_rwa([0,1])
        >>> print(H.has_rwa)
        True
        
        The Hamiltonian the obtain the second time is affected by the changes
        performed outside the Molecule class
        
        >>> H1 = mol.get_Hamiltonian()
        >>> print(H1.has_rwa)
        True
        
        The newly obtained Hamiltonina IS the Hamiltonian obtained earlier.
        
        >>> H1 is H
        True
        
        """    
    
        if (self.HH is not None) and (not recalculate):
            
            return self.HH
    
        if multi:
            
            # create vibrational signatures for each electronic state
            vsignatures = []
            for kk in range(self.nel):
                vibmax = []
                for mk in range(self.nmod):
                    vibmax.append(self.modes[mk].get_nmax(kk))

                signatures = numpy.ndindex(tuple(vibmax))
                vsignatures.append(signatures)
    
            # the list of vibrational signatures 
            self.vibsignatures = vsignatures
    
    
            # list o signatures of all states
            self.all_states = []
            ks = 0
            ke = 0
            for vsig_it in vsignatures:
                for vsig in vsig_it:
                    elvibstate = tuple([ke, vsig])
                    self.all_states.append(elvibstate)
                    ks += 1
                ke += 1
                    
            # total number of states
            self.totstates = ks
              
            #
            # building the Hamiltonian
            #
            ham = numpy.zeros((ks,ks), dtype=REAL)
            
            
            
            
            # FIXME: creation of the coordinate operators will go to 
            # the place where SystemBathInteraction is created
            #
            # coordinate operator for each mode
            #
            coor_ops = []
            for kk in range(self.nel):
                in_state = []
                coor_ops.append(in_state)
                for ii in range(self.nmod):
                    coor = numpy.zeros((ks,ks), dtype=REAL)
                    in_state.append(coor)
                    #print("State:", kk," - mode:", ii, coor.shape)
                    
                    
                
            #
            # for each state and mode, we create vibrational Hamiltonian
            #
            hh_components = []
            qq_components = []
            # loop over electronic states
            for i in range(self.nel):
                el_state = []
                qq_state = []
                hh_components.append(el_state)
                qq_components.append(qq_state)
                if self.nmod > 0:
                    # loop over modes
                    for j in range(self.nmod):
                        
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
                    
                        qq = (1.0/numpy.sqrt(2.0))*(ad+aa) 
                        
                        hh = en*(numpy.dot(ad,aa) 
                                 - dd*qq
                                 + dd*dd*ones/2.0) + (en/2.0)*ones
                        
                        el_state.append(hh)
                        
                        qq = qq - dd*ones
                        qq_state.append(qq)
                       
                    
                # if there are no modes
                else:
                    hh = numpy.zeros((1,1),dtype=REAL)
                    hh[0,0] = self.elenergies[i]
                    
    
                    el_state.append(hh)
                    
                    
 
            # case - NO MODES
            # FIXME: faster code for the no modes case

                
            
            
            # case - AT LEAST ONE MODE
            
            # loop over all states
            ks1 = 0
            for st1 in self.all_states:
                n = st1[0]
                vibn = st1[1]
                
                # loop over all stataes
                ks2 = 0
                for st2 in self.all_states:
                    m = st2[0]
                    vibm = st2[1]
                    
                    #
                    # electronic states
                    #
                    if n == m:
                        
                        el_state = hh_components[n]
                        qq_state = qq_components[n]
                        
                        # electronic part of the energy
                        if ks1 == ks2:
                            ham[ks1, ks2] += self.elenergies[n]
                    
                        for k in range(self.nmod):
                            overl = self._overlap_other(vibn,vibm,k)
                            
                            hh = el_state[k]
                            kn = vibn[k]
                            km = vibm[k]
                            ham[ks1, ks2] += hh[kn,km]*overl
                            
                            #
                            # coordinate operators
                            #
                            qq = qq_state[k]
                            coor = coor_ops[n][k]
                            coor[ks1, ks2] += qq[kn,km]*overl
                            
                    #   
                    # coupling elements
                    #
                    else:
                        
                        if self._diabatic_initialized:
                            dmx = self.diabatic_matrix[n][m]
                            # number of defined couplings
                            Ncoup = len(dmx)
                            if Ncoup > 0:
                                # loop over the coupling definitions
                                for nc in range(Ncoup):
                                    modc = dmx[nc]
                                    val = modc[0]   # coupling constant
                                    coors = modc[1] # which modes contribute
                                    # loop over modes
                                    for ci in range(self.nmod):
                                        
                                        # other modes than ci have to be
                                        # in the same states 
                                        overl = self._overlap_other(vibn,
                                                                   vibm,ci)
                                        # FIXME: bilinear coupling
                                        # this prevents bilinear coupling
                                        if overl == 1.0:
                                            
                                            # FIXME: we only allow linear
                                            #        contribution
                                            if coors[ci] == 1:
                                                if vibn[ci] == vibm[ci]+1:
                                                    ham[ks1, ks2] += \
                                                  val*numpy.sqrt(vibn[ci]/2.0)
                                                elif vibn[ci] == vibm[ci]-1:
                                                    ham[ks1, ks2] += \
                                                  val*numpy.sqrt(vibm[ci]/2.0)
                        
                    
                    ks2 += 1
                    
                ks1 += 1
    
            #
            # here we store all the coordinate operators
            #
            self.coor_operators = coor_ops
            
            if (self.HH is None) or recalculate:
                with energy_units("int"):
                    HH = Hamiltonian(data=ham)
                    
                    # set ground state energy to zero
                    with eigenbasis_of(HH):
                        e00 = HH.data[0,0]
                        
                    HH.data -= numpy.eye(HH.dim)*e00
                    
                    if self.has_rwa:
                        # set Hamiltonian to RWA
                        vib_rwa_indices = \
                            self._calculate_rwa_indices(self.el_rwa_indices)

                        HH.set_rwa(vib_rwa_indices)
                    
                self.HH = HH  
                
                # set information about Rotating Wave Approximation
                
            return self.HH
    
  
###############################################################################
# old single mode version - will be removed in future
###############################################################################
   
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
                hh = numpy.zeros((1,1),dtype=REAL)
                hh[0,0] = self.elenergies[i]
                lham[i] = hh
                ldim[i] = 1
            
        # dimension of the complete Hamiltonian
        totdim = numpy.sum(ldim)
        
        # this will be the Hamiltonian data
        ham = numpy.zeros((totdim,totdim),dtype=REAL)
        
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
                        hj = numpy.zeros((ldim[i],ldim[j]),dtype=REAL)
                        # FIXME: this works only if the frequencies
                        # of the oscillators are the same
                        for k in range(min([ldim[i],ldim[j]])):
                            hj[k,k] = J
                        ham[lb[i]:ub[i],lb[j]:ub[j]] = hj
                        ham[lb[j]:ub[j],lb[i]:ub[i]] = hj.T
                        
                        
        #
        # diabatic coupling matrix
        #
        if self._diabatic_initialized:
            for i in range(self.nel):
                for j in range(self.nel):
                    if i != j:
                        
                        coups = self.diabatic_matrix[i][j]
                        
                        # FIXME
                        # we accept only the first element of the record
                        if len(coups) > 0:
                            rec = coups[0]
                            val = rec[0]
                            modes = rec[1]
                        
                            # FIXME
                            # we work with one mode only
                            if len(modes) == 1:
                                
                                # FIXME
                                # we allow ony linear dependence
                                if modes[0] == 1:
                                    n = 0
                                    for a in range(lb[i],ub[i]):
                                        m = 0
                                        for b in range(lb[j],ub[j]):
                                            if n == m+1:
                                                ham[a,b] = \
                                                    val*numpy.sqrt(n/2.0)
                                            if n == m-1:
                                                ham[a,b] = \
                                                    val*numpy.sqrt((n+1)/2.0)
                                        
                                            m += 1
                                        n += 1
                        
            
        return Hamiltonian(data=ham)
    """   
        
    def _sub_matrix_bounds(self,ldim):
        lbound = 0
        ub = numpy.zeros(self.nel,dtype=int)
        lb = numpy.zeros(self.nel,dtype=int)
        # loop over electronic states
        for i in range(self.nel):
            ubound = lbound + ldim[i]
            #ham[lbound:ubound,lbound:ubound] = lham[i]
            ub[i] = ubound  
            lb[i] = lbound
            lbound = ubound   
        return lb,ub
        
    def _ham_dimension(self):
        """Returns internal information about the dimension of the Hamiltonian
        
        Works only for one mode per molecule.
        Will be removed when multi mode molecule is fully implemented
        
        """

        # list of Hamiltonian dimensions
        ldim = [None]*self.nel

        # loop over electronic states
        for i in range(self.nel):
        
            Nvib = 1
            
            # loop over modes
            for j in range(self.nmod):
                # FIXME: enable more than one mode
                #if j > 0: # limits the number of modes to 1
                #    raise Exception("Not yet implemented") 
                    
                # number of vibrational states in this electronic state
                Nvib = Nvib*self.modes[j].get_nmax(i)
                
            ldim[i] = Nvib
            
        # dimension of the complete Hamiltonian
        totdim = numpy.sum(ldim)

        return totdim, ldim        
               
    
        
    def get_TransitionDipoleMoment(self, multi=True):
        """Returns the transition dipole moment operator
        
        """

        if multi:
            
            try:
                totdim = self.totstates
            except:
                HH = self.get_Hamiltonian()
                totdim = HH.dim

            # This will be the operator data
            dip = numpy.zeros((totdim,totdim,3),dtype=REAL)  
            
            ks1 = 0
            for st1 in self.all_states:
                n = st1[0]
                vibn = st1[1]
                
                ks2 = 0
                for st2 in self.all_states:
                    m = st2[0]
                    vibm = st2[1]  
                    
                    dp = self.dmoments[n,m,:]
                    ovrl = self._overlap_all(vibn,vibm)
                    
                    if ovrl > 0.0:
                        dip[ks1, ks2,:] = dp
                    
                    ks2 += 1
                    
                ks1 += 1
                
            
            return TransitionDipoleMoment(data=dip)

        #
        # older single mode version
        #
        """
        totdim,ldim = self._ham_dimension()

        # This will be the operator data
        dip = numpy.zeros((totdim,totdim,3),dtype=REAL) 
        
        lb,ub = self._sub_matrix_bounds(ldim)
        
        # FIXME: sofar only the lowest state of all is the start of 
        # optical transitions
        for k in range(1,self.nel):
            # FIXME: all just for one mode
            
            # get dipole moment vector
            dp = self.dmoments[0,k,:]
 
            dd = numpy.zeros((ldim[0],ldim[k]),dtype=REAL)        
        
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
    
        """
           
    
    def get_SystemBathInteraction(self): #, timeAxis):
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
                    # for each bath, save the state of the 
                    # transition g -> j
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
        #nmd = nob - ntr

        # number of different instances of correlation functions
        nof = cf.get_number_of_unique_elements()
        
        # number of different baths nob = number of transition environments +
        # number of mode environments
        cfm = None #
        uq = cf.get_unique_elements()
        for i in range(nof):
            el = uq[i] #cf.get_element(i)
            wr = where[el]
            
            if cfm is None:
                timeAxis = el.axis
                cfm = CorrelationFunctionMatrix(timeAxis,nob,nof)

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
            KK = numpy.zeros((totdim,totdim),dtype=REAL)

            state = d[n]

            states_before = 0
            for k in range(state):
                states_before += ldim[k]
            states_inc = states_before +ldim[state]
            # fill 1 on diagonal corresponding to an electronic state "state"
            KK[states_before:states_inc,
               states_before:states_inc] = numpy.diag(
                               numpy.ones(ldim[state],dtype=REAL))
            
            sys_operators.append(KK)

        
        #
        #  FIXME: mode environments must be hadled through Mode class
        #
        #
        #  Add mode interaction operators 
        #
        coor_ops = self.coor_operators
        if self._mode_env_initialized:
            
            for ii in range(self.nmod):
                
                #mod = self.modes[ii]
                
                # if the mode has environment, we set corresponding operators
                for nn in range(self.nel):
                    if self._has_mode_env[ii,nn]:
                        cop = coor_ops[nn][ii]
                        sys_operators.append(cop)
                
                        # get correlation function corresponding to 
                        # the system-bath interaction
                        cf = self._mode_env.get_element(ii, nn)
        
        sbi = SystemBathInteraction(sys_operators,cfm)
        
        return sbi
        
    
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
            raise Exception("Correlation function not specified.")
            
        if not self._mode_env_initialized:
            self._has_mode_env = numpy.zeros((self.nmod,self.nel),
                                             dtype=bool)
            self._mode_env = unique_array(self.nmod,self.nel)
            self._mode_env_initialized = True
            
        if isinstance(elstate, int):
            self._mode_env.set_element(mode,elstate,corfunc)
            self._has_mode_env[mode,elstate] = True
            
        #
        # if elstate is set to "ALL" we give the same environmemnt to all
        # electronic states
        #
        elif elstate == "all" or elstate == "ALL":
            for kk in range(self.nel):
               self._mode_env.set_element(mode,kk,corfunc)
               self._has_mode_env[mode,kk] = True
        else:
            raise Exception("Unknown elstate.")
        
        
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

                        # FIXME: what should be here???
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
        
