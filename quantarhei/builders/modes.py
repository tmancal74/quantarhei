# -*- coding: utf-8 -*-
"""
    Class representing intramolecular vibrations.
       
    The Mode class is a user class through which the user defines an intra
    molecular mode.
        
    The vibrational mode is supposed to be an intramolecular mode of some
    molecule. This class is therefore aware of its Molecule class. Once it
    knows in which Molecule object it lives, it creates instances of the class
    Submode (as many as there are electronic states in the Molecule). Submods 
    hold the parameters of the mode respective to a give electronic state
    of the monomer
    
    These parameters have to be set after the mode is registered in
    the Molecule, and therefore there is an issue of consistency of the class. 
    Currently, the consistency is completely in the hands of the user.
    
    Class Details
    -------------
    
    
"""

import numpy

from ..utils import Integer
from ..utils import Bool

from ..core.managers import UnitsManaged, energy_units
from ..core.wrappers import deprecated

from ..core.saveable import Saveable

from .submodes import SubMode
        
        
# Mode is UnitsManaged with components of different dimensions
# Mode is not BasisManaged
class Mode(UnitsManaged, Saveable):
    """ Vibrational mode
    
        
    Parameters
    ----------
    
    omega : float
        vibrational frequency
        
        
       
    Examples
    --------
    
    >>> import quantarhei as qr
    >>> mol = qr.Molecule([0.0, 1.0])
    >>> md = Mode(frequency=0.2)
    >>> mol.add_Mode(md)
    
    >>> print(md.get_energy(1))
    0.2
    
    This class is units management aware
    
    >>> import quantarhei as qr
    >>> with qr.energy_units("1/cm"):
    ...     md = Mode(250.0)
    >>> mol.add_Mode(md)
    >>> print(md.get_energy(0))
    0.047091289182721326
    
    >>> mol.get_number_of_modes()
    2
    
    """
    
    # number of electronic states (same as the monomer)
    nel = Integer('nel')      
    # is the monomer for this mode set?
    monomer_set = Bool('monomer_set')
    
    def __init__(self, frequency=1.0):
            
        # this holds a list of submods
        self.submodes = []
        # ground state submode is created (with shift = 0 
        # and a default no. of states)
        freq = self.convert_energy_2_internal_u(frequency)
        with energy_units("int"):
            self.submodes.append(SubMode(freq))
        
        # monomer is not set at creation
        self.monomer_set = False
        # no electronic states set or asigned
        self.nel = 0

        
    def set_Molecule(self, monomer):
        """Assigns this mode to a given monomer.
        
        When set, the mode knows on how many electronic states
        it is supposed to live. This method is called by the Molecule's
        `add_Mode` method.
        
        Parameters
        ----------
        
        monomer : quantarhei.Molecule
            Molecule object to which this Mode will be assigned
            
        
        Examples
        --------
        
        `set_Molecule` should not be called directly, except of some (hard
        to imagine) special cases. The `add_Mode` method of the Molecule
        class does some extra work to keep consistent record of the Modes in
        the Molecule. Here, the difference between `set_Molecule` and 
        `add_Mode` method of the Molecule class is demonstrated.
        
        >>> import quantarhei as qr
        >>> mol1 = qr.Molecule([0.0, 2.0])
        >>> mol2 = qr.Molecule([0.0, 2.0])
        >>> mod1 = qr.Mode(0.2)
        >>> mod2 = qr.Mode(0.2)
        
        >>> mol1.add_Mode(mod1)
        >>> mod2.set_Molecule(mol2)
        
        The number of modes recorded is consistent
        
        >>> print(mol1.get_number_of_modes())
        1
        >>> print(len(mod1.submodes))
        2
        
        The number of modes in this case is not consistent. The Mode knows that
        it has SubModes, but Molecule has no mode recorded.
        
        >>> print(mol2.get_number_of_modes())
        0
        >>> print(len(mod2.submodes))
        2
        
        
        
        """

        self.nel = monomer.nel
        self.monomer = monomer
        
        # Must have as many submodes as electronic states in the monomer
        with energy_units("int"):
            for k in range(1, self.nel):
                # submodes are created with the same frequency as
                # in the groundstate and with zero shift
                self.submodes.append(SubMode(self.submodes[0].omega))


        self.monomer_set = True
        
    #
    # Setters
    #
    
    @deprecated
    def set_frequency(self, N, omega):
        """Sets vibrational frequency
        
        Usage of this method is deprecated, use `set_energy` instead
        
        Parameters
        ----------
        
        N : int
            Index of the electronic state
            
        omega : float
            Vibrational frequency
            
        
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mod = mol.get_Mode(0)
        
        This function is deprecated and it throws a warning text when it is run
        
        >>> mod.set_frequency(1, 1.3) # doctest: +ELLIPSIS
        function  <function Mode.set_frequency at ...>  is deprecated
        
        >>> mod.get_energy(1)
        1.3
        
        """
        
        #sbm = self.submodes[N]
        #sbm.omega = omega
        self.set_energy(N, omega)

        
    def set_energy(self, N, omega):
        """Sets energy/frequency of the mode relative to a molecular state
        
        
        Parameters
        ----------
        
        N : int
            Index of the electronic state for we set frequency
            
        omega : float
            Energy of the vibrational quantum / frequency of the oscillator
            
            
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mod = mol.get_Mode(0)
        >>> mod.set_energy(1, 1.3)
        >>> mod.get_energy(1)
        1.3
    
        >>> with qr.energy_units("1/cm"):
        ...     mod.set_energy(1, 1300.0)
        >>> mod.get_energy(1)
        0.2448747037501509
        
        """
        sbm = self.submodes[N]
        sbm.omega = self.convert_energy_2_internal_u(omega)  

          
    def set_shift(self, N, shift):
        """Sets the potential energy surface shift with respect to ground state
        
        The shift is dimensionless. frequency*(shift^2)/2 is the reorganization
        energy.
        
        Parameters
        ----------
        
        N : int
            Index of the electronic state for we set frequency
            
        shift : float
            Shift of the PES with respect to ground state
            
            
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mod = mol.get_Mode(0)
        >>> mod.set_shift(1, 0.1)
        >>> mod.get_shift(1)
        0.1
    
        
        """
        sbm = self.submodes[N]
        sbm.shift = shift

        
    def set_nmax(self, N, nmax):
        """Sets maximum quantum number of the mode in an electronic state N
        
        
        Parameters
        ----------
        
        N : int
            Index of the electronic state
            
        nmax : int
            Maximum quantum number to be set
            
        
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mod = mol.get_Mode(0)
        
        Hamiltonian of a two-level molecule with one vibrational mode, for 
        which we take 5 levels in each electronic state, has 10 levels
        
        >>> mod.set_nmax(0, 5)
        >>> mod.set_nmax(1, 5)
        >>> H = mol.get_Hamiltonian()
        >>> print(H.dim)
        10
        
        """
        sbm = self.submodes[N]
        sbm.nmax = nmax

        
    def set_HR(self, N, hr):
        """Sets Huang-Rhys factor of the PES on the Nth electronic state
        
        
        Parameters
        ----------
        
        N : int
            Index of the state for which we set. Setting for N=0 gives
            exception
            
        hr : float
            Huang-Rhys factor. Dimensionless quantity (ratio of reorganization
            energy and vibrational quantum)
        
        
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mod = mol.get_Mode(0)
        >>> mod.set_HR(1, 1.0)
        
        HR is the shift square devide by 2
        
        >>> freq = mod.get_energy(1)
        >>> sft = mod.get_shift(1)
        >>> hr = (sft**2)/2.0
        >>> print("{0:.2f}".format(hr))
        1.00

        Exception will be thrown if HR factor is set to the ground-state
        
        >>> mod.set_HR(0, 2.0)
        Traceback (most recent call last):
            ...
        Exception: Makes no sense to set HR in the ground state
        
        """
        if N==0:
            raise Exception("Makes no sense to set HR in the ground state")
            
        sh = numpy.sqrt(2.0*hr)
        self.set_shift(N, sh)
        
        
        
    #
    # Getters
    #
    
    @deprecated
    def get_frequency(self, N):
        """Returns vibrational frequency
        
        Usage of this method is deprecated, use `get_energy` instead
        
        Parameters
        ----------
        
        N : int
            Index of the electronic state
            
            
        
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mod = mol.get_Mode(0)
        >>> mod.set_energy(1, 1.3) 
        
        This function is deprecated and it throws a warning text when it is run
                
        >>> mod.get_frequency(1) # doctest: +ELLIPSIS
        function  <function Mode.get_frequency at ...>  is deprecated
        1.3
        
        """
        return self.get_energy(N)


    def get_energy(self, N, no_conversion=True):
        """Returns frequency of the mode corresponding to Nth electronic state
        
        
        Parameters
        ----------
        
        N : int
            Index of the electronic state
            
        no_conversion : bool
            If set to True, the function does not react to units management
            
            
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mod = mol.get_Mode(0)
        >>> mod.get_energy(1)
        1.0

        >>> with qr.energy_units("1/cm"):
        ...     mod.get_energy(1)  # default is `no_conversion=True`
        1.0
        
        >>> with qr.energy_units("1/cm"):
        ...     mod.get_energy(1, no_conversion=False)  # default is `no_conversion=True`        
        5308.837458876145
        
        
        """
        
        if no_conversion:
            return self.submodes[N].omega
        else:
            return self.convert_energy_2_current_u(self.submodes[N].omega)

    
    def get_shift(self, N):
        """Returns the shift of the PES of the mode at Nth electronic state
        
        
        Parameters
        ----------
        
        N : int
            Index of the electronic state
            
            
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mod = mol.get_Mode(0)
        >>> mod.get_shift(0)
        0.0
        
        >>> mod.get_shift(1)
        0.0

        >>> mod.set_shift(1, 1.0)
        >>> mod.get_shift(1)
        1.0
        
        """
        return self.submodes[N].shift

    
    def get_nmax(self, N):
        """Returns maximum quantum number of the mode at electronic state N
        
        Parameters
        ----------
        
        N : int
            Index of the electronic state
            
            
        Examples
        --------

        Default value of `nmax` is 2
        
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mod = mol.get_Mode(0)
        >>> mod.get_nmax(0)
        2
        >>> mod.get_nmax(1)
        2
        
        
        """
        return self.submodes[N].nmax

        
    def get_HR(self, N):
        """Returns Huang-Rhys factor of vib. mode in Nth electronic state
        
        
        Parameters
        ----------
        
        N : int
            Index of the electronic state
            
            
        Examples
        --------
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mod = mol.get_Mode(0)
        >>> mod.get_HR(0)
        0.0

        >>> mod.set_HR(1, 1.0)
        >>> print("{0:.2f}".format(mod.get_HR(1)))
        1.00
        
        >>> mod.get_HR(0)
        0.0
        
        """
        return (self.submodes[N].shift**2)/2.0

    
    def set_all(self, N, param):
        """Sets all the parameters of an intramolecular mode at once
        
        
        Parameters
        ----------
        
        N : int
            Index of the electronic states for which we set the mode parameters
            
        param : list like
            List of the mode parameters in the following order corresponding
            to the properties of the SubMode class, i.e. `omega`, `shift`,
            `nmax`. The parameters `omega` is units managed.
        
        
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mod = mol.get_Mode(0)
        >>> mod.set_all(1, [1.2, 0.8, 3])
            
        
        >>> print("{0:1.2f}".format(mod.get_HR(1)))
        0.32
        
        >>> mod.get_shift(1)
        0.8
        
        >>> mod.get_energy(1)
        1.2
        
        >>> mod.get_nmax(1)
        3
        
        Energy is units management sensitive
        
        >>> with qr.energy_units("1/cm"):
        ...     mod.set_all(1, [1200, 0.8, 3])
        >>> mod.get_energy(1)
        0.22603818807706239
        
        """
        sbm = self.submodes[N]
        sbm.omega = self.convert_energy_2_internal_u(param[0])
        sbm.shift = param[1]
        sbm.nmax  = param[2]

        
    def get_SubMode(self, N):
        """Returns the SubMode class of a give electronic state of the molecule
        
        
        Parameters
        ----------
        
        N : int
            Index of the electronic state
            
            
        Examples
        --------

        >>> import quantarhei as qr
        >>> mol = qr.TestMolecule("two-levels-1-mode")
        >>> mod = mol.get_Mode(0)
        >>> sm = mod.get_SubMode(1)

        >>> print(sm.omega, sm.shift, sm.nmax)
        1.0 0.0 2
        
        
        """
        return self.submodes[N]
