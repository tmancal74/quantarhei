# -*- coding: utf-8 -*-
"""
    Intramolecular vibrations are represented by **two** classes:
       
    **Mode** class is a user class through which the user defines an intra
    molecular mode
    
    **SubMode** class is an internal class thorough which the mode keeps
    track of the fact that in different electronic states of the molecule,
    a given vibrational mode has different parameters.
    
    
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

from ..utils import Float
from ..utils import Integer
from ..utils import Bool

from ..core.managers import UnitsManaged, energy_units
from ..core.wrappers import deprecated

from ..core.saveable import Saveable

class SubMode(UnitsManaged, Saveable):
    """ Instance of a vibrational mode relative to a give electronic state 
    
    When a mode is set on a Molecule object, it has to be indepedently 
    set on each electronic state of the molecule. We keep track of individual
    parameters of the vibrational mode for each electronic state using
    this class. SubMode itself has no idea about this.
    
    
    Examples
    --------
    
    >>> sm = SubMode()
    >>> print(sm.nmax)
    2
    >>> print(sm.shift)
    0.0
    >>> print(sm.omega)
    1.0
    
    This class is aware of energy units
    >>> import quantarhei as qr
    >>> with energy_units("1/cm"):
    ...     sm = SubMode()
    >>> print(sm.omega)
    0.0001883651567308853
    
    """
    
    omega = Float('omega')
    shift = Float('shift')
    nmax  = Integer('nmax')
    
    def __init__(self, omega=1.0, shift=0.0, nmax=2):
        self.omega = self.convert_energy_2_internal_u(omega)
        self.shift = shift
        self.nmax  = nmax
        
        
        
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
        try:
            self.nel = monomer.nel
            self.monomer = monomer
            
            # Must have as many submodes as electronic states in the monomer
            with energy_units("int"):
                for k in range(1, self.nel):
                    # submodes are created with the same frequency as
                    # in the groundstate and with zero shift
                    self.submodes.append(SubMode(self.submodes[0].omega))
        except:
            raise

        self.monomer_set = True
        
    """ Setters """
    @deprecated
    def set_frequency(self,N,omega):
        #sbm = self.submodes[N]
        #sbm.omega = omega
        self.set_energy(N, omega)
        
    def set_energy(self, N, omega):
        sbm = self.submodes[N]
        sbm.omega = self.convert_energy_2_internal_u(omega)  
          
    def set_shift(self, N, shift):
        sbm = self.submodes[N]
        sbm.shift = shift
        
    def set_nmax(self, N, nmax):
        sbm = self.submodes[N]
        sbm.nmax = nmax
        
    def set_HR(self, N, hr):
        """Sets Huang-Rhys factor
        
        """
        if N==0:
            raise Exception("Makes no sense to set HR in the ground state")
            
        sh = numpy.sqrt(2.0*hr)
        self.set_shift(N, sh)
        
        
        
    """ Getters """
    @deprecated
    def get_frequency(self,N):
        #return self.submodes[N].omega
        return self.get_energy(N)

    def get_energy(self, N, no_conversion=True):
        if no_conversion:
            return self.submodes[N].omega
        else:
            return self.convert_energy_2_current_u(self.submodes[N].omega)

    
    def get_shift(self, N):
        return self.submodes[N].shift
    
    def get_nmax(self, N):
        return self.submodes[N].nmax
        
    def get_HR(self, N):
        return (self.submodes[N].shift**2)/2.0
    
    def set_all(self, N, param):
        sbm = self.submodes[N]
        sbm.omega = self.convert_energy_2_internal_u(param[0])
        sbm.shift = param[1]
        sbm.nmax  = param[2]
        
    def get_SubMode(self, N):
        return self.submodes[N]
