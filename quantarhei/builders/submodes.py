# -*- coding: utf-8 -*-
"""
    The SubMode class is an internal class through which the mode keeps
    track of the fact that in different electronic states of the molecule,
    a given vibrational mode has different parameters.
    
    
"""

from ..utils import Float
from ..utils import Integer

from ..core.managers import UnitsManaged

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
    >>> with qr.energy_units("1/cm"):
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
        
