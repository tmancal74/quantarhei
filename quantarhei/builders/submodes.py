# -*- coding: utf-8 -*-
"""
    The SubMode class is an internal class through which the mode keeps
    track of the fact that in different electronic states of the molecule,
    a given vibrational mode has different parameters.
    
    
"""
import numpy

from ..utils import Float
from ..utils import Integer

from ..core.managers import UnitsManaged

from ..core.saveable import Saveable

#from .opensystem import OpenSystem
from ..qm.hilbertspace.hamiltonian import Hamiltonian 
from .. import REAL

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
        


class HarmonicMode(SubMode): #, OpenSystem):
    """Renaming the SubMode to be used as a standalone Harmonic oscillator mode"""

    def __init__(self, omega=1.0, shift=0.0, nmax=2):
        super().__init__(omega=omega, shift=shift, nmax=nmax)
        self._built = False


    def build(self, nmax=None):
        """Building all necessary quantities """

        # if provided, we reset nmax
        if nmax is not None:
            self.nmax = nmax

        N = self.nmax
        HH = numpy.zeros((N,N), dtype=REAL)

        for nn in range(N):
            HH[nn,nn] = self.omega*(nn + 0.5)

        HH -= HH[0,0]
        
        self.HH = HH

        self.HamOp = Hamiltonian(data=HH)

        self._built = True


    def get_Hamiltonian(self):
        """Returns the system Hamiltonian 
        
        """
        if self._built:
            return self.HamOp
        else:
            raise Exception("The Mode has to be built first.")
        



class AnharmonicMode(HarmonicMode):
    """
    
    
    
    """
    def __init__(self, omega=1.0, shift=0.0, nmax=2):
        """
        
        In this case, omega is the difference between levels 1 and 0

        """
        super().__init__(omega=omega, shift=shift, nmax=nmax)

        # unharmonicity
        self.xi = 0.0
        self.om0 = self.omega*(1.0/(1.0-2.0*self.xi))
        self.dom = self.om0 - self.omega


    def set_anharmonicity(self, xi):
        """Sets the ocillator anharmonicity
        
        """
        self.xi = xi
        # corresponding harmonic frequency
        self.om0 = self.omega*(1.0/(1.0-2.0*self.xi))
        self.dom = self.om0 - self.omega

        
    
    def get_anharmonicity(self, dom=None):
        """Returns the ahnarmonicity set previously, or calculates anharmonicity from frequency difference
        
        """
        if dom is None:
            return self.xi 
        else:

            dom_int = self.convert_energy_2_internal_u(dom)

            return dom_int/(2.0*(self.omega + dom_int))
        

         


    def build(self, nmax=None, xi=None):
        """Building all necessary quantities """

        # if provided, we reset nmax
        if nmax is not None:
            self.nmax = nmax

        if xi is not None:
            self.set_unharmonicity(xi)

        N = self.nmax
        HH = numpy.zeros((N,N), dtype=REAL)

        # energy levels are known exactly
        for nn in range(N):
            HH[nn,nn] = self.om0*(nn + 0.5) - self.om0*self.xi*((nn + 0.5)**2)
            if nn == 0:
                hzer = HH[0,0]
            HH[nn,nn] -= hzer
        
        self.HH = HH

        self.HamOp = Hamiltonian(data=HH)

        self._built = True