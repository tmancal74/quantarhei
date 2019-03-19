# -*- coding: utf-8 -*-
"""
    Class to support tests on Aggregate class

    This class provides an easy access to a prefilled, but not initiallized
    object of the Aggregate class. These can be used for testing and 
    demonstration. There are several types of aggregates available, 
    distinguished by the `name` argument of the constructor.
        
        
    TestAggregates Provided 
    -----------------------
    
    Below we list the valid values for the `name` argument of the TestAggregate
    constructor:
                
    dimer-2 :
        Dimer of two-level molecules, with positions in space
        and transition dipole moments specified. No environment
        is defined.
        
    dimer-2-env :
        Dimer of two-level molecules, with positions in space and
        transition dipole moments specified. For each molecule
        we define energy gap correlation function (energy gao 
        correlation functions on different sites are not correlated). 

    homodimer-2 :
        Homo-dimer of two-level molecules (molecules with the same energy 
        gaps), with positions in space
        and transition dipole moments specified. No environment
        is defined.
        
    homodimer-2-env :
        Homo-dimer of two-level molecules, with positions in space and
        transition dipole moments specified. For each molecule
        we define energy gap correlation function (energy gao 
        correlation functions on different sites are not correlated).
        
    Class Details
    -------------

"""

import numpy

from .aggregates import Aggregate
from .molecules import Molecule
from .modes import Mode
from ..core.units import convert
from ..core.time import TimeAxis
from ..qm.corfunctions.correlationfunctions import CorrelationFunction
from ..core.managers import energy_units 

class TestAggregate(Aggregate):
    """Class to support tests on Aggregate class
    
    
    Parameters
    ----------
    
    name : str
        Name characterizing the test aggregate.
        
        
    Examples
    --------

    General dimers
    
    >>> # Dimer of two-level systems 
    >>> tagg = TestAggregate(name="dimer-2")
    >>> tagg.build()
    >>> tagg.has_SystemBathInteraction()
    False

    
    >>> # Dimer of two-level systems with an environment
    >>> tagg = TestAggregate(name="dimer-2-env")
    >>> tagg.build()
    >>> tagg.has_SystemBathInteraction()
    True
    
    Homo-dimers
    
    >>> # Dimer of two-level systems 
    >>> tagg = TestAggregate(name="homodimer-2")
    >>> tagg.build()
    >>> tagg.has_SystemBathInteraction()
    False

    
    >>> # Dimer of two-level systems with an environment
    >>> tagg = TestAggregate(name="homodimer-2-env")
    >>> tagg.build()
    >>> tagg.has_SystemBathInteraction()
    True

    >>> # Trimer of two-level systems without an environment
    >>> tagg = TestAggregate(name="trimer-2")
    >>> tagg.build()
    >>> tagg.has_SystemBathInteraction()
    False
    
    """
    
    def __init__(self, name=None):
        """ Some more doctests
        
        >>> TestAggregate()
        Traceback (most recent call last):
            ...
        Exception: Aggregate name not specified
        
        
        """
        
        if name is None:
            raise Exception("Aggregate name not specified")
            
            
        #
        # Test dimer
        #
        if name == "dimer-2-env":
            
            m1, m2 = self._molecules(N=2, nst=2)
            
            # set their environment
            time = TimeAxis(0, 1000, 1.0)
            cpar = dict(ftype="OverdampedBrownian", reorg=20,
                        cortime=100, T=300)
            with energy_units("1/cm"):
                cfce = CorrelationFunction(time, cpar)
                
            m1.set_transition_environment((0, 1), cfce)
            m2.set_transition_environment((0, 1), cfce)
                
            super().__init__(molecules=[m1, m2])

        elif name == "trimer-2-env":
            
            m1, m2, m3 = self._molecules(N=3, nst=2)
            
            # set their environment
            time = TimeAxis(0, 1000, 1.0)
            cpar = dict(ftype="OverdampedBrownian", reorg=20,
                        cortime=100, T=300)
            with energy_units("1/cm"):
                cfce = CorrelationFunction(time, cpar)
                
            m1.set_transition_environment((0, 1), cfce)
            m2.set_transition_environment((0, 1), cfce)
            m3.set_transition_environment((0, 1), cfce)
                
            super().__init__(molecules=[m1, m2, m3])
            
        elif name == "dimer-2":
            
            m1, m2 = self._molecules(N=2, nst=2)
           
            # set their environment
            # nothing here
                
            super().__init__(molecules=[m1, m2])
            
        elif name == "trimer-2":
            
            m1, m2, m3 = self._molecules(N=3, nst=2, homo=False)
            
            super().__init__(molecules=[m1, m2, m3])
            
        elif name == "homodimer-2-env":
            
            m1, m2 = self._molecules(N=2, nst=2, homo=True)
            
            # set their environment
            time = TimeAxis(0, 1000, 1.0)
            cpar = dict(ftype="OverdampedBrownian", reorg=20,
                        cortime=100, T=300)
            with energy_units("1/cm"):
                cfce = CorrelationFunction(time, cpar)
                
            m1.set_transition_environment((0, 1), cfce)
            m2.set_transition_environment((0, 1), cfce)
                
            super().__init__(molecules=[m1, m2]) 
            
        elif name == "homodimer-2":
            
            m1, m2 = self._molecules(N=2, nst=2, homo=True)
           
            # set their environment
            # nothing here
                
            super().__init__(molecules=[m1, m2])   
            
        elif name == "dimer-2-vib":
            
            m1, m2 = self._molecules(N=2, nst=2)
           
            with energy_units("1/cm"):
                mod1 = Mode(100.0)
                m1.add_Mode(mod1)
                mod1.set_HR(1,0.1)
                
                mod2 = Mode(100.0)
                m2.add_Mode(mod2)
                mod2.set_HR(1,0.1)
                
            super().__init__(molecules=[m1, m2])            
  
    
    def _molecules(self, N, nst, homo=False):
        """Creates molecules to be filled into Aggregate
        
        Testing that None is returned for wrong arguments
        
        >>> tagg = TestAggregate("dimer-2")
        >>> mols = tagg._molecules(3, 5)
        >>> print(mols)
        None
        
        """
        
        if (N == 2) and (nst == 2):
            
            nstates = nst
            
            # check inputs
            if nstates != 2:
                raise Exception()
                
            # set parameters
            gap1 = convert(12000,"1/cm", to="int")
            
            energies1 = numpy.zeros(nstates)
            for s in range(nstates):
                energies1[s] = s*gap1
                
            if homo:
                gap2 = convert(12000,"1/cm", to="int")
            else:
                gap2 = convert(12300,"1/cm", to="int")
                
            energies2 = numpy.zeros(nstates)
            for s in range(nstates):
                energies2[s] = s*gap2
                
            # molecules
            m1 = Molecule(elenergies=energies1)
            m2 = Molecule(elenergies=energies2)
            
            # set transition dipole moments
            dip1 = [0.0, 2.0, 0.0]
            dip2 = [0.0, 1.3, 1.4]
            m1.set_dipole(0, 1, dip1)
            m2.set_dipole(0, 1, dip2)
            
            #set molecular positions
            r1 = [0.0, 0.0, 0.0]
            r2 = [5.0, 0.0, 0.0]
            m1.position = r1
            m2.position = r2
 
            return [m1, m2]
        
        elif (N == 3) and (nst == 2):
        
            nstates = nst
            
            # check inputs
            if nstates != 2:
                raise Exception()
                
            # set parameters
            gap1 = convert(12000, "1/cm", to="int")
            energies1 = numpy.zeros(nstates)
            for s in range(nstates):
                energies1[s] = s*gap1
                
            if homo:
                gap2 = convert(12000, "1/cm", to="int")
                gap3 = gap2
            else:
                gap2 = convert(12300, "1/cm", to="int")
                gap3 = convert(12350, "1/cm", to="int")
            energies2 = numpy.zeros(nstates)
            energies3 = numpy.zeros(nstates)
            for s in range(nstates):
                energies2[s] = s*gap2
                energies3[s] = s*gap3
                
            # molecules
            m1 = Molecule(elenergies=energies1)
            m2 = Molecule(elenergies=energies2)
            m3 = Molecule(elenergies=energies3)
            
            # set transition dipole moments
            dip1 = [0.0, 2.0, 0.0]
            dip2 = [0.0, 1.3, 1.4]
            dip3 = [1.0, 1.2, 0.0]
            m1.set_dipole(0, 1, dip1)
            m2.set_dipole(0, 1, dip2)
            m3.set_dipole(0, 1, dip3)
            
            #set molecular positions
            r1 = [0.0, 0.0, 0.0]
            r2 = [5.0, 0.0, 0.0]
            r3 = [0.0, 0.0, 5.0]
            m1.position = r1
            m2.position = r2
            m3.position = r3
 
            return [m1, m2, m3]
                        
        else:
            
            return None
        