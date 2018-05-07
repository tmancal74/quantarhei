# -*- coding: utf-8 -*-

import numpy

from .aggregates import Aggregate
from .molecules import Molecule
from .. import convert
from ..core.time import TimeAxis
from ..qm.corfunctions.correlationfunctions import CorrelationFunction
from ..core.managers import energy_units 

class TestAggregate(Aggregate):
    
    def __init__(self, name=None):
        """Class to support tests on Aggregate class
        
        
        Parameters
        ----------
        
        name : str
            Name characterizing the test aggregate.
            
            
        Names
        -----
            
        dimer-2 :
            Dimer of two-level molecules, with positions in space
            and transition dipole moments specified. No environment
            is defined.
            
        dimer-2-env:
            Dimer of two-level molecules, with positions in space and
            transition dipole moments specified. For each molecule
            we define energy gap correlation function (energy gao 
            correlation functions on different sites are not correlated). 
            
        
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
            
        elif name == "dimer-2":
            
            m1, m2 = self._molecules(N=2, nst=2)
           
            # set their environment
            # nothing here
                
            super().__init__(molecules=[m1, m2])
  
    
    def _molecules(self, N, nst):
        
        if (N == 2) and (nst == 2):
            
            nstates = nst
            
            # check inputs
            if (nstates < 2) or (nstates > 3):
                raise Exception()
                
            # set parameters
            gap = convert(12000,"1/cm", to="int")
            
            energies = numpy.zeros(nstates)
            for s in range(nstates):
                energies[s] = s*gap
                
            # molecules
            m1 = Molecule(elenergies=energies)
            m2 = Molecule(elenergies=energies)
            
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
        
        else:
            
            return None