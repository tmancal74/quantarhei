# -*- coding: utf-8 -*-

import numpy

from .aggregates import Aggregate
from .molecules import Molecule
from .. import convert

class TestAggregate(Aggregate):
    
    def __init__(self, name=None, params=None):
        """Class to support tests on Aggregate class
        
        
        Parameters
        ----------
        
        name : str
            Name of the test aggregate
            
        params : dict
            Dictionary of TestAggregate parameters
            
        Names
        -----
            
        dimer {nstates=2/3, environment=True/False}
        """
        
        if name is None:
            raise Exception("Aggregate name not specified")
            
            
        #
        # Test dimer
        #
        if name == "dimer":
            
            nstates = params["nstates"]
            env     = params["environment"]
            
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
            
            # set their environment
            if env:
                pass
                
            super().__init__(molecules=[m1, m2])
            
            