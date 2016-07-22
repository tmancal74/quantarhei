# -*- coding: utf-8 -*-

from ..utils import derived_type
from ..builders import Molecule 
from ..core.time import TimeAxis

class AbsSpect:
    """Linear absorption spectrum 
    
    Linear absorption spectrum of a molecule or an aggregate of molecules.
    
    Parameters
    ----------
    timeaxis : quantarhei.TimeAxis
        TimeAxis on which the calculation will be performed
        
    system : quantarhei.Molecule or quantathei.Aggregate
        System for which the absorption spectrum will be calculated.
        
        
    """

    TimeAxis = derived_type("TimeAxis",TimeAxis)
    system = derived_type("system",[Molecule])
    
    def __init__(self,timeaxis,system=None,dynamics="secular"):
        # protected properties
        self.TimeAxis = timeaxis
        self.system = system
        
        #FIXME: properties to be protected
        self.dynamics = dynamics
        
        # unprotected properties
        self.data = None
        
    def build(self,rwa=0.0):
        """ Calculates the absorption spectrum 
        
        
        """
        if self.system is not None:
            if isinstance(self.system,Molecule):
                self._calculate_Molecule(rwa)                
#            elif isinstance(self.system,aggregate):
#                self._calculate_aggregate(rwa)
        else:
            raise Exception("System to calculate spectrum for is not defined")
        


    def _calculateMolecule(self,rwa):
        
        if self.system._has_system_bath_coupling:
            raise Exception("Not yet implemented")
        else: 
            # calculating stick spectra  

            stick_width = 1.0/0.1
            
            
        