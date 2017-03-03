# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt

from ..core.time import TimeAxis
from ..builders.aggregates import Aggregate
from ..builders.molecules import Molecule

from ..utils import derived_type

class DFunction2:
    """Descrete two-dimensional function 
    
    """
    def __init__(self, x=None, y=None, z=None):
        pass
    
    def save(self, filename):
        pass
    
    def load(self, filename):
        pass
    
    def plot(self,**kwargs):
        pass    


class TwoDSpectrumBase(DFunction2):
    """Basic container for two-dimensional spectrum
    
    
    """
    
    def __init__(self):
        super().__init__()
        self.data = None
        self.xaxis = None
        self.yaxis = None



    def set_axis_1(self, axis):
        self.xaxis = axis
        
    def set_axis_3(self, axis):
        self.yaxis = axis
        
    def set_data(self, data):
        self.data = data
        
    def save(self, filename):
        super().save(filename)
    
    def load(self, filename):
        super().load(filename)
    
    def plot(self,**kwargs):
        pass
    
    
class TwoDSpectrumContainer(TwoDSpectrumBase):
    
    def __init__(self):
        pass
    
    
    
class TwoDSpectrumCalculator:
    """Calculator of the 2D spectrum
    
    
    Enables setting up parameters of 2D spectrum calculation for later
    evaluation. The method `calculate` returns TwoDSpectrumContainer
    with a 2D spectrum.
    
    Parameters
    ----------
    
    
    """

    t1axis = derived_type("t1axis",TimeAxis)
    t2axis = derived_type("t2axis",TimeAxis)
    t3axis = derived_type("t3axis",TimeAxis)
    
    system = derived_type("system",[Molecule,Aggregate])
    
    def __init__(self, t1axis, t2axis, t3axis,
                 system=None,
                 dynamics="secular",
                 relaxation_tensor=None,
                 rate_matrix=None,
                 effective_hamiltonian=None):
        self.t1axis = t1axis
        self.t2axis = t2axis
        self.t3axis = t3axis
        
        #FIXME: properties to be protected
        self.dynamics = dynamics
        
        # unprotected properties
        self.data = None
        
        self._relaxation_tensor = None
        self._rate_matrix = None
        self._relaxation_hamiltonian = None
        self._has_relaxation_tensor = False
        if relaxation_tensor is not None:
            self._relaxation_tensor = relaxation_tensor
            self._has_relaxation_tensor = True
        if effective_hamiltonian is not None:
            self._relaxation_hamiltonian = effective_hamiltonian
        if rate_matrix is not None:
            self._rate_matrix = rate_matrix
            self._has_rate_matrix = True
            
        self._have_aceto = False
       
        
    def calculate(self):
        """Returns 2D spectrum
        
        Calculates and returns TwoDSpectrumContainer containing 2D spectrum
        based on the parameters specified in this object.
        
        
        """

        try:
            import aceto.nr3td as nr3td 
            
            self._have_aceto = True
            
        except:
            #
            # FIXME: There should be an optional warning and a fall back onto
            # quantarhei.implementations.aceto module
            #
            raise Exception("Aceto not available")
            
            from ..implementations.aceto import nr3td
            
            self._have_aceto = False
    
    
        if self._have_aceto:
            
            # calculate 2D spectrum using aceto library

        
            ret = TwoDSpectrumContainer()
            
        else:
            
            # fall bakc on quantarhei's own implementation
        
            ret = TwoDSpectrumContainer()
            
        
        return ret
    
    
    
        
        