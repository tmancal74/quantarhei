# -*- coding: utf-8 -*-

import numpy
import scipy.interpolate

from ...core.dfunction import DFunction
from ...core.units import kB_intK

from ...core.managers import UnitsManaged
from ...core.managers import energy_units
from ...core.frequency import FrequencyAxis 
from ...core.time import TimeAxis
from .correlationfunctions import CorrelationFunction

class SpectralDensity(DFunction, UnitsManaged):
    
    
    def __init__(self,axis, params):
        
        #FIXME: valueAxis attribute should be changed to axis
        #FIXME: attribute frequencyAxis should be deleted
        if isinstance(axis,TimeAxis):
            # protect the frequency axis creation from units management
            with energy_units("int"):
                faxis = axis.get_FrequencyAxis()
            self.valueAxis = faxis
        else:
            self.valueAxis = axis
            
        self.frequencyAxis = self.valueAxis
         
        self._splines_initialized = False 
        
        try:
            ftype = params["ftype"]
            if ftype in CorrelationFunction.allowed_types:
                self.ftype = ftype
            else:
                raise Exception("Unknown Correlation Function Type")
                 
            # we need to save the defining energy units
            self.energy_units = self.manager.get_current_units("energy")
            # because params are in the defining units
            self.params = params

        except:
            raise Exception
            
            
        if (self.ftype == "OverdampedBrownian"):
            
            T    = params["T"]
            tc   = params["cortime"]
            lamb = self.manager.iu_energy(params["reorg"],
                                          units=self.energy_units)    
            
            # Masubara frequencies are needed only in time domain
            #try:
            #    N = params["matsubara"]
            #except:
            #    N = 10
            
            # Temperature is also not needed here, but it is stored
            #kBT = kB_intK*T            

            # protect calculation from units management
            with energy_units("int"):     
                w = self.frequencyAxis.data   
                cc = (2.0*lamb/tc)*w/(w**2 + (1.0/tc)**2)
            
            #self.data = cc
            self._make_me(self.frequencyAxis,cc)
            
            # this is in internal units
            self.lamb = lamb
            self.T    = T
            
            
        else:
            raise Exception("Unknown correlation function type of"+
                            "type domain combination.")
                            

    def get_temperature(self):
        """Returns the temperature of the correlation function
        
        """
        return self.T

        
    def copy(self):
        """Creates a copy of the current correlation function"""
        return SpectralDensity(self.frequencyAxis,self.params)                            
        
        
    def get_CorrelationFunction(self):
        
        ta = self.frequencyAxis.get_TimeAxis()        
        return CorrelationFunction(ta,self.params)
        