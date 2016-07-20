# -*- coding: utf-8 -*-

from .valueaxis import ValueAxis
from .frequency import FrequencyAxis

import numpy

class TimeAxis:
    """ Class representing time in time dependent calculations.
    

          
    Parameters
    ----------
    startValue : float
        start of the TimeAxis
        
    noPoints : int
        number of steps

    step : float
        time step
        
    Attributes
    ----------
    data : float array 
        Holds the values of time, it is equivalent to the atribute
        TimeAxis.time
        
    time : float array
        equivalent to the attribute TimeAxis.values
        
        
    Examples
    --------
        
    >>> ta = TimeAxis(0.0,100,0.1)
        
        
    """
    
    def __init__(self,startValue,noPoints,step):
    
        ValueAxis.__init__(self,startValue=startValue,
                           noPoints=noPoints,step=step)
        self.time = self.data
        self.dt = self.step
        self.length = self.noPoints
        self.tmin = self.startValue
        self.tmax = self.endValue        
        
        
    def get_FrequencyAxis(self):
        """ Returns corresponding FrequencyAxis object 
        
        
        """
        pass
        
        
        
        
        