# -*- coding: utf-8 -*-

from .valueaxis import ValueAxis
#from .time import TimeAxis

class FrequencyAxis:
    """ Class representing frequency axis of calculations
    

          
    Parameters
    ----------
    startValue : float
        start of the FrequencyAxis
        
    noPoints : int
        number of steps

    step : float
        time step
        
    Attributes
    ----------
    data : float array 
        Holds the values of time, it is equivalent to the atribute
        TimeAxis.time
        
    frequency : float array
        equivalent to the attribute TimeAxis.values
        
        
    Examples
    --------
        
    >>> ta = FrequencyAxis(0.0,100,0.1)
        
        
    """
    
    def __init__(self,startValue,noPoints,step):
    
        ValueAxis.__init__(self,startValue=startValue,
                           noPoints=noPoints,step=step)
        self.frequency = self.data
        self.domega = self.step
        self.length = self.noPoints
        self.omin = self.startValue
        self.omax = self.endValue    
        
    def get_TimeAxis(self):
        pass
    
    
    