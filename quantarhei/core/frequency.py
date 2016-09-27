# -*- coding: utf-8 -*-

from .valueaxis import ValueAxis
from .managers import EnergyUnitsManaged
from .managers import energy_units
from ..utils.types import UnitsManagedReal

import numpy

class FrequencyAxis(ValueAxis,EnergyUnitsManaged):
    """ Class representing frequency axis of calculations
    

          
    Parameters
    ----------
    startValue : float
        start of the FrequencyAxis
        
    noPoints : int
        number of steps

    step : float
        time step
        
    atype : string {"complete","upper-half"}
        Axis type
        
    Attributes
    ----------
    data : float array 
        Holds the values of time, it is equivalent to the atribute
        TimeAxis.time
        
    frequency : float array
        equivalent to the attribute TimeAxis.values
        
        
    Examples
    --------
        
    The default type of the `FrequencyAxis` is `complete`. See the discussion
    of the types in the `TimeAxis` documentation. 
    
    >>> wa = FrequencyAxis(0.0,100,0.05)
    >>> ta = wa.get_TimeAxis()
    >>> print(ta.length)
    100
    
    The type `upper-half` only refers to the corresponding TimeAxix. Everything
    about the FrequencyAxis remains the same as with `complete`.
    
    >>> wa = FrequencyAxis(0.0, 100, 0.05, atype = "upper-half")
    >>> ta = wa.get_TimeAxis()
    >>> print(ta.length)
    50
    
    #>>> print(ta.step,2.0*numpy.pi/(100*wa.domega))
    >>> print(numpy.allclose(ta.step,2.0*numpy.pi/(100*wa.domega)))
    True
    
    For `complete`, everything should work also for an odd number of points
    
    >>> wa = FrequencyAxis(0.0,99,0.05)
    >>> ta = wa.get_TimeAxis()
    >>> print(ta.length)
    99
    
    But `upper-half` throughs an exception, because by definition its number of
    points is `2*N`, where `N` is an integer.
    
    >>> wa = FrequencyAxis(0.0,99,0.05,atype="upper-half")
    >>> ta = wa.get_TimeAxis()
    Traceback (most recent call last):
    ...
    Exception: Cannot create upper-half TimeAxis from an odd number of points
    


    Relation between TimeAxis and FrequencyAxis
    -------------------------------------------
    
    Complete FrequencyAxis and even number of points
    
    >>> wa = FrequencyAxis(0.0,10,0.1,atype="complete")
    >>> frequencies = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    >>> print(numpy.allclose(wa.frequency,frequencies))
    True
    
    >>> ta = wa.get_TimeAxis()
    >>> times = 2.0*numpy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(10,0.1))
    >>> print(numpy.allclose(ta.time,times))
    True
        
    >>> print(numpy.allclose(ta.dt,times[1]-times[0]))
    True
    
    >>> wb = ta.get_FrequencyAxis()
    >>> print(numpy.allclose(wb.frequency,frequencies))
    True
    
    >>> tb = wb.get_TimeAxis()
    >>> print(numpy.allclose(tb.time,times))
    True
    
    Complete FrequencyAxis and odd number of points
    
    >>> wa = FrequencyAxis(0.0,11,0.1,atype="complete")
    >>> frequencies = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    >>> print(numpy.allclose(wa.frequency,frequencies))
    True
    
    >>> ta = wa.get_TimeAxis()
    >>> times = 2.0*numpy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(11,0.1))
    >>> print(numpy.allclose(ta.time,times))
    True
        
    >>> print(numpy.allclose(ta.dt,times[1]-times[0]))
    True
    
    >>> wb = ta.get_FrequencyAxis()
    >>> print(numpy.allclose(wb.frequency,frequencies))
    True
    
    >>> tb = wb.get_TimeAxis()
    >>> print(numpy.allclose(tb.time,times))
    True
        
    Upper-half FrequencyAxis and even number of points
    
    >>> wa = FrequencyAxis(0.0,10,0.1,atype="upper-half")
    >>> frequencies = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    >>> print(numpy.allclose(wa.frequency,frequencies))
    True
    
    >>> ta = wa.get_TimeAxis()
    >>> times = 2.0*numpy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(10,0.1))
    >>> print(numpy.allclose(ta.time,times[5:10]))
    True
        
    >>> print(numpy.allclose(ta.dt,times[1]-times[0]))
    True
    
    >>> wb = ta.get_FrequencyAxis()
    >>> print(numpy.allclose(wb.frequency,frequencies))
    True
    
    >>> tb = wb.get_TimeAxis()
    >>> print(numpy.allclose(tb.time,times[5:10]))
    True
        
    """
    
    data = UnitsManagedReal("data")
    
    def __init__(self,startValue,noPoints,step,atype='complete',t0=0.0):
   
        ValueAxis.__init__(self,startValue=startValue,
                           noPoints=noPoints,step=step)
                           
        self.frequency = self.data
        
        # these have to be converted manually into internal units
        self.domega = self.convert_2_internal_u(self.step)
        self.length = self.noPoints
        self.omin = self.convert_2_internal_u(self.startValue)
        self.omax = self.convert_2_internal_u(self.endValue)    
        self.t0 = t0

        self.allowed_atypes = ["upper-half","complete"]
        if atype in self.allowed_atypes:         
            self.atype = atype
        else:
            raise Exception("Unknown frequency axis type")
            
    def get_TimeAxis(self):
        """Returns the corresponding TimeAxis object
        
        """
        from .time import TimeAxis
        
        with energy_units("int"):
        
            if self.atype == 'complete':
                tt = numpy.fft.fftshift(
                numpy.fft.fftfreq(self.length,self.domega/(2.0*numpy.pi)))
            
                step = tt[1]-tt[0]
                start = self.t0  + tt[0]
                nosteps = self.length
            
                w0 = self.frequency[self.length//2] 
            
            
            elif self.atype == 'upper-half':
            
                if not self.length % 2 == 0:
                    raise Exception("Cannot create upper-half TimeAxis"
                    + " from an odd number of points")
                
                tt = numpy.fft.fftshift(
                (2.0*numpy.pi)*numpy.fft.fftfreq(self.length,self.domega))
            
                start = tt[int(self.length/2)] + self.t0
                nosteps = int(self.length/2)
                step = tt[1]-tt[0]
            
                w0 = self.frequency[self.length//2]
        
            else:
                raise Exception("Unknown frequency axis type")

        
        return TimeAxis(start,nosteps,step,atype=self.atype,w0=w0)        
    
    
    