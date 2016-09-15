# -*- coding: utf-8 -*-

from .valueaxis import ValueAxis

import numpy

class TimeAxis(ValueAxis):
    """ Class representing time in time dependent calculations.
    
    The `TimeAxis` class stands in a close relation to `FrequencyAxis`
    class of the `quantarhei` package. `FrequencyAxis` represents the
    frequencies one obtains in the Fast Fourier transform of a function of
    the `TimeAxis`. By default, `TimeAxis` is of the type `upper-half` which
    means that by specifying the `startValue`, `noPoints` and `step` we 
    represent the upper half of the interval `<startValue-noPoints*step,
    startValue+noPoints*step>`. The Fourier transform of a time dependent 
    object defined on the `TimeAxis` will then have twice as many points as 
    the `TimeAxis`. This is usefull when the time dependent object has some 
    special symmetries. One example is the so-called quantum bath correlation
    function which fulfills the relation (in LaTeX) 
    
    C(-t) = C^{*}(t)
    
          
    Parameters
    ----------
    startValue : float
        start of the TimeAxis
        
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
        
    time : float array
        equivalent to the attribute TimeAxis.values
        
        
    Examples
    --------
        
    Default TimeAxis is of the 'upper-half' type
    
    >>> ta = TimeAxis(0.0,100,0.1)
    >>> ta.atype
    'upper-half'
    
    It is defined between the `startValue` and `startvalue+(noPoints-1)*step`
    
    >>> '%.3f' % ta.tmin
    '0.000'
    
    >>> '%.3f' % ta.tmax
    '9.900'
 
    However, when we ask for the corresponding FrequencyAxis, we get an 
    array of frequencies which is twice as long and contains `2*noPoints` 
    points
    
    >>> wa = ta.get_FrequencyAxis()
    >>> print(wa.length)
    200

    The frequency step is therefore twice shorter than one would normally 
    expect.
    
    >>> print(wa.step == wa.domega)
    True
    >>> import numpy
    >>> print(numpy.allclose(wa.step,2.0*numpy.pi/(0.1*200)))
    True
    
    This definition of the `TimeAxis` is used in conjunction with 
    the symmetries of the corelation functions, lineshape functions
    and response functions.
    
    
    In some case you want the 'complete' type of the TimeAxis. 
    
    >>> ta = TimeAxis(0.0,100,1.0,atype="complete")
    >>> ta.atype
    'complete'
    
    In this case the relation of the `TimeAxis` and the `FrequencyAxis` is
    more straightforward. No the `TimeAxis` represents the interval on
    which a time dependent object is defined completely. The number of
    points in the corresponding `FrequencyAxis` is the same as in the
    `TimeAxis` and the frequency step is as expected from a normal FFT.
    
    >>> wa = ta.get_FrequencyAxis()
    >>> print(wa.length)
    100
    
    >>> print(wa.step == wa.domega)
    True
    
    >>> import numpy
    >>> print(numpy.allclose(wa.step,2.0*numpy.pi/(1.0*100)))
    True
    
    
    No other types than `complete` and `upper-half` are defined at the moment.
    The following with throw an Exception
    
    >>> ta = TimeAxis(0.0,100,1.0,atype="lower-half")
    Traceback (most recent call last):
    ...
    Exception: Unknown time axis type
    
    
    Relation between TimeAxis and FrequencyAxis
    -------------------------------------------
    
    Complete TimeAxis and even number of points
    
    >>> ta = TimeAxis(0.0,10,0.1,atype="complete")
    >>> times = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    >>> print(numpy.allclose(ta.time,times))
    True
    
    >>> wa = ta.get_FrequencyAxis()
    >>> freques = 2.0*numpy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(10,0.1))
    >>> print(numpy.allclose(wa.frequency,freques))
    True
        
    >>> print(numpy.allclose(wa.domega,freques[1]-freques[0]))
    True
    
    >>> tb = wa.get_TimeAxis()
    >>> print(numpy.allclose(tb.time,times))
    True
    
    >>> wb = tb.get_FrequencyAxis()
    >>> print(numpy.allclose(wb.frequency,freques))
    True
    
    Complete TimeAxis and odd number of points    

    >>> ta = TimeAxis(0.0,11,0.1,atype="complete")    
    >>> times = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    >>> print(numpy.allclose(ta.time,times))
    True
    
    >>> wa = ta.get_FrequencyAxis()
    >>> freques = 2.0*numpy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(11,0.1))
    >>> print(numpy.allclose(wa.frequency,freques))
    True
        
    >>> print(numpy.allclose(wa.domega,freques[1]-freques[0]))
    True
    
    >>> tb = wa.get_TimeAxis()
    >>> print(numpy.allclose(tb.time,times))
    True
    
    >>> wb = tb.get_FrequencyAxis()
    >>> print(numpy.allclose(wb.frequency,freques))
    True
    
    
    Upper-half TimeAxis and even number of points
    
    >>> ta = TimeAxis(0.0,10,0.1)
    >>> print(ta.atype=="upper-half")
    True
    
    >>> times = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    >>> print(numpy.allclose(ta.time,times))
    True
    
    >>> wa = ta.get_FrequencyAxis()
    >>> freques = 2.0*numpy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(20,0.1))
    >>> print(numpy.allclose(wa.frequency,freques))
    True
        
    >>> print(numpy.allclose(wa.domega,freques[1]-freques[0]))
    True
    
    >>> tb = wa.get_TimeAxis()
    >>> print(numpy.allclose(tb.time,times))
    True
    
    >>> wb = tb.get_FrequencyAxis()
    >>> print(numpy.allclose(wb.frequency,freques))
    True
    
    
    Upper-half TimeAxis and odd number of points    

    >>> ta = TimeAxis(0.0,11,0.1)    
    >>> times = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    >>> print(numpy.allclose(ta.time,times))
    True
    
    >>> wa = ta.get_FrequencyAxis()
    >>> freques = 2.0*numpy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(2*11,0.1))
    >>> print(numpy.allclose(wa.frequency,freques))    
    True
        
    >>> print(numpy.allclose(wa.domega,freques[1]-freques[0]))    
    True
    
    >>> tb = wa.get_TimeAxis()
    >>> print(numpy.allclose(tb.time,times))
    True
    
    >>> wb = tb.get_FrequencyAxis()
    >>> print(numpy.allclose(wb.frequency,freques))
    True
    
    """
    
    def __init__(self,startValue,noPoints,step,
                 atype="upper-half",w0=0.0): #,complete_with='zero'):

        
        ValueAxis.__init__(self,startValue=startValue,
                           noPoints=noPoints,step=step)
        self.time = self.data
        self.dt = self.step
        self.length = self.noPoints
        self.tmin = self.startValue
        self.tmax = self.endValue
        self.w0 = w0
        
        self.allowed_atypes = ["upper-half","complete"]
        if atype in self.allowed_atypes:         
            self.atype = atype
        else:
            raise Exception("Unknown time axis type")
            
        #self.allowed_upper_half = ["zero","symmetric","anti-symmetric",
        #                           "symmetric-conjugate",
        #                           "anti-symmetric-conjugate"]
        
        
    def get_FrequencyAxis(self):
        """ Returns corresponding FrequencyAxis object 
        
        
        """
        from .frequency import FrequencyAxis
        
        if self.atype == 'complete':
        
            om = numpy.fft.fftshift(
            (2.0*numpy.pi)*numpy.fft.fftfreq(self.length,self.dt))
            
            step = om[1]-om[0]
            start = om[0] + self.w0 

                
            nosteps = len(om)
            t0 = self.time[self.length//2] #+self.length%2]#(self.tmax + self.tmin)/2.0
            
        elif self.atype == 'upper-half':
            om = numpy.fft.fftshift(
                (2.0*numpy.pi)*numpy.fft.fftfreq(2*self.length,self.dt))             
            start = om[0] + self.w0
            step = om[1]-om[0]
            nosteps = len(om)    
            t0 = self.tmin
        else:
            raise Exception("Unknown time axis type")             
        
        return FrequencyAxis(start,nosteps,step,atype=self.atype,t0=t0)
       
       
       
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
        
        
        
        
        