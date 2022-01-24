# -*- coding: utf-8 -*-
"""
    Class representing frequency axis of calculations
    
    
    Examples
    --------

    The default type of the `FrequencyAxis` is `complete`. See the discussion
    of the types in the `TimeAxis` documentation.

    >>> wa = FrequencyAxis(0.0,100,0.05)
    >>> ta = wa.get_TimeAxis()
    >>> print(ta.length)
    100

    The type `upper-half` only refers to the corresponding TimeAxis. Everything
    about the FrequencyAxis remains the same as with `complete`.

    >>> wa = FrequencyAxis(0.0, 100, 0.05, atype = "upper-half")
    >>> ta = wa.get_TimeAxis()
    >>> print(ta.length)
    50

    #>>> print(ta.step,2.0*numpy.pi/(100*wa.step))
    
    >>> print(numpy.allclose(ta.step,2.0*numpy.pi/(100*wa.step)))
    True

    For `complete`, everything should work also for an odd number of points

    >>> wa = FrequencyAxis(0.0,99,0.05)
    >>> ta = wa.get_TimeAxis()
    >>> print(ta.length)
    99

    But `upper-half` throws an exception, because by definition its number of
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
    >>> print(numpy.allclose(wa.data,frequencies))
    True

    >>> ta = wa.get_TimeAxis()
    >>> times = 2.0*numpy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(10,0.1))
    >>> print(numpy.allclose(ta.data,times))
    True

    >>> print(numpy.allclose(ta.step,times[1]-times[0]))
    True

    >>> wb = ta.get_FrequencyAxis()
    >>> print(numpy.allclose(wb.data,frequencies))
    True

    >>> tb = wb.get_TimeAxis()
    >>> print(numpy.allclose(tb.data,times))
    True

    Complete FrequencyAxis and odd number of points

    >>> wa = FrequencyAxis(0.0,11,0.1,atype="complete")
    >>> frequencies = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    >>> print(numpy.allclose(wa.data,frequencies))
    True

    >>> ta = wa.get_TimeAxis()
    >>> times = 2.0*numpy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(11,0.1))
    >>> print(numpy.allclose(ta.data,times))
    True

    >>> print(numpy.allclose(ta.step,times[1]-times[0]))
    True

    >>> wb = ta.get_FrequencyAxis()
    >>> print(numpy.allclose(wb.data,frequencies))
    True

    >>> tb = wb.get_TimeAxis()
    >>> print(numpy.allclose(tb.data,times))
    True

    Upper-half FrequencyAxis and even number of points

    >>> wa = FrequencyAxis(0.0,10,0.1,atype="upper-half")
    >>> frequencies = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    >>> print(numpy.allclose(wa.data,frequencies))
    True

    >>> ta = wa.get_TimeAxis()
    >>> times = 2.0*numpy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(10,0.1))
    >>> print(numpy.allclose(ta.data,times[5:10]))
    True

    >>> print(numpy.allclose(ta.step,times[1]-times[0]))
    True

    >>> wb = ta.get_FrequencyAxis()
    >>> print(numpy.allclose(wb.data,frequencies))
    True

    >>> tb = wb.get_TimeAxis()
    >>> print(numpy.allclose(tb.data,times[5:10]))
    True

    
    Class Details
    -------------

"""
import numpy

from .valueaxis import ValueAxis
from .managers import EnergyUnitsManaged
from .managers import energy_units
from ..utils.types import UnitsManagedRealArray
from ..utils.types import UnitsManagedReal


class FrequencyAxis(ValueAxis, EnergyUnitsManaged):
    """Class representing frequency axis of calculations 

    Parameters
    ----------
    start : float
        start of the FrequencyAxis

    length : int
        number of steps

    step : float
        time step

    atype : string {"complete","upper-half"}
        Axis type

    """

    data = UnitsManagedRealArray("data")
    start = UnitsManagedReal("start")
    step = UnitsManagedReal("step")

    def __init__(self, start=0.0, length=1, step=1.0,
                 atype='complete', time_start=0.0):

        #if step > 0:
        if True:
            self.step = step
        #else:
        #    raise Exception("Parameter step has to be > 0")
        self.start = start
        self.length = length

        super().__init__(start=start,
                         length=length, step=step)

        # This would be the alternative if calling super()__init__ would
        # break the units management
        #self.data = numpy.linspace(start,
        #                           start+(length-1)*step, length,
        #                           dtype=numpy.float)

        self.time_start = time_start

        self.allowed_atypes = ["upper-half", "complete"]
        if atype in self.allowed_atypes:
            self.atype = atype
        else:
            raise Exception("Unknown frequency axis type")

    def copy(self):
        axis = FrequencyAxis(self.start, self.length, self.step,
                             atype=self.atype, time_start=self.time_start)
        return axis
        
        
    def get_TimeAxis(self):
        """Returns the corresponding TimeAxis object

        """
        from .time import TimeAxis

        with energy_units("int"):

            if self.atype == 'complete':

                times = numpy.fft.fftshift(
                    numpy.fft.fftfreq(self.length, self.step/(2.0*numpy.pi)))

                step = times[1]-times[0]
                start = self.time_start  + times[0]
                nosteps = self.length

                frequency_start = self.data[self.length//2]


            elif self.atype == 'upper-half':

                if (self.length % 2) != 0:
                    raise Exception("Cannot create upper-half TimeAxis"
                                    + " from an odd number of points")

                times = numpy.fft.fftshift(
                    (2.0*numpy.pi)*numpy.fft.fftfreq(self.length, self.step))

                start = times[int(self.length/2)] + self.time_start
                nosteps = int(self.length/2)
                step = times[1]-times[0]

                frequency_start = self.data[self.length//2]

            else:
                raise Exception("Unknown frequency axis type")


        return TimeAxis(start, nosteps, step,
                        atype=self.atype, frequency_start=frequency_start)
