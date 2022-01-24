# -*- coding: utf-8 -*-
"""
    Class representing time in time dependent calculations.
    
    User level function of the Quantarhei package. To be used as:
        
    >>> import quantarhei as qr   
    >>> time = qr.TimeAxis()

    The `TimeAxis` class stands in a close relation to `FrequencyAxis`
    class of the `quantarhei` package. `FrequencyAxis` represents the
    frequencies one obtains in the Fast Fourier transform of a function of
    the `TimeAxis`. By default, `TimeAxis` is of the type `upper-half` which
    means that by specifying the `start`, `length` and `step` we
    represent the upper half of the interval `<start-length*step,
    start+(length-1)*step>`. The Fourier transform of a time dependent
    object defined on the `TimeAxis` will then have twice as many points as
    the `TimeAxis`. This is usefull when the time dependent object has some
    special symmetries. One example is the so-called quantum bath correlation
    function which fulfills the relation (in LaTeX)

    C(-t) = C^{*}(t)


    Examples
    --------

    Default TimeAxis is of the 'upper-half' type

    >>> ta = TimeAxis(0.0,100,0.1)
    >>> ta.atype
    'upper-half'

    It is defined between the `start` and `start+(length-1)*step`

    >>> '%.3f' % ta.min
    '0.000'

    >>> '%.3f' % ta.max
    '9.900'

    However, when we ask for the corresponding FrequencyAxis, we get an
    array of frequencies which is twice as long and contains `2*length`
    points

    >>> wa = ta.get_FrequencyAxis()
    >>> print(wa.length)
    200

    The frequency step is therefore twice shorter than one would normally
    expect.

    >>> print(wa.step == wa.step)
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
    more straightforward. Now, the `TimeAxis` represents the interval on
    which a time dependent object is defined completely. The number of
    points in the corresponding `FrequencyAxis` is the same as in the
    `TimeAxis` and the frequency step is as expected from a normal FFT.

    >>> wa = ta.get_FrequencyAxis()
    >>> print(wa.length)
    100

    >>> print(wa.step == wa.step)
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
    >>> print(numpy.allclose(ta.data,times))
    True

    >>> wa = ta.get_FrequencyAxis()
    >>> freques = 2.0*numpy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(10,0.1))
    >>> print(numpy.allclose(wa.data,freques))
    True

    >>> print(numpy.allclose(wa.step,freques[1]-freques[0]))
    True

    >>> tb = wa.get_TimeAxis()
    >>> print(numpy.allclose(tb.data,times))
    True

    >>> wb = tb.get_FrequencyAxis()
    >>> print(numpy.allclose(wb.data,freques))
    True

    Complete TimeAxis and odd number of points

    >>> ta = TimeAxis(0.0,11,0.1,atype="complete")
    >>> times = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    >>> print(numpy.allclose(ta.data,times))
    True

    >>> wa = ta.get_FrequencyAxis()
    >>> freques = 2.0*numpy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(11,0.1))
    >>> print(numpy.allclose(wa.data,freques))
    True

    >>> print(numpy.allclose(wa.step,freques[1]-freques[0]))
    True

    >>> tb = wa.get_TimeAxis()
    >>> print(numpy.allclose(tb.data,times))
    True

    >>> wb = tb.get_FrequencyAxis()
    >>> print(numpy.allclose(wb.data,freques))
    True


    Upper-half TimeAxis and even number of points

    >>> ta = TimeAxis(0.0,10,0.1)
    >>> print(ta.atype=="upper-half")
    True

    >>> times = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    >>> print(numpy.allclose(ta.data,times))
    True

    >>> wa = ta.get_FrequencyAxis()
    >>> freques = 2.0*numpy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(20,0.1))
    >>> print(numpy.allclose(wa.data,freques))
    True

    >>> print(numpy.allclose(wa.step,freques[1]-freques[0]))
    True

    >>> tb = wa.get_TimeAxis()
    >>> print(numpy.allclose(tb.data,times))
    True

    >>> wb = tb.get_FrequencyAxis()
    >>> print(numpy.allclose(wb.data,freques))
    True


    Upper-half TimeAxis and odd number of points

    >>> ta = TimeAxis(0.0,11,0.1)
    >>> times = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    >>> print(numpy.allclose(ta.data,times))
    True

    >>> wa = ta.get_FrequencyAxis()
    >>> freques = 2.0*numpy.pi*numpy.fft.fftshift(numpy.fft.fftfreq(2*11,0.1))
    >>> print(numpy.allclose(wa.data,freques))
    True

    >>> print(numpy.allclose(wa.step,freques[1]-freques[0]))
    True

    >>> tb = wa.get_TimeAxis()
    >>> print(numpy.allclose(tb.data,times))
    True

    >>> wb = tb.get_FrequencyAxis()
    >>> print(numpy.allclose(wb.data,freques))
    True


    Class Details
    -------------



"""

import numpy

from .valueaxis import ValueAxis
from .managers import energy_units

class TimeDependent:
    pass

class TimeAxis(ValueAxis):
    """ Class representing time in time dependent calculations.


    Parameters
    ----------
    start : float
        start of the TimeAxis

    length : int
        number of steps

    step : float
        time step

    atype : string {"complete","upper-half"}
        Axis type


    """

    def __init__(self, start=0.0, length=1, step=1.0,
                 atype="upper-half", frequency_start=0.0):


        ValueAxis.__init__(self, start=start,
                           length=length, step=step)

        self.frequency_start = frequency_start

        self.allowed_atypes = ["upper-half", "complete"]
        if atype in self.allowed_atypes:
            self.atype = atype
        else:
            raise Exception("Unknown time axis type")


    def shift_to_zero(self):
        """Shifts the values so that the first one is zero
        
        """
        if self.start != self.data[0]:
            raise Exception("Inconsistent data")
            
        if self.start > 0.0:
            self.data[:] = self.data[:] - self.start
            self.start = 0.0
            
            
    def get_FrequencyAxis(self):
        """ Returns corresponding FrequencyAxis object


        """
        from .frequency import FrequencyAxis

        if self.atype == 'complete':

            frequencies = numpy.fft.fftshift(
                (2.0*numpy.pi)*numpy.fft.fftfreq(self.length, self.step))

            step = frequencies[1]-frequencies[0]
            start = frequencies[0] + self.frequency_start

            nosteps = len(frequencies)
            time_start = self.data[self.length//2]

        elif self.atype == 'upper-half':

            frequencies = numpy.fft.fftshift(
                (2.0*numpy.pi)*numpy.fft.fftfreq(2*self.length, self.step))

            start = frequencies[0] + self.frequency_start
            step = frequencies[1] - frequencies[0]
            nosteps = len(frequencies)
            time_start = self.min

        else:
            raise Exception("Unknown time axis type")

        # this creation has to be protected from units management
        with energy_units("int"):
            faxis = FrequencyAxis(start, nosteps, step,
                                  atype=self.atype, time_start=time_start)

        return faxis


if __name__ == "__main__":
    import doctest
    doctest.testmod()
