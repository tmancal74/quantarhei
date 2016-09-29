# -*- coding: utf-8 -*-
"""
    Quantarhei package (http://www.github.com/quantarhei)

    ValueAxis module


"""

import numpy

from ..utils import Float
from ..utils import Integer

class ValueAxis:
    """ Linear array of values which are used as variables of numerical
    functions and parameter dependent matrices



    Properties
    ----------

    start : float
        Starting value of the axis

    step : float
        Step between consecutive values of the axis

    length : int
        Length of the `data` array

    data : numpy.array of float
        Values of the axis


    Methods
    -------

    min()
        returns the same value as the attribute `start`

    max()
        returns the last value of the `data` attribute
        (`data[self.length-1]`)

    locate(val)
        returns the index of the element of `data` attribute which has the
        closes lower value than `val`

    nearest(val)
        returns the index of the nearest neighbor to `val`


    Examples
    --------

    Creation of the `ValueAxis` object

    >>> va = ValueAxis(0.0,100,1.0)

    Its standard attributes are the following:

    >>> print(va.length)
    100
    >>> print(va.step)
    1.0
    >>> print(va.start)
    0.0

    The attribute `data` is an array of floats

    >>> print(va.data[3:5])
    [ 3.  4.]

    Attributes `min` and `max` are provided for convenience

    >>> print(va.min)
    0.0
    >>> print(va.max)
    99.0

    You can locate an index of a specific value. The lower neighbor and
    the difference to it is returned.

    >>> i,diff = va.locate(16.3)
    >>> print(i)
    16
    >>> print(diff)
    0.3

    The following returns the nearest index on the axis

    >>> i = va.nearest(16.3)
    >>> print(i)
    16
    >>> i = va.nearest(16.7)
    >>> print(i)
    17
    >>> i = va.nearest(16.5)
    >>> print(i)
    17

    """

    step = Float("step")
    start = Float("start")
    length = Integer("length")

    def __init__(self, start=0.0, length=1, step=1.0):

        if step > 0:
            self.step = step
        else:
            raise Exception("Parameter step has to be > 0")
        self.start = start
        self.length = length

        self.data = numpy.linspace(start,
                                   start+(length-1)*step, length,
                                   dtype=numpy.float)


    @property
    def min(self):
        """Returns the minimum value on the axis

        """
        return self.start

    @property
    def max(self):
        """Returns the maximum value on the axis

        """
        return self.data[self.length-1]


    def locate(self, val):
        """
        Returns the index of the lower neigbor of the ``val``

        Returns the index of the lower neigbor of the value x and the
        remaining distance to the value of ``val``

        Parameters
        ----------

        val : float
            A value within the min and max values of the axis

        """

        # nearest smaller neighbor index
        nsni = numpy.int(numpy.floor((val-self.start)/self.step))

        # if n0 is within bounds calculate distance
        # from the lower neighbor
        if (nsni >= 0) and (nsni < self.length):
            dval = val-self.data[nsni]
            return nsni, dval
        else:
            raise Exception("Value out of bounds")



    def nearest(self, val):
        """
        Returns the index of the nearest neighbor to ``val``

        Returns the index of the nearest neighbor of ``val``. If ``val``
        is out of the upper bounds within the ``self.step`` from the upper
        limit, the lower neighbor index is returned.

        Parameters
        ----------

        val : float
            A value within the min and max values of the axis

        """

        # nearest smaller neighbor index
        nsni = numpy.int(numpy.floor((val-self.start)/self.step))


        if (nsni >= 0) and (nsni < self.length):

            # if n0 is with bounds calculate difference
            # from the lower neighbor
            diff1 = numpy.abs(val-self.data[nsni])

            # if the upper neighbor is within bounds calculate difference
            # from the upper neighbor
            if nsni+1 < self.length:
                diff2 = numpy.abs(val - self.data[nsni+1])
            else:
                diff2 = 5*self.step

            # return the closer neighbor
            if diff1 < diff2:
                return nsni
            else:
                if nsni+1 < self.length:
                    return nsni + 1
                else:
                    # if the upper neigbor is out of bounds
                    # return the lower one
                    return nsni

        else:
            raise Exception("Value out of bounds")
