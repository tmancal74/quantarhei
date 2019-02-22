


# -*- coding: utf-8 -*-
"""
    Linear array of values which are used as variables of numerical functions 
    and parameter dependent matrices



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
    >>> print("{0:.1f}".format(diff))
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

    Mutual compatibility of the two ValueAxis objects can be tested as follows
    
    >>> va1 = ValueAxis(0.0, 2000, 1.0)
    >>> va2 = ValueAxis(0.0, 1000, 1.0)
    >>> va2.is_subsection_of(va1)
    True
    >>> va1.is_subsection_of(va2)
    False
    >>> va1.is_extension_of(va2)
    True
    >>> va2.is_subsection_of(va1)
    True
    >>> va2.is_extension_of(va1)
    False
    >>> va1.is_equal_to(va2)
    False
    >>> va1.is_equal_to(va1)
    True
    >>> va3 = ValueAxis(0.0, 2000, 1.0)
    >>> va3.is_equal_to(va1)
    True

    Class Details
    -------------

"""

import numpy

from ..utils import Float
from ..utils import Integer
from .saveable import Saveable
from .. import REAL

class ValueAxis(Saveable):
    """Linear array of values which are used as variables of numerical functions 
    and parameter dependent matrices 
    
    Parameters
    ----------
    
    start : float
        Beginning of the axis
        
    lenght : int
        Number of points in the DFunction
        
    step : float
        Step size on the axis
        
    
    """

    step = Float("step")
    start = Float("start")
    length = Integer("length")

    def __init__(self, start=0.0, length=1, step=1.0):

        #if step > 0:
        if True:
            self.step = step
        #else:
        #    raise Exception("Parameter step has to be > 0")
        self.start = start
        self.length = length

        self.data = numpy.linspace(start,
                                   start+(length-1)*step, length,
                                   dtype=REAL)


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


    def is_equal_to(self, axis):
        """Returns True if the axis is equal to this ValueAxis
        
        """
        return ((self.start == axis.start) and (self.step == axis.step)
                and (self.length == axis.length))

        
    def __eq__(self, other):
        return self.is_equal_to(other)


    def is_extension_of(self, axis):
        """Returns True if the axis is contained in this ValueAxis
        
        """
        ret = True
        ret = ret and (self.start <= axis.start)
        ret = ret and (self.step == axis.step)
        ret = ret and (self.max >= axis.max)
        
        found_one_point = False
        N = self.length
        k = 0
        while (not found_one_point) and (k < N):
            point = self.data[k]
            if point in axis.data:
                found_one_point = True
            k += 1
            
        ret = ret and found_one_point
        
        return ret
    
    def is_subsection_of(self, axis):
        """Returns True if the axis contains this ValueAxis
        
        """
        ret = True
        ret = ret and (self.start in axis.data)
        ret = ret and (self.step == axis.step)
        ret = ret and (self.max in axis.data)
        
        return ret
    
    def is_subset_of(self, axis):
        """Returns True if the ValueAxis is a subset of axis
        
        We check if all values of this axis are also values of the submitted 
        axis object.
        
        Examples
        --------
        
        >>> ta1 = ValueAxis(0.0, 1000, 1.0)
        >>> ta2 = ValueAxis(0.0, 490, 2.0)
        >>> ta2.is_subset_of(ta1)
        True
        
        >>> ta1 = ValueAxis(0.0, 1000, 1.0)
        >>> ta2 = ValueAxis(3.0, 900, 1.0)
        >>> ta2.is_subset_of(ta1)
        True

        >>> ta1 = ValueAxis(0.0, 1000, 1.0)
        >>> ta2 = ValueAxis(3.0, 400, 2.0)
        >>> ta2.is_subset_of(ta1)
        True

        >>> ta1 = ValueAxis(0.0, 1000, 1.12345)
        >>> ta2 = ValueAxis(3*1.12345, 400, 2*1.12345)
        >>> ta2.is_subset_of(ta1)
        True
        
        >>> ta1 = ValueAxis(0.0, 1000, 1.12345)
        >>> ta2 = ValueAxis(3.0, 400, 2*1.12345)
        >>> ta2.is_subset_of(ta1)
        False
        
        >>> ta1 = ValueAxis(0.0, 1000, 1.0)
        >>> ta2 = ValueAxis(3.0, 500, 2.0)
        >>> ta2.is_subset_of(ta1)
        False

        """
        ret = True
        
        Nst = round(self.step/axis.step)
        ret = ret and (Nst*axis.step == self.step)
        ret = ret and ((self.start in axis.data) and (self.start < axis.max))
        ret = ret and ((self.max in axis.data))
        
        return ret
    
    
    
    def is_superset_of(self, axis):
        """Returns True if the ValueAxis is a superset of axis
        
        We check if all values of this axis are also values of the submitted 
        axis object.
        
        Examples
        --------
        
        >>> ta1 = ValueAxis(0.0, 1000, 1.0)
        >>> ta2 = ValueAxis(0.0, 490, 2.0)
        >>> ta1.is_superset_of(ta2)
        True
        
        >>> ta1 = ValueAxis(0.0, 1000, 1.0)
        >>> ta2 = ValueAxis(3.0, 900, 1.0)
        >>> ta1.is_superset_of(ta2)
        True

        >>> ta1 = ValueAxis(0.0, 1000, 1.0)
        >>> ta2 = ValueAxis(3.0, 400, 2.0)
        >>> ta1.is_superset_of(ta2)
        True

        >>> ta1 = ValueAxis(0.0, 1000, 1.12345)
        >>> ta2 = ValueAxis(3*1.12345, 400, 2*1.12345)
        >>> ta1.is_superset_of(ta2)
        True
        
        >>> ta1 = ValueAxis(0.0, 1000, 1.12345)
        >>> ta2 = ValueAxis(3.0, 400, 2*1.12345)
        >>> ta1.is_superset_of(ta2)
        False
        
        >>> ta1 = ValueAxis(0.0, 1000, 1.0)
        >>> ta2 = ValueAxis(3.0, 500, 2.0)
        >>> ta1.is_superset_of(ta2)
        False
        
        """
        ret = True
        
        Nst = round(axis.step/self.step)
        ret = ret and (Nst*self.step == axis.step)
        ret = ret and ((axis.start in self.data) and (axis.start < self.max))
        ret = ret and ((axis.max in self.data))        
        
        return ret
    
    