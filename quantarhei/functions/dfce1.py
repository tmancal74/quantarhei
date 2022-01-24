# -*- coding: utf-8 -*-

import numpy
import scipy.signal as signal

from .. import DFunction

class Cos(DFunction):
    """Discrete cosine defined on a value axis
    
    Parameters
    ----------
    
    x : ValueAxis
        ValueAxis (TimeAxis, FrequencyAxis) on which the function is defined
        
    omega : float
        Circular frequency 
        
    phi : float
        Phase of the cosine; default is zero
        
        
    Examples
    --------
    
    >>> import quantarhei as qr
    >>> omega = 2.0*3.14159/20.0
    >>> t = 12.3
    >>> time = qr.TimeAxis(0.0, 100, 20.0/100)
    >>> cos = Cos(time, omega)
    >>> a = numpy.cos(omega*t)
    >>> b = cos.at(t, approx="spline")
    >>> numpy.testing.assert_almost_equal(a, b, decimal=7)

    """
    def __init__(self, x, omega, phi=0.0):
        super().__init__()
        self._make_me(x, numpy.cos(omega*x.data + phi))
        
        
class Tukey(DFunction):
    """Tukey window, also known as a tapered cosine window
    
    Parameters
    ----------
    
    x : ValueAxis
        ValueAxis (TimeAxis, FrequencyAxis) on which the function is defined
        
    r : float
        Fraction of the window inside the cosine tapered region
        
    sym : bool, optional
        If True, the window is symmetric. If False, we take only upper half
        of the window. The fraction `r` is always with respect to the whole
        symmetric window.
    
    
    Examples
    --------
    
    >>> import quantarhei as qr
    >>> r = 0.3
    >>> time = qr.TimeAxis(0.0, 1000, 1.0)
    >>> tuw = Tukey(time, r)
    >>> numpy.testing.assert_equal(tuw.at(500.0), 1.0)
    >>> print("tukey at 100.0: {:6.4f}".format(tuw.at(950.0)))
    tukey at 100.0: 0.2414
    
    >>> tuw = Tukey(time, r, sym=False)
    >>> numpy.testing.assert_equal(tuw.at(100.0), 1.0)
    >>> print("tukey at 100.0: {:6.4f}".format(tuw.at(950.0)))
    tukey at 100.0: 0.0170
        
    """
    
    def __init__(self, x, r, sym=True, x_offset=0.0):
        super().__init__()
        L = x.length
        
        if sym:
            data = signal.tukey(L, r, sym=sym)
            
        else:
            data = signal.tukey(2*L, r*2, sym=sym)
            data = data[L:]
            
        ii = 0
        for xi in x.data:
            if xi < x_offset:
                data[ii] = 0.0
            ii += 1
                
        self._make_me(x, data)     