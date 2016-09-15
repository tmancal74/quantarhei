# -*- coding: utf-8 -*-

import scipy.interpolate
import numpy

from .valueaxis import ValueAxis 

#FIXME Check the posibility to set a derivative of the spline at the edges
#FIXME Enable vectorial arguments and values
class DFunction:
    """
    Discrete function with interpolation
    
    Discrete representation of a function with several modes of interpolation. 
    Once defined, the function values are obtained by the method at(x), which
    takes the function argument and optionally a specification of the 
    interpolation type. See examples below.
    
    The linear interpolation is the default initially. Once the function is
    interpolated by splines (which happens why you call it with 
    approx="spline"), the default switches to "spline". You can always enforce
    the type of interpolation by specifying it explicitely by the `approx`
    argument.
    
    Parameters
    ----------
    
    x : ValueAxis (such as TimeAxis, FrequencyAxis etc.)
        Array of the values of the argument of the discrete function
        
    y : numpy.ndarray
        Array of the function values


    Attributes
    ----------
    
    allowed_interp_types : tuple
        Lists allowed interpolation types. Currently `default`, `linear` 
        and `spline`.
        
    Methods
    -------
    
    at(x,approx="default")
        Returns the value of the function at a given value of argument `x`. The
        default interpolation is linear, until the spline interpolation is
        initialized by calling the method with approx = "spline". From then
        on, the default is spline.
       
       
    Examples
    --------
    
    >>> import numpy
    >>> x = numpy.linspace(0.0,95.0,20)
    >>> y = numpy.exp(-x/30.0)
    >>> u = ValueAxis(0.0,len(x),x[1]-x[0])
    >>> f = DFunction(u,y)
    >>> "%.4f" % f.axis.step
    '5.0000'
    
    Values at the points where the function was defined are exact
    
    >>> "%.4f" % f.at(0.0)
    '1.0000'
    >>> "%.4f" % f.at(5.0)
    '0.8465'
    
    This is the exact value between the discrete points on which the function 
    was defined
    
    >>> "%.4f" % numpy.exp(-2.0/30.0)
    '0.9355'

    Default linear approximation leads to a difference at the second digit
    
    >>> "%.4f" % f.at(2.0)
    '0.9386'

    Spline approximation is much better
    
    >>> "%.4f" % f.at(2.0,approx='spline')
    '0.9355'
    
    """

    allowed_interp_types = ("linear","spline","default")    
    
    def __init__(self,x,y):

        if isinstance(x,ValueAxis):
            self.axis = x
        else:
            raise Exception("First argument has to be of a ValueAxis type")

        if isinstance(y,numpy.ndarray):
            
            if len(y.shape) == 1:
                if y.shape[0] != self.axis.noPoints:
                    raise Exception("Wrong number of elements"
                    + " in 1D numpy.ndarray")
                """
                Set values of the function
                """
                self.data = y

            else:
                raise Exception("Second argument has to be"
                + " one-dimensional numpy.ndarray")

        else:
            raise Exception("Second argument has to be"
            + " one-dimensional numpy.ndarray")
         
        self._splines_initialized = False         
         
        
    
    def at(self,x,approx="default"):
        """Returns the function value at the argument `x`
        
        Parameters
        ----------
        
        x : number
            Function argument
            
        approx : string {"default","linear","spline"}
            Type of interpolation 
        
        """
        
        if not approx in self.allowed_interp_types:
            raise Exception("Unknown interpolation type")
            
        if approx == "default":
            if not self._splines_initialized:
                approx = "linear"
            else:
                approx = "spline"                
                
        if approx == "linear":
            return self._get_linear_approx(x)
        elif approx == "spline":
            return self._get_spline_approx(x)

    """
    
    Implementations of various interpolation types

    """        
    def _get_linear_approx(self,x):
        """Returns linear interpolation of the function
        
        """
        n,dval = self.axis.locate(x)
        if n+1 >= self.axis.noPoints:
            val = self.data[n] \
            + dval/self.axis.step*(self.data[n]-self.data[n-1])
        else:
            val = self.data[n] \
            + dval/self.axis.step*(self.data[n+1]-self.data[n])
        
        return val
        
    def _get_spline_approx(self,x):
        """Returns spline interpolation of the function
        
        """
        if not self._splines_initialized:
            self._set_splines()
        return self._spline_value(x)
            
    def _set_splines(self):
        """Calculates the spline representation of the function 
        
        
        """
        self._spline = \
               scipy.interpolate.UnivariateSpline(\
                  self.axis.data,self.data,s=0)
        self._splines_initialized = True
        #print("Calculating splines")
        
    def _spline_value(self,x):
        """Returns the splie interpolated value of the function
        
        """
        return self._spline(x)

