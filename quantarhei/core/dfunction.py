# -*- coding: utf-8 -*-
import scipy.interpolate
import numpy
import numbers
        
import matplotlib.pyplot as plt

from .valueaxis import ValueAxis
from .time import TimeAxis
from .frequency import FrequencyAxis 

#FIXME Check the posibility to set a derivative of the spline at the edges
#FIXME Enable vectorial arguments and values
#FIXME Make sure we can interpolate complex functions
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
    
    at(x, approx="default")
        Returns the value of the function at a given value of argument `x`. The
        default interpolation is linear, until the spline interpolation is
        initialized by calling the method with approx = "spline". From then
        on, the default is spline.
        
    get_Fourier_transform()
        Returns a Fourier transformed DFunction
       
       
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
    
    Fourier transform of a DFunction
    
    >>> from .time import TimeAxis
    >>> dt = 0.1; Ns = 10000
    >>> t = TimeAxis(-(Ns//2)*dt,Ns,dt,atype="complete")
    >>> gg = 1.0/30.0
    >>> y = numpy.exp(-numpy.abs(t.time)*gg)
    >>> f = DFunction(t,y)
    >>> F = f.get_Fourier_transform()
    >>> print(numpy.allclose(F.at(0.0,approx="spline"),2.0/gg,rtol=1.0e-5))
    True
    
    >>> print(numpy.allclose(F.at(1.0),2*gg/(gg**2 + 1.0**2),rtol=1.0e-3))
    True
    
    >>> print(numpy.allclose(F.at(0.15),2*gg/(gg**2 + 0.15**2),rtol=1.0e-4))
    True
    
    >>> t = TimeAxis(0,Ns,dt)
    >>> print(t.atype == "upper-half")
    True
    
    >>> gg = 1.0/30.0
    >>> y = numpy.exp(-numpy.abs(t.time)*gg)
    >>> f = DFunction(t,y)
    >>> F = f.get_Fourier_transform()
    >>> print(numpy.allclose(F.at(0.0,approx="spline"),2.0/gg,rtol=1.0e-5))
    True
    
    >>> print(numpy.allclose(F.at(1.0),2*gg/(gg**2 + 1.0**2),rtol=1.0e-3))
    True
    
    >>> print(numpy.allclose(F.at(0.15),2*gg/(gg**2 + 0.15**2),rtol=1.0e-4))
    True

    DFunction can be complex valued    
    
    >>> y = numpy.sin(t.time/10.0) + 1j*numpy.cos(t.time/10.0)
    >>> fi = DFunction(t,y)
    >>> print(numpy.allclose(fi.at(13.2,approx="spline"),(numpy.sin(13.2/10.0) \
        + 1j*numpy.cos(13.2/10.0)),rtol=1.0e-7))    
    True
    
    """

    allowed_interp_types = ("linear","spline","default")    
    
    def __init__(self,x,y):
        
        self._has_imag = False

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
                
                if isinstance(self.data[0],numbers.Real):
                    self._has_imag = False
                else:
                    self._has_imag = True    

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
        self._spline_r = \
               scipy.interpolate.UnivariateSpline(
                  self.axis.data,numpy.real(self.data),s=0)
        if self._has_imag:
            self._spline_i = \
               scipy.interpolate.UnivariateSpline(
                  self.axis.data,numpy.imag(self.data),s=0)
                  
        self._splines_initialized = True
        #print("Calculating splines")
        
    def _spline_value(self,x):
        """Returns the splie interpolated value of the function
        
        """
        if self._has_imag:
            ret = self._spline_r(x) + 1j*self._spline_i(x)
        else:
            ret = self._spline_r(x)
        return ret
        
        
    """
 
    Fast Fourier transform    

    """        
        
    def get_Fourier_transform(self):
        """Returns Fourier transform of the DFunction
        
        """
        
        t = self.axis
        y = self.data
        
        if isinstance(t,TimeAxis):
            
            w = t.get_FrequencyAxis()
            
            if t.atype == "complete":
                Y = t.length*numpy.fft.fftshift(numpy.fft.ifft(
                numpy.fft.fftshift(y)))*t.dt
            elif t.atype == "upper-half":
                yy = numpy.zeros(w.length,dtype=y.dtype)
                yy[0:w.length//2] = y
                for k in range(0,t.length-1):
                    yy[w.length-k-1] = numpy.conj(y[k+1])
                Y = 2.0*t.length*numpy.fft.fftshift(numpy.fft.ifft(yy))*t.dt
                
            F = DFunction(w,Y)
            
            
        elif isinstance(t,FrequencyAxis):
            
            w = t
            t = w.get_TimeAxis()
            
            if w.atype == "complete":
                Y = w.length*numpy.fft.fftshift(numpy.fft.ifft(
                numpy.fft.fftshift(y)))*w.domega/(numpy.pi*2.0) 
                
            F = DFunction(t,Y)
            
            
        else:
            pass

                
        return F
    
    
    def get_inverse_Fourier_transform(self):
        """Returns inverse Fourier transform of the DFunction
        
        """
        
        t = self.axis
        y = self.data
        
        if isinstance(t,TimeAxis):
            
            w = t.get_FrequencyAxis()
            
            if t.atype == "complete":
                Y = numpy.fft.fftshift(numpy.fft.fft(
                numpy.fft.fftshift(y)))*t.dt
            elif t.atype == "upper-half":
                yy = numpy.zeros(w.length,dtype=y.dtype)
                yy[0:w.length//2] = y
                for k in range(0,t.length-1):
                    yy[w.length-k-1] = numpy.conj(y[k+1])
                Y = numpy.fft.fftshift(numpy.fft.fft(yy))*t.dt
                
            F = DFunction(w,Y)
            
            
        elif isinstance(t,FrequencyAxis):
            
            w = t
            t = w.get_TimeAxis()
            
            if w.atype == "complete":
                Y = numpy.fft.fftshift(numpy.fft.fft(
                numpy.fft.fftshift(y)))*w.domega/(numpy.pi*2.0) 
                
            F = DFunction(t,Y)
            
            
        else:
            pass

                
        return F
        
        
    def plot(self, title=None,
             title_font=None,
             axis=None,
             xlabel=None,
             ylabel=None,
             label_font=None,
             text=None,
             text_font=None,
             real_only=True,
             show=True,
             color=None):
        """Plotting of the DFunction's data against the ValueAxis 
        
        """


        if color is not None:
            if len(color) == 1:
                clr = [color,color]
            else:
                clr = [color[0],color[1]]
            
        if isinstance(self.data[0],numbers.Complex):
            if color is not None:
                plt.plot(self.axis.data,numpy.real(self.data),clr[0])
                if not real_only:
                    plt.plot(self.axis.data,numpy.imag(self.data),clr[1])
            else:
                plt.plot(self.axis.data,numpy.real(self.data))
                if not real_only:
                    plt.plot(self.axis.data,numpy.imag(self.data))
        else:
            if color is not None:
                plt.plot(self.axis.data,self.data,clr[0])
            else:
                plt.plot(self.axis.data,self.data)
            
        if axis is not None:
            plt.axis(axis)
            
        if title is not None:
            plt.title(title)
            
        if text is not None:
            if text_font is not None:
                plt.text(text[0],text[1],
                     text[2], fontdict=text_font)
            else:
                plt.text(text[0],text[1],text[2])     
        
        if label_font is not None:
            font = label_font
        else:
            font={'size':20}
            
        if xlabel is not None:
            xl = '$\omega$ [fs$^{-1}$]'
            
        if isinstance(self.axis,FrequencyAxis):
            xl = '$\omega$ [fs$^{-1}$]'
            yl = '$F(\omega)$'
        if isinstance(self.axis,TimeAxis):
            xl = '$t$ [fs]'
            yl = '$f(t)$'
        
        if xlabel is not None:
            xl = xlabel
        if ylabel is not None:
            yl = ylabel

        plt.xlabel(xl,**font)
        plt.ylabel(yl,**font)   
            
        if show:
            plt.show()

        
