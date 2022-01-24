# -*- coding: utf-8 -*-
"""
    Discrete function with interpolation
    
    User level function of the Quantarhei package. To be used as:
        
    >>> import quantarhei as qr   
    >>> f = qr.DFunction()

    Discrete representation of a function with several modes of interpolation.
    Once defined, the function values are obtained by the method at(x), which
    takes the function argument and optionally a specification of the
    interpolation type. See examples below.

    The linear interpolation is the default initially. Once the function is
    interpolated by splines (which happens why you call it with
    approx="spline"), the default switches to "spline". You can always enforce
    the type of interpolation by specifying it explicitely by the `approx`
    argument.


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
    >>> y = numpy.exp(-numpy.abs(t.data)*gg)
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
    >>> y = numpy.exp(-numpy.abs(t.data)*gg)
    >>> f = DFunction(t,y)
    >>> F = f.get_Fourier_transform()
    >>> print(numpy.allclose(F.at(0.0,approx="spline"),2.0/gg,rtol=1.0e-5))
    True

    >>> print(numpy.allclose(F.at(1.0),2*gg/(gg**2 + 1.0**2),rtol=1.0e-3))
    True

    >>> print(numpy.allclose(F.at(0.15),2*gg/(gg**2 + 0.15**2),rtol=1.0e-4))
    True

    DFunction can be complex valued

    >>> y = numpy.sin(t.data/10.0) + 1j*numpy.cos(t.data/10.0)
    >>> fi = DFunction(t,y)
    >>> print(numpy.allclose(fi.at(13.2,approx="spline"),\
    (numpy.sin(13.2/10.0) + 1j*numpy.cos(13.2/10.0)),rtol=1.0e-7))
    True


    Class Details
    -------------

"""

import os

import scipy.interpolate
import numpy
import numbers

import matplotlib.pyplot as plt

from .valueaxis import ValueAxis
from .time import TimeAxis
from .frequency import FrequencyAxis
from .saveable import Saveable
from .datasaveable import DataSaveable
from .. import REAL


#FIXME Check the posibility to set a derivative of the spline at the edges
#FIXME Enable vectorial arguments and values
class DFunction(Saveable, DataSaveable):
    """Discrete function with interpolation

    Parameters
    ----------

    x : ValueAxis (such as TimeAxis, FrequencyAxis etc.)
        Array of the values of the argument of the discrete function

    y : numpy.ndarray
        Array of the function values

    """

    allowed_interp_types = ("linear", "spline", "default")

    def __init__(self, x=None, y=None):


        self._has_imag = None
        self._is_empty = False

        if not ((x is None) and (y is None)):
            self._make_me(x, y)
        else:
            self._is_empty = True

        self._splines_initialized = False

    def _make_me(self, x, y):
        """Creates the DFunction internals
        
        
        """

        self._has_imag = False

        if isinstance(x, ValueAxis):
            self.axis = x
        else:
            raise Exception("First argument has to be of a ValueAxis type")

        if isinstance(y, numpy.ndarray):

            if len(y.shape) == 1:
                if y.shape[0] != self.axis.length:
                    raise Exception("Wrong number of elements"
                    + " in 1D numpy.ndarray")

                #Set values of the function
                self.data = y

                if isinstance(self.data[0], numbers.Real):
                    self._has_imag = False
                else:
                    self._has_imag = True

            else:
                raise Exception("Second argument has to be"
                + " one-dimensional numpy.ndarray")

        else:
            raise Exception("Second argument has to be"
            + " one-dimensional numpy.ndarray")


    def _add_me(self, x, y):
        """Adds data to the DFunction
        
        """
        # if the DFunction is not initialized, call _make_me
        if self._has_imag is None:
            self._make_me(x,y)
            
        # otherwise do the job of adding
        else:
            
            if isinstance(x, ValueAxis):
                xaxis = x
            else:
                raise Exception("First argument has to be of a ValueAxis type")
    
            if isinstance(y, numpy.ndarray):
    
                if len(y.shape) == 1:
                    if y.shape[0] != xaxis.length:
                        raise Exception("Wrong number of elements"
                        + " in 1D numpy.ndarray")
    
                    #Set values of the function
                    data = y
    
                    if not isinstance(data[0], numbers.Real):
                        self._has_imag = True
    
                else:
                    raise Exception("Second argument has to be"
                    + " one-dimensional numpy.ndarray")
    
            else:
                raise Exception("Second argument has to be"
                + " one-dimensional numpy.ndarray")

            # check axis
            if self.axis == xaxis:
                # add data
                self.data += data
            else: 
                raise Exception("On addition, axis objects have to be"
                                +" identical")

    
    def change_axis(self, axis):
        """Replaces the axis object with a compatible one, zero pads or trims the values
        
        """
        if self.axis.is_equal_to(axis):
            # simple replacement
            self.axis = axis
            
        elif self.axis.is_extension_of(axis):
            
            # we trim to the new axis
            ndata = numpy.zeros(axis.length, dtype=self.data.dtype)
            for ii in range(axis.length):
                ndata[ii] = self.at(axis.data[ii])
                
            self.__init__(x=axis, y=ndata)
                
        elif self.axis.is_subsection_of(axis):
            # we zero pad the values
            ndata = numpy.zeros(axis.length, dtype=self.data.dtype)
            for ii in range(axis.length):
                x = axis.data[ii]
                if (x >= self.axis.min) and (x <= self.axis.max):
                    ndata[ii] = self.at(axis.data[ii])
                else:
                    ndata[ii] = 0.0
                
            self.__init__(x=axis, y=ndata)
        
        else:
            raise Exception("Incompatible axis")
        

    def at(self, x, approx="default"):
        """Returns the function value at the argument `x`
        
        Returns the value of the function at a given value of argument `x`. The
        default interpolation is linear, until the spline interpolation is
        initialized by calling the method with approx = "spline". From then
        on, the default is spline.

        Parameters
        ----------

        x : number
            Function argument

        approx : string {"default","linear","spline"}
            Type of interpolation

        """

        if approx not in self.allowed_interp_types:
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

    #
    #
    # Implementations of various interpolation types
    #
    #

    def _get_linear_approx(self, x_in):
        """Returns linear interpolation of the function

        """
        # FIXME: we need qr.FLOAT or more flexible, here
        try:
            ln = len(x_in)
        except:
            ln = 0
        
        if ln > 0:
                
            val = numpy.zeros(len(x_in), dtype=self.data.dtype)
            k_i = 0
            for x in x_in:
#                n,dval = self.axis.locate(x)
#                if n+1 >= self.axis.length:
#                    val[k_i] = self.data[n] \
#                    + dval/self.axis.step*(self.data[n]-self.data[n-1])
#                else:
#                    val[k_i] = self.data[n] \
#                    + dval/self.axis.step*(self.data[n+1]-self.data[n])
                val[k_i] = self._approx_point(x)
                k_i += 1
        else:
            val = self._approx_point(x_in)
            
        return val


    def _approx_point(self, x):
        
        n, dval = self.axis.locate(x)
        if n+1 >= self.axis.length:
            val = self.data[n] \
            + dval/self.axis.step*(self.data[n]-self.data[n-1])
        else:
            val = self.data[n] \
            + dval/self.axis.step*(self.data[n+1]-self.data[n])
        return val
    
        
    def _get_spline_approx(self, x):
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
                  self.axis.data, numpy.real(self.data),s=0)
        if self._has_imag:
            self._spline_i = \
               scipy.interpolate.UnivariateSpline(
                  self.axis.data, numpy.imag(self.data),s=0)

        self._splines_initialized = True
        #print("Calculating splines")

    def _spline_value(self, x):
        """Returns the splie interpolated value of the function

        """
        if self._has_imag:
            ret = self._spline_r(x) + 1j*self._spline_i(x)
        else:
            ret = self._spline_r(x)
        return ret

    def __add__(self, other):
        """Adding two DFunctions together
        
        """
        f = DFunction(self.axis, self.data.copy())
        
        if other.axis == self.axis:
            f.data += other.data
        else:
            raise Exception("axis attribute of both"
                            + " functions has to be identical")
            
        return f
    
    def add_to_data(self, other):
        pass
    

    def apply_to_data(self, func):
        """Applies a submitted function to the data
        
        """
        self.data = func(self.data)
        if self._splines_initialized:
            self._splines_initiated = False


    #
    #
    # Fast Fourier transform
    #
    #

    def get_Fourier_transform(self, window=None):
        """Returns Fourier transform of the DFunction

        """

        t = self.axis
        y = self.data

        if isinstance(t, TimeAxis):
            
            if window is None:
                winfce = DFunction(self.axis, 
                                   numpy.ones(self.axis.length, dtype=REAL))
            else:
                winfce = window
                
            y = y*winfce.data

            w = t.get_FrequencyAxis()

            if t.atype == "complete":

                Y = t.length*numpy.fft.fftshift(numpy.fft.ifft(
                numpy.fft.fftshift(y)))*t.step

            elif t.atype == "upper-half":

                yy = numpy.zeros(w.length, dtype=y.dtype)
                yy[0:w.length//2] = y
                # fill the negative side of the axis according to the symmetry
                # FIXME: No choice of symmetry provided !!!!
                for k in range(0, t.length-1):
                    yy[w.length-k-1] = numpy.conj(y[k+1])

                Y = 2.0*t.length*numpy.fft.fftshift(numpy.fft.ifft(yy))*t.step

            else:
                raise Exception("Unknown axis type"
                                +" (must be complete or upper-half)")

            F = DFunction(w, Y)


        elif isinstance(t, FrequencyAxis):

            w = t
            t = w.get_TimeAxis()

            Y = w.length*numpy.fft.fftshift(numpy.fft.ifft(
                numpy.fft.fftshift(y)))*w.step/(numpy.pi*2.0)

            if w.atype == "complete":

                F = DFunction(t, Y)

            elif w.atype == "upper-half":

                Y = Y[t.length:2*t.length]
                F = DFunction(t, Y)

            else:
                raise Exception("Unknown axis type"
                                +" (must be complete or upper-half)")


        else:
            pass


        return F


    def get_inverse_Fourier_transform(self):
        """Returns inverse Fourier transform of the DFunction

        """

        t = self.axis
        y = self.data

        if isinstance(t, TimeAxis):

            w = t.get_FrequencyAxis()

            if t.atype == "complete":

                #Y = t.length*numpy.fft.fftshift(numpy.fft.ifft(
                #    numpy.fft.fftshift(y)))*t.step
                Y = numpy.fft.fftshift(numpy.fft.fft(
                    numpy.fft.fftshift(y)))*t.step

            elif t.atype == "upper-half":

                yy = numpy.zeros(w.length, dtype=y.dtype)
                yy[0:w.length//2] = y
                # fill the negative side of the axis according to the symmetry
                # FIXME: No choice of symmetry provided !!!!
                for k in range(0, t.length-1):
                    yy[w.length-k-1] = numpy.conj(y[k+1])

                Y = 2.0*numpy.fft.fftshift(numpy.fft.fft(yy))*t.step

            else:
                raise Exception("Unknown axis type"
                                +" (must be complete or upper-half)")

            F = DFunction(w, Y)


        elif isinstance(t, FrequencyAxis):

            w = t
            t = w.get_TimeAxis()

            Y = numpy.fft.fftshift(numpy.fft.fft(
            numpy.fft.fftshift(y)))*w.step/(numpy.pi*2.0)

            if t.atype == "complete":

                F = DFunction(t, Y)

            elif t.atype == "upper-half":

                # this is a temporary fix - whole this part needs
                # rethinking
                y = numpy.zeros(t.length,dtype=numpy.complex128)
                
                # Vlada's suggestion
                #y[0:t.length-1] = Y[t.length+1:2*t.length]
                #y[t.length-1] = 0.0
                
                y[0:t.length] = Y[t.length:2*t.length]
                F = DFunction(t, y)

            else:
                raise Exception("Unknown axis type"
                                +" (must be complete or upper-half)")

        else:
            pass


        return F


    def fit_exponential(self, guess=None):
        """Exponential fit of the function
        
        """
        from scipy.optimize import curve_fit
        
        if guess is None:
            # we make a guess with one exponential of 100 fs decay time
            guess = [self.data[0], 1.0/100.0, 0.0]
            
        popt, pcov = curve_fit(_exp_fcion, self.axis.data, self.data, p0=guess)
        
        return popt
        

    def fit_gaussian(self, N=1, guess=None, plot=False, Nsvf=251):
        from scipy.signal import savgol_filter
        from scipy.interpolate import UnivariateSpline
        """Performs a Gaussian fit of the spectrum based on an initial guess
        
        
        Parameters
        ----------
        
        Nsvf : int
            Length of the Savitzky-Golay filter window (odd integer)
            
            
        """
        x = self.axis.data
        y = self.data
        
        if guess is None:
            
            raise Exception("Guess is required at this time")
            # FIXME: create a reasonable guess
            guess = [1.0, 11000.0, 300.0, 0.2,
                     11800, 400, 0.2, 12500, 300]
            
            #
            # Find local maxima and guess their parameters
            #

            # Fit with a given number of Gaussian functions
            
            if not self._splines_initialized:
                self._set_splines()
            
            # get first derivative and smooth it
            der = self._spline_r.derivative()
            y1 = der(x)
            y1sm = savgol_filter(y1,Nsvf,polyorder=3)
        
            # get second derivative and smooth it
            y1sm_spl_der = UnivariateSpline(x,y1sm,s=0).derivative()(x)
            y2sm = savgol_filter(y1sm_spl_der,Nsvf,polyorder=3)
        
            # find positions of optima by looking for zeros of y1sm
        
        
            # isolate maxima by looking at the value of y2sm
        

            #plt.plot(x, der(x))
            #plt.plot(x, y1sm)
            plt.plot(x, y2sm)
            plt.show()
        
        
        
        def funcf(x, *p):
            return _n_gaussians(x, *p)
        
        # minimize, leastsq,
        from scipy.optimize import curve_fit            
        popt, pcov = curve_fit(funcf, x, y, p0=guess)
        
        if plot:
        
            plt.plot(x,y)
            plt.plot(x,_n_gaussians(x, N, *popt))
            for i in range(N):
                a = popt[3*i]
                print(i, a)
                b = popt[3*i+1]
                c = popt[3*i+2]
                y = _gaussian(x, a, b, c)
                plt.plot(x, y,'-r')
            plt.show()
        
        # FIXME: Create a readable report
        
        return popt #, pcov
            


    def plot(self, fig=None, title=None,
             title_font=None,
             axis=None,
             vmax=None,
             vmin=None,
             xlabel=None,
             ylabel=None,
             label_font=None,
             text=None,
             text_font=None,
             label = None,
             text_loc=[0.05,0.9],
             fontsize="20",
             real_only=True,
             show=False,
             color=None):
        """Plotting of the DFunction's data against the ValueAxis.
        
        
        Parameters
        ----------
        
        title : str
            Title of the plot
            
        title_font : str
            Name of the title font
            

        """
        #
        #  How to treat the figures
        #
        if fig is None:
            fig, ax = plt.subplots(1,1)
        else:
            fig.clear()
            fig.add_subplot(1,1,1)
            ax = fig.axes[0]

        if axis is None:
            if vmax is None:
                hore = numpy.amax(numpy.real(self.data))
            else:
                hore = vmax
            if vmin is None:
                dole = numpy.amin(numpy.real(self.data))
            else:
                dole = vmin

            levo = self.axis.min
            prvo = self.axis.max

        else:
            levo = axis[0]
            prvo = axis[1]
            dole = axis[2]
            hore = axis[3]
        
        #
        # Label
        #
        pos = text_loc
        if label is not None:
            label = label    
            ax.text((prvo-levo)*pos[0]+levo,
                (hore-dole)*pos[1]+dole,
                label,
                fontsize=str(fontsize))
            
        if color is not None:
            if (len(color) == 1) or isinstance(color, str):
                clr = [color, color]
            else:
                clr = [color[0], color[1]]

        if isinstance(self.data[0], numbers.Complex):
            if color is not None:
                plt.plot(self.axis.data, numpy.real(self.data), clr[0])
                if not real_only:
                    plt.plot(self.axis.data, numpy.imag(self.data), clr[1])
            else:
                plt.plot(self.axis.data, numpy.real(self.data))

                if not real_only:
                    plt.plot(self.axis.data, numpy.imag(self.data))

        else:
            if color is not None:
                plt.plot(self.axis.data, self.data,clr[0])
            else:
                plt.plot(self.axis.data, self.data)

        if axis is not None:
            plt.axis(axis)

        if title is not None:
            plt.title(title)

        if text is not None:
            if text_font is not None:
                plt.text(text[0],text[1],
                     text[2], fontdict=text_font)
            else:
                plt.text(text[0], text[1], text[2])

        if label_font is not None:
            font = label_font
        else:
            font={'size':20}

        if xlabel is None:
            xl = ""
        if ylabel is None:
            yl = ""
            
        if xlabel is not None:
            xl = r'$\omega$ [fs$^{-1}$]'

        if isinstance(self.axis, FrequencyAxis):
            units = self.axis.unit_repr_latex()
            xl = r'$\omega$ ['+units+']'
            yl = r'$F(\omega)$'
        if isinstance(self.axis, TimeAxis):
            xl = r'$t$ [fs]'
            yl = r'$f(t)$'

        if xlabel is not None:
            xl = xlabel
        if ylabel is not None:
            yl = ylabel

        if xl is not None:
            plt.xlabel(xl, **font)
        if xl is not None:
            plt.ylabel(yl, **font)

        if show:
            plt.show()
            
            
    def savefig(self, filename):
        """Saves current figure into a file
        
        
        """
        
        fig = plt.gcf()
        fig.savefig(filename, bbox_inches='tight')
        
        

    def _fname_ext(self, filename, ext):
        """
        
        """
        # find out if filename has an extension
        fname, ex = os.path.splitext(filename)
        ex = ex[1:]
        has_ext = not (ex[1:] == '')

        if has_ext and (ext is None):
            fname = filename
            ext = ex
        elif (not has_ext) and (ext is None):
            ext = 'npy'
        elif ext is not None:
            fname = filename + "." + ext
            
        return fname, ext
            

def _exp_fcion(t, *params):
    
    np = len(params)
    ret = 0.0
    nexp = int((np - 1)/2)
    kp = 0
    for kk in range(nexp):
        #print(kp, params[kp])
        #print(kp+1,params[kp+1], t)
        ret += params[kp]*numpy.exp(-params[kp+1]*t)
        kp += 2
    #print(kp, params[kp])
    ret += params[kp]
    
    return ret

            
def _gaussian(x, height, center, fwhm, offset=0.0):
    """Gaussian function with a possible offset
    
    
    Parameters
    ----------
    
    x : float array
        values to calculate Gaussian function at
        
    height : float
        height of the Gaussian at maximum
        
    center : float
        position of maximum
        
    fwhm : float
        full width at half maximum of the Gaussian function
        
    offset : float
        the value at infinity; effectively an offset on the y-axis
        
    
    """
    
    return height*numpy.exp(-(((x - center)**2)*4.0*numpy.log(2.0))/
                            (fwhm**2)) + offset   


def _n_gaussians(x, *params):
    """Sum of N Gaussian functions plus an offset from zero

    Parameters
    ----------
    
    x : float
        values to calculate Gaussians function at        
        
    params : floats
        3*N + 1 parameters corresponding to height, center, fwhm  for each 
        Gaussian and one value of offset
        
    """
    n = len(params)
    k = n//3
    
    if k*3 + 1 == n:       
        res = 0.0
        pp = numpy.zeros(3)
        for i in range(k):
            pp[0:3] = params[3*i:3*i+3]
            #pp[3] = 0.0
            arg = tuple(pp)
            res += _gaussian(x, *arg)
        res += params[n-1] # last parameter is an offset
        return res
            
    else:
        raise Exception("Inconsistend number of parameters")        

