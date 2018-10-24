# -*- coding: utf-8 -*-


import numpy

import matplotlib.pyplot as plt

from ..core.frequency import FrequencyAxis
from ..core.dfunction import DFunction
from ..core.managers import EnergyUnitsManaged
from ..core.units import cm2int


class AbsSpectrumBase(DFunction, EnergyUnitsManaged):
    """Provides basic container for absorption spectrum
    
    """
    
    def __init__(self, axis=None, data=None):
        super().__init__()
        self.axis = axis
        self.data = data
        
    def set_axis(self, axis):
        """Sets axis atribute
        
        Parameters
        ----------
        
        axis : FrequencyAxis object
            Frequency axis object. This object has managed energy units
            
        """
        self.axis = axis
        
    def set_data(self, data):
        """Sets data atribute
        
        Parameters
        ----------
        
        data : array like object (numpy array)
            Sets the data of the absorption spectrum
            
        """
        self.data = data
        
    def add_data(self, data):
        self.data += data
        
    def set_by_interpolation(self, x, y, xaxis="frequency"):
        
        from scipy import interpolate
        
        if xaxis == "frequency":
            
            om = self.convert_2_internal_u(x)
            
        elif xaxis == "wavelength":
            # convert to internal (nano meters) units of wavelength
            
            
            # convert to energy (internal units)
            # to cm
            om = 1.0e-7*x
            # to 1/cm
            om = 1.0/om
            # to 1/fs
            om = om*cm2int
          
        if om[1] > om[2]:
            # reverse order
            om = numpy.flip(om,0)
            y = numpy.flip(y,0)
            
        # equidistant points on the x-axis
        omin = numpy.amin(om)
        omax = numpy.amax(om)
        length = om.shape[0]
        step = (omax-omin)/length
        
        # new frequency axis
        waxis = FrequencyAxis(omin, length, step)
        
        # spline interpolation 
        tck = interpolate.splrep(om, y, s=0)
        ynew = interpolate.splev(waxis.data, tck, der=0)
        
        # setting the axis and data
        self.axis = waxis
        self.data = ynew
        
    
    def clear_data(self):
        """Sets spectrum data to zero
        
        """
        shp = self.data.shape
        self.data = numpy.zeros(shp, dtype=numpy.float64)

    def normalize2(self,norm=1.0):
        """Normalizes spectrum to a given value
        
        """
        mx = numpy.max(self.data)
        self.data = norm*self.data/mx

    def normalize(self):
        """Normalization to one
        
        """
        self.normalize2(norm=1.0)
        
    def subtract(self, val):
        """Subtracts a value from the spectrum to shift its base line
        
        """
        self.data -= val
        

    def add_to_data(self, spect):
        """Performs addition on the data.
        
        Expects a compatible object holding absorption spectrum
        and adds its data to the present absorption spectrum.
        
        Parameters
        ----------
        
        spect : spectrum containing object
            This object should have a compatible axis and some data
        
        """

        
        if self.axis is None:
            self.axis = spect.axis.copy()
            
        if not numpy.allclose(spect.axis.data, self.axis.data):
            numpy.savetxt("spect_data_wrong.dat", spect.axis.data)
            numpy.savetxt("self_data_wrong.dat", self.axis.data)
            raise Exception("Incompatible axis")
            
        if self.data is None:
            self.data = numpy.zeros(len(spect.data),
                                    dtype=spect.axis.data.dtype)
        
        self.data += spect.data
        
        
    def load_data(self, filename, ext=None, replace=False):
        """Load the spectrum from a file
        
        Uses the load method of the DFunction class to load the absorption
        spectrum from a file. It sets the axis type to 'frequency', otherwise
        no changes to the inherited method are applied.
        
        Parameters
        ----------
        
        """
        super().load_data(filename, ext=ext, axis='frequency', replace=replace)

    #save method is inherited from DFunction 
    
        
        
    def plot(self, **kwargs):
        """ Plotting absorption spectrum using the DFunction plot method
        
        """
        if "ylabel" not in kwargs:
            ylabel = r'$\alpha(\omega)$ [a.u.]'
            kwargs["ylabel"] = ylabel
            
        fig = super().plot(**kwargs)
        if fig is not None:
            return fig


        
    def gaussian_fit(self, N=1, guess=None, plot=False, Nsvf=251):
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
            return _n_gaussians(x, N, *p)
        
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
        
        return popt, pcov
        
#    def convert_to_energy(self, eaxis, units):
#        """
#        
#        """
#        
#        if units == "nm":
#            x = self.axis.data
#            y = self.data
#            
#            # to cm
#            x = 1.0e-7*x
#            # to 1/cm
#            x = 1.0/x
#            # to rad/fs
#            x = x*cm2int
#            
#            xn = numpy.zeros(x.shape, dtype=x.dtype)
#            yn = numpy.zeros(y.shape, dtype=y.dtype) 
#            
#            for i in range(len(x)):
#                xn[i] = x[len(x)-i-1]
#                yn[i] = y[len(x)-i-1]
#                
#            # spline it
#            
#            # evaluate at points if eaxis
#
            
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


def _n_gaussians(x, N, *params):
    """Sum of N Gaussian functions plus an offset from zero

    Parameters
    ----------
    
    x : float
        values to calculate Gaussians function at        

    N : int
        number of Gaussians
        
    params : floats
        3*N + 1 parameters corresponding to height, center, fwhm  for each 
        Gaussian and one value of offset
        
    """
    n = len(params)
    k = n//3
    
    if (k*3 == n) and (k == N):
        
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

