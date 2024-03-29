# -*- coding: utf-8 -*-
import numpy

from ..core.frequency import FrequencyAxis
from ..core.dfunction import DFunction
from ..core.datasaveable import DataSaveable
from ..core.managers import EnergyUnitsManaged
from ..core.units import cm2int




class AbsSpectrumBase(DFunction, EnergyUnitsManaged, DataSaveable):
    """Provides basic container for absorption spectrum
    
    
    Examples
    --------
    
    
    
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
        
    #def add_data(self, data):
    #    self.data += data
        
    def set_by_interpolation(self, x, y, xaxis="frequency"):
        """Sets the data by interpolation with splines
        
        When the spectrum is defined in wavelength, it is converted to 
        an internal representation in frequency.
        
        Examples
        --------
        
        >>> from quantarhei import REAL
        >>> from quantarhei import energy_units
        >>> abs = AbsSpectrumBase()
        >>> x = numpy.array([600.0 + 5.0*ii for ii in range(100)], dtype=REAL)
        >>> y = numpy.exp(-(x-800.0)**2/(50**2))
        >>> abs.set_by_interpolation(x, y, xaxis="wavelength")
        >>> with energy_units("1/cm"):
        ...     print("%6.3f, %6.3f" % (abs.axis.min, abs.axis.max))
        9132.420, 16591.324
        
        """
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
        
        
        Examples
        --------
        
        Standard usage
        
        >>> from quantarhei import energy_units
        >>> abs = AbsSpectrumBase()
        >>> other = AbsSpectrumBase()
        >>> with energy_units("1/cm"):
        ...     w = FrequencyAxis(10000.0, 1000, 5.0)
        ...     y = numpy.exp(-(w.data-12500.0)**2/(200.0**2))
        >>> other.set_axis(w)
        >>> other.set_data(y)
        >>> abs.add_to_data(other)
        
        Axes have to be compatible (i.e. the same)
        
        >>> yetanother = AbsSpectrumBase()
        >>> yetanother.set_axis(FrequencyAxis(8000.0, 1000, 4.0))
        >>> abs.add_to_data(yetanother)
        Traceback (most recent call last):
            ...
        Exception: Incompatible axis
        
        An empty AbsSpectrumBase can be filled by add_to_data() method
        Axis is taken from the class that we add
        
        >>> yetanother.set_axis(w)
        >>> yetanother.add_to_data(abs)
        >>> onemore = AbsSpectrumBase()
        >>> onemore.add_to_data(other)
        
        >>> numpy.testing.assert_allclose(onemore.data,abs.data) 

        
        """

        
        if self.axis is None:
            self.axis = spect.axis.copy()
            
        if not numpy.allclose(spect.axis.data, self.axis.data):
            #numpy.savetxt("spect_data_wrong.dat", spect.axis.data)
            #numpy.savetxt("self_data_wrong.dat", self.axis.data)
            raise Exception("Incompatible axis")
            
        if self.data is None:
            self.data = numpy.zeros(len(spect.data),
                                    dtype=spect.axis.data.dtype)
        
        self.data += spect.data
        
    #save method is inherited from DFunction 
    
    def save_data(self, filename):
        """Saves the data of this absorption spectrum
        
        """
        super().save_data(filename, with_axis=self.axis)


    def load_data(self, filename):
        """Loads data from file into this absorption spectrum
        
        """
        if self.axis is None:
            raise Exception("The property `axis` has to be defined")
        super().load_data(filename, with_axis=self.axis)

        
    def plot(self, **kwargs):
        """ Plotting absorption spectrum using the DFunction plot method
        
        """
        if "ylabel" not in kwargs:
            ylabel = r'$\alpha(\omega)$ [a.u.]'
            kwargs["ylabel"] = ylabel
            
        fig = super().plot(**kwargs)
        if fig is not None:
            return fig


        
