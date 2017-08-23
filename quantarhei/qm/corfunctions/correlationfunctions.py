# -*- coding: utf-8 -*-
"""
    Quantarhei package (http://www.github.com/quantarhei)

    correlationfunctions module


"""

import numpy
import scipy.interpolate as interp

from ...core.dfunction import DFunction
from ...core.units import kB_intK
from ...core.managers import UnitsManaged
from ...core.managers import energy_units
from ...core.time import TimeAxis
from ...core.frequency import FrequencyAxis


class CorrelationFunction(DFunction, UnitsManaged):
    """Provides typical Bath correlation function types.

    Most important types of bath or energy gap correlation functions are
    provided. Where possible, the correlation function is calculated
    from the parameters from analytical formulae. Where such formulae are
    not available, correlation function is calculated by transformation
    of the spectral density.

    Parameters
    ----------

    axis : TimeAxis
        TimeAxis object specifying the time interval on which the
        correlation function is defined.

    params : dictionary
        A dictionary of the correlation function parameters

    values : optional
        Correlation function can be set by specifying values at all times

    Methods
    -------

    is_analytical()
        Returns `True` if the correlation function is calculated from an
        analytical formula, `False` otherwise.

    copy()
        Returns a copy of the CorrelationFunction object

    get_temperature()
        Returns the temperature of the correlation function

    get_reorganization_energy()
        Returns the reorganization energy parameters of
        the correlation function

    measure_reorganization_energy()
        Calculates reorganization energy from the shape of the correlation
        function

    get_FTCorrelationFunction()
        Returns the Fourier transform of the correlation function

    get_EvenFTCorrelationFunction()
        Returns the Fourier transform of the real part of the correlation
        function

    get_OddFTCorrelationFunction()
        Returns the Fourier transform of the imaginary part of the correlation
        function

    get_SpectralDensity()
        Returns numerically calculated spectral density


    Types of correlation function provided
    --------------------------------------
    OverdampedBrownian-HighTemperature :
        OverdampedBrownian oscillator in high temperature limit

    OverdampedBrownian :
        General overdampedBrownian oscillator

    Examples
    --------

    >>> from quantarhei import TimeAxis
    >>> params = dict(ftype="OverdampedBrownian", cortime=100, reorg=20, T=300)
    >>> time = TimeAxis(0.0,1000,1.0)
    >>> with energy_units("1/cm"):
    ...     cf = CorrelationFunction(time,params)

    >>> with energy_units("1/cm"):
    ...     print(cf.get_reorganization_energy())
    20.0

    Reorganization energy of a correlation function can be calculated from the
    shape of the spectral density by integrating over it. The accuracy
    of such estimation depends on numerics, hence the relative tolerance of
    only 1.0e-4 below

    >>> lamb_definition = cf.get_reorganization_energy()
    >>> lamb_measured = cf.measure_reorganization_energy()
    >>> print(numpy.allclose(lamb_definition, lamb_measured, rtol=1.0e-4))
    True

    """

    allowed_types = ("OverdampedBrownian-HighTemperature",
                     "OverdampedBrownian", 
                     "UnderdampedBrownian",
                     "Value-defined")

    analytical_types = ("OverdampedBrownian-HighTemperature",
                        "OverdampedBrownian")

    def __init__(self, axis, params, values=None):
        super().__init__()

        if not isinstance(axis, TimeAxis):
            taxis = axis.get_TimeAxis()
            self.axis = taxis
        else:
            self.axis = axis

        self._is_composed = False
        
        self.lamb = -1.0
        self.temperature = -1.0
        self.cutoff_time = None

        try:
            ftype = params["ftype"]
            if ftype in CorrelationFunction.allowed_types:
                self.ftype = ftype
            else:
                raise Exception("Unknown Correlation Function Type")

            # we need to save the defining energy units
            self.energy_units = self.manager.get_current_units("energy")
            # because params are set in these units
            self.params = params

        except:
            raise Exception()

        if self.ftype == "OverdampedBrownian-HighTemperature":

            self._make_overdamped_brownian_ht(params, values=values)

        elif self.ftype == "OverdampedBrownian":

            self._make_overdamped_brownian(params, values=values)
            
        elif self.ftype == "UnderdampedBrownian":
            
            self._make_underdamped_brownian(params, values=values)

        elif self.ftype == "Value-defined":

            self._make_value_defined(params, values)

        else:
            raise Exception("Unknown correlation function type of"+
                            "type domain combination.")

    def _matsubara(self, kBT, ctime, nof):
        msf = 0.0
        nut = 2.0*numpy.pi*kBT
        time = self.axis.data
        for i in range(0, nof):
            n = i+1
            msf += nut*n*numpy.exp(-nut*n*time)/((nut*n)**2-(1.0/ctime)**2)
        return msf

    def _make_overdamped_brownian(self, params, values=None):

        temperature = params["T"]
        ctime = params["cortime"]
        lamb = self.manager.iu_energy(params["reorg"],
                                      units=self.energy_units)

        if "matsubara" in params.keys():
            nmatsu = params["matsubara"]
        else:
            nmatsu = 10

        kBT = kB_intK*temperature
        time = self.axis.data

        if values is not None:

            cfce = values
        else:
            cfce = (lamb/(ctime*numpy.tan(1.0/(2.0*kBT*ctime))))\
                *numpy.exp(-time/ctime) \
                - 1.0j*(lamb/ctime)*numpy.exp(-time/ctime)

            cfce += (4.0*lamb*kBT/ctime) \
                *self._matsubara(kBT, ctime, nmatsu)

        self._make_me(self.axis, cfce)
        
        self.lamb = lamb
        self.temperature = temperature
        self.cutoff_time = 5.0*ctime


    def _make_overdamped_brownian_ht(self, params, values=None):
        temperature = params["T"]
        ctime = params["cortime"]
        # use the units in which params was defined
        lamb = self.manager.iu_energy(params["reorg"],
                                      units=self.energy_units)
        kBT = kB_intK*temperature
        time = self.axis.data

        if values is not None:
            cfce = values
        else:
            cfce = 2.0*lamb*kBT*(numpy.exp(-time/ctime)
                                 - 1.0j*(lamb/ctime)*numpy.exp(-time/ctime))

        self._make_me(self.axis, cfce)

        self.lamb = lamb
        self.temperature = temperature
        self.cutoff_time = 5.0*ctime

    def _make_underdamped_brownian(self, params, values=None):
        from .spectraldensities import SpectralDensity
        
        temperature = params["T"]
        ctime = params["gamma"]
        #omega = params["freq"]
        
        # use the units in which params was defined
        lamb = self.manager.iu_energy(params["reorg"],
                                      units=self.energy_units)
        #kBT = kB_intK*temperature
        time = self.axis #.data

        if values is not None:
            cfce = values
        else:
            # Make it via SpectralDensity
            fa = SpectralDensity(time, params)
            
            cf = fa.get_CorrelationFunction(temperature=temperature)
            
            cfce = cf.data
                   #2.0*lamb*kBT*(numpy.exp(-time/ctime)
                   #              - 1.0j*(lamb/ctime)*numpy.exp(-time/ctime))

        self._make_me(self.axis, cfce)

        self.lamb = lamb
        self.temperature = temperature
        self.cutoff_time = 5.0/ctime  
        
        
    def _make_value_defined(self, params, values):
        
        lamb = self.convert_energy_2_internal_u(params["reorg"])
        temp = params["T"]
        
        if "cutoff-time" in params.keys():
            tcut = params["cutoff-time"]
        else:
            tcut = self.axis.max

        if values is not None:
            if len(values) == self.axis.length:
                self.data = self.convert_energy_2_internal_u(values)
            else:
                raise Exception("Incompatible values")
        else:
            raise Exception("Valued-defined correlation function without values")
 
        self.lamb = lamb
        self.temperature = temp
        self.cutoff_time = tcut
           
    #
    # Aritmetic operations
    #
    
    def __add__(self, other):
        """Addition of two correlation functions
        
        """
        t1 = self.axis
        t2 = other.axis
        if t1 == t2:
            
            params = {}
            params["ftype"] = "Value-defined"
            params["reorg"] = self.lamb
            params["T"] = self.temperature
            params["cutoff-time"] = self.cutoff_time
                      
            f = CorrelationFunction(t1, params=params, values=self.data)
            f.add_to_data(other)
            
        else:
            raise Exception("In addition, functions have to share"
                            +" the same TimeAxis object")
            
        return f
    
    def __iadd__(self, other):
        """Inplace addition of two correlation functions
        
        """  
        self.add_to_data(other)
        return self
    
            
    def add_to_data(self, other):
        """Addition of data from a specified CorrelationFunction to this object
        
        """
        t1 = self.axis
        t2 = other.axis
        if t1 == t2:
            
            self.data += other.data
            self.lamb += other.lamb  # reorganization energy is additive
            # cutoff time is take as the longer one of the two
            self.cutoff_time = max(self.cutoff_time, other.cutoff_time)
            
            # FIXME: storage of the component params etc. has to be handled
            self.params = dict(ftype="Value-defined")
            self.params["reorg"] = self.lamb
            self.params["cutoff-time"] = self.cutoff_time
            self.params["T"] = self.temperature
            self._is_composed = True
            self._is_empty = False
            

        else:
            raise Exception("In addition, functions have to share"
                            +" the same TimeAxis object")
        
    def reorganization_energy_consistent(self, rtol=1.0e-3):
        """Checks if the reorganization energy is consistent with the data
        
        Calculates reorganization energy from the data and checks if it
        is within specified tolerance from the expected value
        """
        
        lamb1 = self.measure_reorganization_energy()
        lamb2 = self.convert_energy_2_current_u(self.lamb)
        if (abs(lamb1 - lamb2)/(lamb1+lamb2)) < rtol:
            return True
        else:
            return False

    def is_analytical(self):
        """Returns `True` if analytical

        Returns `True` if the CorrelationFunction object is constructed
        by analytical formula. Returns `False` if the object was constructed
        by numerical transformation from spectral density.
        """

        return bool(self.params["ftype"] in self.analytical_types)


    def get_temperature(self):
        """Returns the temperature of the correlation function

        """
        return self.temperature

    def get_reorganization_energy(self):
        """Returns the reorganization energy of the correlation function

        """
        return self.convert_energy_2_current_u(self.lamb)

    def measure_reorganization_energy(self):
        """Calculates the reorganization energy of the correlation function

        Calculates the reorganization energy of the correlation function by
        integrating its imaginary part.

        """
        #with energy_units("int"):
        primitive = c2h(self.axis, self.data)
        lamb = -numpy.imag(primitive[self.axis.length-1])
        return self.convert_energy_2_current_u(lamb)


    def copy(self):
        """Creates a copy of the current correlation function

        """
        with energy_units(self.energy_units):
            cfce = CorrelationFunction(self.axis, self.params)
        return cfce


    def get_SpectralDensity(self):
        """ Returns a corresponding SpectralDensity object

        Returns a SpectralDensity corresponding to this CorrelationFunction


        """

        from .spectraldensities import SpectralDensity

        # protect this from external units
        with energy_units("int"):
            frequencies = self.axis.get_FrequencyAxis()
            vals = self.get_OddFTCorrelationFunction().data

        # params are saved in user defined units
        with energy_units(self.energy_units):
            # FIXME: This has to be done numerically
            spectd = SpectralDensity(frequencies, self.params, values=vals)

        return spectd


    def get_FTCorrelationFunction(self):
        """Returns a Fourier transform of the correlation function

        Returns a Fourier transform of the correlation function in form
        of an instance of a special class ``FTCorrelationFunction``

        """
        with energy_units(self.energy_units):
            ftcf = FTCorrelationFunction(self.axis, self.params)
        return ftcf

    def get_OddFTCorrelationFunction(self):
        """Returns a odd part of the Fourier transform of correlation function

        Returns the odd part of a Fourier transform of the correlation
        function in form of an instance of a special class
        ``OddFTCorrelationFunction``

        """
        with energy_units(self.energy_units):
            if self.params["ftype"] == "Value-defined":
                oftcf = OddFTCorrelationFunction(self.axis, self.params, 
                                                 values=self.data)
            else:
                oftcf = OddFTCorrelationFunction(self.axis, self.params)
        return oftcf

    def get_EvenFTCorrelationFunction(self):
        """Returns a even part of the Fourier transform of correlation function

        Returns the even part of a Fourier transform of the correlation
        function in form of an instance of a special class
        ``EvenFTCorrelationFunction``

        """
        with energy_units(self.energy_units):
            eftcf = EvenFTCorrelationFunction(self.axis, self.params)
        return eftcf


class FTCorrelationFunction(DFunction, UnitsManaged):
    """Fourier transform of the correlation function

    Numerically calculated Fourier transform of the correlation function

    Parameters
    ----------

    axis: TimeAxis
        Time interval from which the frequency interval is calculated

    params: dictionary
        Dictionary of the correlation function parameters

    """

    def __init__(self, axis, params, values=None):
        super().__init__()

        if not isinstance(axis, FrequencyAxis):
            faxis = axis.get_FrequencyAxis()
            self.axis = faxis
        else:
            self.axis = axis
            

        ftype = params["ftype"]
        if ftype in CorrelationFunction.allowed_types:
            self.ftype = ftype
        else:
            raise Exception("Unknown Correlation Function Type")

        #if (values is not None) and (ftype != "Value-defined"):
        #    raise Exception("Only value defined ftype can have values argument")

        # we need to save the defining energy units
        self.energy_units = self.manager.get_current_units("energy")

        self.params = params

        if values is None:
            # We create CorrelationFunction and FTT it
            with energy_units(self.energy_units):
                cfce = CorrelationFunction(axis, params)
            # data have to be protected from change of units
            with energy_units("int"):
                ftvals = cfce.get_Fourier_transform()
                self.data = ftvals.data
            #self.axis = ftvals.axis
        else:
            # This is not protected from change of units!!!!
            self.data = values
            #self.axis = cfce.axis.get_FrequencyAxis()





class OddFTCorrelationFunction(DFunction, UnitsManaged):
    """Odd part of the Fourier transform of the correlation function

    Numerically calculated odd part Fourier transform of the correlation
    function. Calculated as  Fourier transform of the imaginary part of the
    correlation function.

    Parameters
    ----------

    axis: TimeAxis
        Time interval from which the frequency interval is calculated

    params: dictionary
        Dictionary of the correlation function parameter


    Examples
    --------

    >>> ta = TimeAxis(0.0,1000,1.0)
    >>> params = dict(ftype="OverdampedBrownian",reorg=20,cortime=100,T=300)
    >>> with energy_units("1/cm"):
    ...    ocf = OddFTCorrelationFunction(ta,params)
    ...    print(numpy.allclose(ocf.at(-100), -ocf.at(100)))
    True

    """

    def __init__(self, axis, params, values=None):
        super().__init__()

        if not isinstance(axis, FrequencyAxis):
            faxis = axis.get_FrequencyAxis()
            self.axis = faxis
        else:
            self.axis = axis     

        ftype = params["ftype"]
        if ftype in CorrelationFunction.allowed_types:
            self.ftype = ftype
        else:
            raise Exception("Unknown Correlation Function Type")

        # we need to save the defining energy units
        self.energy_units = self.manager.get_current_units("energy")
        self.params = params
            
        # We create CorrelationFunction and FTT it
        with energy_units(self.energy_units):
            if params["ftype"] == "Value-defined":
                if values is None:
                    raise Exception()
                else:
                    cfce = CorrelationFunction(axis, params, values=values)
            else:
                cfce = CorrelationFunction(axis, params)

        cfce.data = 1j*numpy.imag(cfce.data)

        # data have to be protected from change of units
        with energy_units("int"):
            ftvals = cfce.get_Fourier_transform()
            self.data = numpy.real(ftvals.data)

        self.axis = ftvals.axis



class EvenFTCorrelationFunction(DFunction, UnitsManaged):
    """Even part of the Fourier transform of the correlation function

    Numerically calculated even part Fourier transform of the correlation
    function. Calculated as  Fourier transform of the real part of the
    correlation function.

    Parameters
    ----------

    axis: TimeAxis
        Time interval from which the frequency interval is calculated

    params: dictionary
        Dictionary of the correlation function parameter

    Examples
    --------

    >>> ta = TimeAxis(0.0,1000,1.0)
    >>> params = dict(ftype="OverdampedBrownian",reorg=20,cortime=100,T=300)
    >>> with energy_units("1/cm"):
    ...    ecf = EvenFTCorrelationFunction(ta,params)
    ...    print(numpy.allclose(ecf.at(-100), ecf.at(100)))
    True

    """

    def __init__(self, axis, params):
        super().__init__()

        try:
            ftype = params["ftype"]
            if ftype in CorrelationFunction.allowed_types:
                self.ftype = ftype
            else:
                raise Exception("Unknown Correlation Function Type")

            # we need to save the defining energy units
            self.energy_units = self.manager.get_current_units("energy")


        except:
            raise Exception

        # We create CorrelationFunction and FTT it
        with energy_units(self.energy_units):
            cfce = CorrelationFunction(axis, params)

        cfce.data = numpy.real(cfce.data)

        self.params = params

        # data have to be protected from change of units
        with energy_units("int"):
            ftvals = cfce.get_Fourier_transform()
            self.data = numpy.real(ftvals.data)

        self.axis = ftvals.axis


#FIXME: these functions can go to DFunction
def c2g(timeaxis, coft):
    """ Converts correlation function to lineshape function

    Explicit numerical double integration of the correlation
    function to form a lineshape function.

    Parameters
    ----------

    timeaxis : cu.oqs.time.TimeAxis
        TimeAxis of the correlation function

    coft : complex numpy array
        Values of correlation function given at points specified
        in the TimeAxis object


    """

    time = timeaxis
    preal = numpy.real(coft)
    pimag = numpy.imag(coft)
    splr = interp.UnivariateSpline(time.data,
                                   preal, s=0).antiderivative()(time.data)
    splr = interp.UnivariateSpline(time.data,
                                   splr, s=0).antiderivative()(time.data)
    spli = interp.UnivariateSpline(time.data,
                                   pimag, s=0).antiderivative()(time.data)
    spli = interp.UnivariateSpline(time.data,
                                   spli, s=0).antiderivative()(time.data)
    goft = splr + 1j*spli
    return goft


def c2h(timeaxis, coft):
    """ Integrates correlation function in time with an open upper limit

    Explicit numerical integration of the correlation
    function to form a precursor to the lineshape function.

    Parameters
    ----------

    timeaxis : TimeAxis
        TimeAxis of the correlation function

    coft : complex numpy array
        Values of correlation function given at points specified
        in the TimeAxis object
    """

    time = timeaxis
    preal = numpy.real(coft)
    pimag = numpy.imag(coft)
    splr = interp.UnivariateSpline(time.data,
                                   preal, s=0).antiderivative()(time.data)
    spli = interp.UnivariateSpline(time.data,
                                   pimag, s=0).antiderivative()(time.data)
    hoft = splr + 1j*spli

    return hoft

def h2g(timeaxis, coft):
    """ Integrates and integrated correlation function

    Explicit numerical integration of the correlation
    function to form a precursor to the lineshape function.

    Parameters
    ----------

    timeaxis : TimeAxis
        TimeAxis of the correlation function

    coft : complex numpy array
        Values of correlation function given at points specified
        in the TimeAxis object
    """
    return c2h(timeaxis, coft)
