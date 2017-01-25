# -*- coding: utf-8 -*-
"""
    Quantarhei package (http://www.github.com/quantarhei)

    spectraldensities module


"""
import numpy

from ...core.dfunction import DFunction
from ...core.managers import UnitsManaged
from ...core.managers import energy_units
from ...core.time import TimeAxis
#from ...core.frequency import FrequencyAxis
from .correlationfunctions import CorrelationFunction
from .correlationfunctions import FTCorrelationFunction
from ...core.units import kB_int

#from .correlationfunctions import c2h

class SpectralDensity(DFunction, UnitsManaged):
    """This class represents the so-called spectral density

    Parameters
    ----------

    axis : TimeAxis, FrequencyAxis
        ValueAxis object specifying the frequency range directly or through
        Fourier transform frequencies corresponding to a TimeAxis

    params : dictionary
        Parameters of the spectral density


    Methods
    -------

....is_analytical()
        Returns `True` if the spectral density is calculated from an analytical
        formula, `False` otherwise.

    copy()
        Makes a copy of the SpectralDensity object

    get_temperature()
        Returns temperature assigned to the spectral density or raises an
        exception if it was not set

    get_reorganization_energy()
        Returns the reorganization energy parameters of
        the spectral density

    measure_reorganization_energy()
        Calculates reorganization energy from the shape of the spectral
        density

    get_CorrelationFunction()
        Returns numerically calculated correlation function, based on the
        spectral density

    get_FTCorrelationFunction()
        Returns a numerically calculated Fourier transform of the correlation
        function


    Examples
    --------

    `SpectralDensity` object can be ctreated with the same parameters as
    `CorrelationFunction`. The temperature can be set, but it is not
    a compulsory parameter.

    >>> from quantarhei import TimeAxis
    >>> params = dict(ftype="OverdampedBrownian", cortime=100, reorg=20, T=300)
    >>> time = TimeAxis(0.0,1000,1.0)
    >>> with energy_units("1/cm"):\
           sd = SpectralDensity(time, params)

    >>> cf = sd.get_CorrelationFunction()

    >>> print(sd.get_temperature())
    300

    If we create the same without temperature

    >>> parwoT = dict(ftype="OverdampedBrownian", cortime=100, reorg=20)
    >>> with energy_units("1/cm"):\
           sdwoT = SpectralDensity(time, parwoT)

    everything is alright. `CorrelationFunction`, however, cannot be created
    as above, because temperature must be known. Attempt to create `CorrelationFunction`
    as above would lead to an exception. We have to specify temperature as
    a parameter to the `get_CorrelationFunction` method.

    >>> cf = sdwoT.get_CorrelationFunction(temperature=300)

    Reorganization of the spectral density is an input parameter which
    can be obtained by calling the corresponding method

    >>> with energy_units("1/cm"):
    ...     print(sdwoT.get_reorganization_energy())
    20.0

    At the same time, reorganization energy can be calculated from the
    shape of the spectral density by integrating over it. The accuracy
    of such estimation depends on numerics, hence the relative tolerance of
    only 1.0e-2 below

    >>> lamb_definition = sd.get_reorganization_energy()
    >>> lamb_measured = sd.measure_reorganization_energy()
    >>> print(numpy.allclose(lamb_definition, lamb_measured, rtol=1.0e-2))
    True

    """

    analytical_types = ("OverdampedBrownian")

    def __init__(self, axis, params, values=None):
        super().__init__()

        if isinstance(axis, TimeAxis):
            # protect the frequency axis creation from units management
            with energy_units("int"):
                faxis = axis.get_FrequencyAxis()
            self.axis = faxis
        else:
            self.axis = axis

        self._splines_initialized = False

        try:
            ftype = params["ftype"]
            if ftype in CorrelationFunction.allowed_types:
                self.ftype = ftype
            else:
                raise Exception("Unknown Correlation Function Type")

            # we need to save the defining energy units
            self.energy_units = self.manager.get_current_units("energy")
            # because params are in the defining units
            self.params = params

        except:
            raise Exception

        if "T" in params.keys():
            self.temperature = params["T"]

        if self.ftype == "OverdampedBrownian":

            self._make_overdamped_brownian()

        elif self.ftype == "UnderdampedBrownian":

            self._make_underdamped_brownian(params)
            
        elif self.ftype == "Value-defined":

            self._make_value_defined(values=values)

        else:
            raise Exception("Unknown correlation function type or"+
                            " type domain combination.")

    def _make_overdamped_brownian(self, values=None):
        """ Sets the Overdamped Brownian oscillator spectral density

        """
        ctime = self.params["cortime"]
        lamb = self.manager.iu_energy(self.params["reorg"],
                                      units=self.energy_units)

        # protect calculation from units management
        with energy_units("int"):
            omega = self.axis.data
            cfce = (2.0*lamb/ctime)*omega/(omega**2 + (1.0/ctime)**2)

        if values is not None:
            self._make_me(self.axis, values)
        else:
            self._make_me(self.axis, cfce)

        # this is in internal units
        self.lamb = lamb

        self.lim_omega = numpy.zeros(2)
        self.lim_omega[0] = 0.0
        self.lim_omega[1] = ctime

    def _make_underdamped_brownian(self, params, values=None):
         
        #temperature = params["T"]
        ctime = params["gamma"]
        # use the units in which params was defined
        omega0 = self.manager.iu_energy(params["freq"],
                                      units=self.energy_units)
        lamb = self.manager.iu_energy(params["reorg"],
                                      units=self.energy_units)
        #kBT = kB_intK*temperature
        #time = self.axis.data 
        
        # protect calculation from units management
        with energy_units("int"):
            omega = self.axis.data
            cfce = (lamb*ctime)*omega/((omega-omega0)**2 + (ctime)**2) \
                  +(lamb*ctime)*omega/((omega+omega0)**2 + (ctime)**2)

        if values is not None:
            self._make_me(self.axis, values)
        else:
            self._make_me(self.axis, cfce)

        # this is in internal units
        self.lamb = lamb            
        self.lim_omega = numpy.zeros(2)
        self.lim_omega[0] = 0.0
        self.lim_omega[1] = 0.0
        
    def _make_value_defined(self, values=None):
        """ Value defined spectral density

        """
        if values is None:
            raise Exception()
            
        self._make_me(self.axis, values)
        self.lamb = self.manager.iu_energy(self.params["reorg"],
                                      units=self.energy_units)

        self.lim_omega = numpy.zeros(2)
        self.lim_omega[0] = 0.0
        self.lim_omega[1] = 0.0

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
            #params["T"] = self.temperature
            #params["cutoff-time"] = self.cutoff_time
                      
            f = SpectralDensity(t1, params=params, values=self.data)
            f.add_to_data(other)
            for i in range(2):
                self.lim_omega[i] += other.lim_omega[i]
            
        else:
            raise Exception("In addition, functions have to share"
                            +" the same FrequencyAxis object")
            
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
            #self.cutoff_time = max(self.cutoff_time, other.cutoff_time)
            
            # FIXME: storage of the component params etc. has to be handled
            self.params = dict(ftype="Value-defined")
            self.params["reorg"] = self.lamb
            #self.params["cutoff-time"] = self.cutoff_time
            #self.params["T"] = self.temperature
            self._is_composed = True
            self._is_empty = False
            

        else:
            raise Exception("In addition, functions have to share"
                            +" the same TimeAxis object")


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
        if "T" in self.params.keys():
            return self.temperature
        else:
            raise Exception("SpectralDensity was not assigned temperature")

    def get_reorganization_energy(self):
        """Returns the reorganization energy of the cspectral density

        """
        return self.convert_energy_2_current_u(self.lamb)

    def measure_reorganization_energy(self):
        """Calculates the reorganization energy of the spectral density

        Calculates the reorganization energy of the spectral density by
        integrating over frequency.
        """
        import scipy.interpolate as interp

        integr = self.data/self.axis.data
        uvspl = interp.UnivariateSpline(self.axis.data, integr, s=0)
        integ = uvspl.integral(0.0, self.axis.max)/numpy.pi

#        ind_of_zero = self.axis.nearest(0.0)
#        cdouble = self.data[ind_of_zero:self.axis.length]
#        omega = self.axis.data[ind_of_zero:self.axis.length]
#        integrand = cdouble/omega
#        length = len(integrand)
#        faxis = FrequencyAxis(0.0,length,self.axis.step)
#        integ = numpy.real(c2h(faxis,integrand)[length-1])/numpy.pi

        return integ


    def copy(self):
        """Creates a copy of the current correlation function

        """
        return SpectralDensity(self.axis, self.params)


    def get_CorrelationFunction(self, temperature=None):
        """Returns correlation function corresponding to the spectral density

        """

        params = self.params.copy()
        if temperature is not None:
            params["T"] = temperature

        time = self.axis.get_TimeAxis()

        # everything has to be protected from change of units
        with energy_units("int"):

            ftcf = self.get_FTCorrelationFunction(temperature=params["T"])
            cftd = ftcf.get_inverse_Fourier_transform()
            cfce = CorrelationFunction(time, params, values=cftd.data)

        return cfce


    def get_FTCorrelationFunction(self, temperature=None):
        """Returns Fourier transformed correlation function

        Fourier transformed correlation function is calculated from the
        analytical formula connecting spectral density and FT correlation
        function.

        Parameters
        ----------

        temperature : optional
            Temperature which can be missing among the spectral density
            parameters


        """
        params = self.params.copy()
        if temperature is not None:
            params["T"] = temperature

        temp = params["T"]
        #params["ftype"] = "Value-defined"

        ind_of_zero, diff = self.axis.locate(0.0)
        atol = 1.0e-7
        twokbt = 2.0*kB_int*temp

        with energy_units("int"):
            # if zero is sufficiently away from any point that is evaluated
            if numpy.abs(diff) > atol:
                # do the evaluation directly
                vals = (1.0 +
                        (1.0/numpy.tanh(self.axis.data/twokbt)))*self.data
            # otherwise
            else:
                data = self.data
                vals = numpy.zeros(self.data.shape)
                # evaluate everything before zero
                omega = self.axis.data[0:ind_of_zero]
                spect = data[0:ind_of_zero]
                auxi = (1.0 + (1.0/numpy.tanh(omega/twokbt)))*spect
                vals[0:ind_of_zero] = auxi
                # then after zero
                omega = self.axis.data[ind_of_zero+1:self.axis.length]
                spect = data[ind_of_zero+1:self.axis.length]
                auxi = (1.0 + (1.0/numpy.tanh(omega/twokbt)))*spect
                vals[ind_of_zero+1:self.axis.length] = auxi
                # and used the limit at zero devided by omega
                vals[ind_of_zero] = 2.0*self.lamb*twokbt*self.lim_omega[1]


        with energy_units(self.energy_units):
            ftc = FTCorrelationFunction(self.axis, params, values=vals)

        return ftc
