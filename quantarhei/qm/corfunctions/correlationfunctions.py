# -*- coding: utf-8 -*-
"""
    Provides typical Bath correlation function types.

    Most important types of bath or energy gap correlation functions are
    provided. Where possible, the correlation function is calculated
    from the parameters from analytical formulae. Where such formulae are
    not available, correlation function is calculated by transformation
    of the spectral density.

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

    Details of Classes Provided
    ---------------------------
    
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


    Parameters
    ----------

    axis : TimeAxis
        TimeAxis object specifying the time interval on which the
        correlation function is defined.

    params : dictionary
        A dictionary of the correlation function parameters

    values : optional
        Correlation function can be set by specifying values at all times


    """

    allowed_types = ("OverdampedBrownian-HighTemperature",
                     "OverdampedBrownian",
                     "OverdampedBrownian_from_Specdens",
                     "UnderdampedBrownian",
                     "Underdamped",
                     "B777",
                     "CP29",
                     "Value-defined"
                     )

    analytical_types = ("OverdampedBrownian-HighTemperature",
                        "OverdampedBrownian")
    
    energy_params = ("reorg", "omega", "freq", "fcp", "g_FWHM", "l_FWHM",\
                     "freq1", "freq2", "gamma")

    def __init__(self, axis=None, params=None , values=None):
        super().__init__()
        
        if (axis is not None) and (params is not None):
            
            # FIXME: values might also need handling according to specified 
            # energy units
    
            # handle axis (it can be TimeAxis or FrequencyAxis)
            if not isinstance(axis, TimeAxis):
                taxis = axis.get_TimeAxis()
                self.axis = taxis
            else:
                self.axis = axis
    
            # handle params
            self.params = []  # this will always be a list of components
            p2calc = []
            try:
                # if this passes, we assume params is a dictionary
                params.keys()
                self._is_composed = False
                p2calc.append(params)
                
            except:
                # othewise we assume it is a list of dictionaries 
                self._is_composed = True
                for p in params:
                    p2calc.append(p)
                
                
            self.lamb = 0.0
            self.temperature = -1.0
            self.cutoff_time = 0.0

            for params in p2calc:
                
                try:
                    ftype = params["ftype"]
                    
                    if ftype not in CorrelationFunction.allowed_types:
                        raise Exception("Unknown CorrelationFunction type")
        
                    # we mutate the parameters into internal units
                    prms = {}
                    for key in params.keys():
                        if key in self.energy_params:
                            prms[key] = self.convert_energy_2_internal_u(params[key])
                        else:
                            prms[key] = params[key]
                            
                except:
                    raise Exception("Dictionary of parameters does not contain "
                                    +" `ftype` key")
                    
                self.params.append(prms)                    
    
    
            if values is None:
                #
                # loop over parameter sets
                #
                for prms in self.params:
                    
#                    try:
#                        ftype = params["ftype"]
#                        
#                        if ftype not in CorrelationFunction.allowed_types:
#                            raise Exception("Unknown CorrelationFunction type")
#            
#                        # we mutate the parameters into internal units
#                        prms = {}
#                        for key in params.keys():
#                            if key in self.energy_params:
#                                prms[key] = self.convert_energy_2_internal_u(params[key])
#                            else:
#                                prms[key] = params[key]
#                                
#                    except:
#                        raise Exception("Dictionary of parameters does not contain "
#                                        +" `ftype` key")                    
        
                    if ftype == "OverdampedBrownian-HighTemperature":
            
                        self._make_overdamped_brownian_ht(prms) #, values=values)
            
                    elif ftype == "OverdampedBrownian":
            
                        self._make_overdamped_brownian(prms) #, values=values)
                        
                    elif ftype == "UnderdampedBrownian":
                        
                        self._make_underdamped_brownian(prms) #, values=values)
                        
                    elif ftype == "Underdamped":
                        
                        self._make_underdamped(params, values=values)
                        
                    elif ftype == "B777":
                        
                        self._make_B777(params, values=values)
                        
                    elif ftype == "CP29":
                        
                        self._make_CP29_spectral_density(params, values=values)
            
                    elif ftype == "Value-defined":
            
                        self._make_value_defined(prms, values)
            
                    else:
                        raise Exception("Unknown correlation function type or"+
                                        "type domain combination.")
                        
                    #self.params.append(prms)
                    
            else:
                
                self._add_me(self.axis, values) 
                # update reorganization energy
                self.lamb = 0.0
                self.temperature = self.params[0]["T"]
                for prms in self.params:
                    self.lamb += prms["reorg"]
                    if self.temperature != prms["T"]:
                        raise Exception("Inconsistent temperature! "
                                        +"Temperatures of all "
                                        +"components have to be the same")
                    
                #FIXME: set cut-off time and temperature
                #self._set_temperature_and_cutoff_time(self.params[0])
                    
                

    def _matsubara(self, kBT, ctime, nof):
        """Matsubara frequency part of the Brownian correlation function
        
        """
        msf = 0.0
        nut = 2.0*numpy.pi*kBT
        time = self.axis.data
        for i in range(0, nof):
            n = i+1
            msf += nut*n*numpy.exp(-nut*n*time)/((nut*n)**2-(1.0/ctime)**2)
        return msf

    def _set_temperature_and_cutoff_time(self, temperature, ctime):
        """Sets the temperature and cutoff time of for the component
        
        """
        
        # Temperatures of all components have to be the same
        # is this the first time that temperature is assigned?
        if self.temperature == -1.0: 
            self.temperature = temperature
        elif self.temperature != temperature:
            raise Exception("Inconsistent temperature! Temperatures of all "
                            +"components have to be the same")
            
        # longest cortime has to be preserved
        new_cutoff_time = 5.0*ctime
        if new_cutoff_time > self.cutoff_time: 
            self.cutoff_time = new_cutoff_time
            
            
    def _make_overdamped_brownian(self, params): #, values=None):
        """Creates the overdamped Brownian oscillator component
        of the correlation function
        
        """
        
        temperature = params["T"]
        ctime = params["cortime"]
        lamb = params["reorg"]

        if "matsubara" in params.keys():
            nmatsu = params["matsubara"]
        else:
            nmatsu = 10

        kBT = kB_intK*temperature
        time = self.axis.data

        #if values is not None:
        #    cfce = values
            
        #else:
        if True:
            
            cfce = (lamb/(ctime*numpy.tan(1.0/(2.0*kBT*ctime))))\
                *numpy.exp(-time/ctime) \
                - 1.0j*(lamb/ctime)*numpy.exp(-time/ctime)

            cfce += (4.0*lamb*kBT/ctime) \
                *self._matsubara(kBT, ctime, nmatsu)
        
        # this is a call to the function inherited from DFunction class 
        self._add_me(self.axis, cfce)
        
        # update reorganization energy
        self.lamb += lamb
        # check temperature and update cutoff time
        self._set_temperature_and_cutoff_time(temperature, 5.0*ctime)      



    def _make_overdamped_brownian_ht(self, params): # , values=None):
        """Creates the high temperature overdamped Brownian oscillator 
        component of the correlation function
        
        """
        temperature = params["T"]
        ctime = params["cortime"]
        lamb = params["reorg"]
        
        kBT = kB_intK*temperature
        time = self.axis.data

        #if values is not None:
        #    cfce = values
        #    
        #else:
        if True:
            cfce = 2.0*lamb*kBT*(numpy.exp(-time/ctime)
                                 - 1.0j*(lamb/ctime)*numpy.exp(-time/ctime))

        # this is a call to the function inherited from DFunction class 
        self._add_me(self.axis, cfce)

        # update reorganization energy
        self.lamb += lamb
        
        # check temperature and update cutoff time
        self._set_temperature_and_cutoff_time(temperature, 5.0*ctime)  
        

    def _make_underdamped_brownian(self, params): #, values=None):
        """Creates underdamped Brownian oscillator component of the correlation
        function
        
        
        """
        from .spectraldensities import SpectralDensity
        
        temperature = params["T"]
        ctime = params["gamma"]
        #omega = params["freq"]
        lamb = params["reorg"]
        
        #kBT = kB_intK*temperature
        time = self.axis #.data

        #if values is not None:
        #    cfce = values
        #    
        #else:
        if True:
            with energy_units("int"):
                # Make it via SpectralDensity
                fa = SpectralDensity(time, params)
            
                cf = fa.get_CorrelationFunction(temperature=temperature)
            
                cfce = cf.data
                   #2.0*lamb*kBT*(numpy.exp(-time/ctime)
                   #              - 1.0j*(lamb/ctime)*numpy.exp(-time/ctime))

        # this is a call to the function inherited from DFunction class 
        self._add_me(self.axis, cfce)

        # update reorganization energy
        self.lamb += lamb
        
        # check temperature and update cutoff time
        self._set_temperature_and_cutoff_time(temperature, 5.0/ctime)  
        
        
    def _make_underdamped(self, params, values=None):
        from .spectraldensities import SpectralDensity
        
        temperature = params["T"]
        ctime = params["gamma"]
        
        # use the units in which params was defined
        lamb = params["reorg"]
        time = self.axis #.data

        if values is not None:
            cfce = values
        else:
            # Make it via SpectralDensity
            fa = SpectralDensity(time, params)
            
            cf = fa.get_CorrelationFunction(temperature=temperature)
            
            cfce = cf.data

         # this is a call to the function inherited from DFunction class 
        self._add_me(self.axis, cfce)

        # update reorganization energy
        self.lamb += lamb
        
        # check temperature and update cutoff time
        self._set_temperature_and_cutoff_time(temperature, 5.0*ctime) 
        
    def _make_B777(self, params, values=None):
        from .spectraldensities import SpectralDensity
        
        temperature = params["T"]
        ctime = params["gamma"]
        
        # use the units in which params was defined
        lamb = self.manager.iu_energy(params["reorg"],
                                      units=self.energy_units)
        time = self.axis #.data

        if values is not None:
            cfce = values
        else:
            # Make it via SpectralDensity
            fa = SpectralDensity(time, params)
            
            cf = fa.get_CorrelationFunction(temperature=temperature)
            
            cfce = cf.data
            
        # this is a call to the function inherited from DFunction class 
        self._add_me(self.axis, cfce)

        # update reorganization energy
        self.lamb += lamb
        
        # check temperature and update cutoff time
        self._set_temperature_and_cutoff_time(temperature, 5.0*ctime)     

        
    def _make_CP29_spectral_density(self, params, values=None):
        from .spectraldensities import SpectralDensity
        
        temperature = params["T"]
        ctime = params["gamma"]
        #omega = params["freq"]
        
        # use the units in which params was defined
        lamb = self.manager.iu_energy(params["reorg"],
                                      units=self.energy_units)
        print('correlation function lamb in int units %f' %lamb)
        time = self.axis #.data

        if values is not None:
            cfce = values
        else:
            # Make it via SpectralDensity
            fa = SpectralDensity(time, params)
            

            cf = fa.get_CorrelationFunction(temperature=temperature)
            
            cfce = cf.data

        # this is a call to the function inherited from DFunction class 
        self._add_me(self.axis, cfce)

        # update reorganization energy
        self.lamb += lamb
        
        # check temperature and update cutoff time
        self._set_temperature_and_cutoff_time(temperature, 5.0*ctime)   
        
    def _make_value_defined(self, params, values):
        
        lamb = params["reorg"]
        temperature = params["T"]
        
        if "cutoff-time" in params.keys():
            ctime = params["cutoff-time"]
        else:
            ctime = self.axis.max

        if values is not None:
            if len(values) == self.axis.length:
                cfce = values
            else:
                raise Exception("Incompatible values")
        else:
            raise Exception("Valued-defined correlation function without values")
 
        # this is a call to the function inherited from DFunction class 
        self._add_me(self.axis, cfce)

        # update reorganization energy
        self.lamb += lamb
        
        # check temperature and update cutoff time
        self._set_temperature_and_cutoff_time(temperature, ctime)  
           
    #
    # Aritmetic operations
    #
    
    def __add__(self, other):
        """Addition of two correlation functions
        
        """
        t1 = self.axis
        t2 = other.axis
        if t1 == t2:
                      
            f = CorrelationFunction(t1, params=self.params)
            f.add_to_data(other)
            
        else:
            raise Exception("In addition, functions have to share"
                            +" the same TimeAxis object")
            
        return f
    
    def __iadd__(self, other):
        """Inplace addition of two correlation functions
        
        """  
        self.add_to_data2(other)
        return self
    
            
    def add_to_data(self, other):
        """Addition of data from a specified CorrelationFunction to this object
        
        """
        t1 = self.axis
        t2 = other.axis
        if t1 == t2:
            
            self.data += other.data
            self.lamb += other.lamb  # reorganization energy is additive
            if other.cutoff_time > self.cutoff_time: 
                self.cutoff_time = other.cutoff_time  
                
            if self.temperature != other.temperature:
                raise Exception("Cannot add two correlation functions on different temperatures")
    
            for p in other.params:
                self.params.append(p)
                
            self._is_composed = True
            self._is_empty = False
            

        else:
            raise Exception("In addition, functions have to share"
                            +" the same TimeAxis object")
 
    def add_to_data2(self, other):
        """Addition of data from a specified CorrelationFunction to this object
        
        """
        if self == other:
            ocor = CorrelationFunction(other.axis,other.params)
        else:
            ocor = other
            
        t1 = self.axis
        t2 = ocor.axis
        if t1 == t2:
            
            self.data += ocor.data
            self.lamb += ocor.lamb  # reorganization energy is additive
            if ocor.cutoff_time > self.cutoff_time: 
                self.cutoff_time = ocor.cutoff_time  
                
            if self.temperature != ocor.temperature:
                raise Exception("Cannot add two correlation functions on different temperatures")
    

            for p in ocor.params:
                self.params.append(p)
            
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

    def get_correlation_time(self):
        """Returns correlation time associated with the first component 
        of the bath correlation function
        """
        return self.params[0]["ctime"]

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
        #with energy_units(self.energy_units):
        cfce = CorrelationFunction(self.axis, self.params)
        return cfce


    def get_SpectralDensity(self, fa=None):
        """
        Returns a SpectralDensity corresponding to this CorrelationFunction.
        If a FrequencyAxis object is included, the SpectralDensity
        object will be returned with that FrequencyAxis instance as its 
        frequency axis.

        """

        from .spectraldensities import SpectralDensity

        # protect this from external units
        with energy_units("int"):
            frequencies = self.axis.get_FrequencyAxis()
            vals = self.get_OddFTCorrelationFunction().data
            
            if fa is not None:
                if numpy.all(numpy.isclose(fa.data, frequencies.data, 1e-5)):
                    time = ta
                else:
                    raise Exception("The provided FrequencyAxis does not "
                                    + "have the same data as the Fourier "
                                    + "transformed axis")

            # FIXME: how to set the limit of SpectralDensity at w->0
            spectd = SpectralDensity(frequencies, self.params, values=vals)

        return spectd


    def get_FTCorrelationFunction(self):
        """Returns a Fourier transform of the correlation function

        Returns a Fourier transform of the correlation function in form
        of an instance of a special class ``FTCorrelationFunction``

        """
        with energy_units("int"):
            ftcf = FTCorrelationFunction(self.axis, self.params)
        return ftcf


    def get_OddFTCorrelationFunction(self):
        """Returns a odd part of the Fourier transform of correlation function

        Returns the odd part of a Fourier transform of the correlation
        function in form of an instance of a special class
        ``OddFTCorrelationFunction``

        """

        with energy_units("int"):
            oftcf = OddFTCorrelationFunction(self.axis, self.params)
        return oftcf


    def get_EvenFTCorrelationFunction(self):
        """Returns a even part of the Fourier transform of correlation function

        Returns the even part of a Fourier transform of the correlation
        function in form of an instance of a special class
        ``EvenFTCorrelationFunction``

        """
        with energy_units("int"):
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
    energy_params = ("reorg", "omega", "freq")
    
    def __init__(self, axis, params, values=None):
        super().__init__()

        if not isinstance(axis, FrequencyAxis):
            faxis = axis.get_FrequencyAxis()
            self.axis = faxis
        else:
            self.axis = axis
            
        # handle params
        self.params = []  # this will always be a list of components
        p2calc = []
        try:
            # if this passes, we assume params is a dictionary
            params.keys()
            self._is_composed = False
            p2calc.append(params)
            
        except:
            # othewise we assume it is a list of dictionaries 
            self._is_composed = True
            for p in params:
                p2calc.append(p)
                
        for params in p2calc:
            
            try:
                ftype = params["ftype"]
                
                if ftype not in CorrelationFunction.allowed_types:
                    raise Exception("Unknown CorrelationFunction type")
    
                # we mutate the parameters into internal units
                prms = {}
                for key in params.keys():
                    if key in self.energy_params:
                        prms[key] = self.convert_energy_2_internal_u(params[key])
                    else:
                        prms[key] = params[key]
                    
            except:
                raise Exception("Dictionary of parameters does not contain "
                                +" `ftype` key")
        self.params = prms

        if values is None:
            # data have to be protected from change of units
            with energy_units("int"):
                cfce = CorrelationFunction(axis, self.params)
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

        # handle params
        self.params = []  # this will always be a list of components
        p2calc = []
        try:
            # if this passes, we assume params is a dictionary
            params.keys()
            self._is_composed = False
            p2calc.append(params)
            
        except:
            # othewise we assume it is a list of dictionaries 
            self._is_composed = True
            for p in params:
                p2calc.append(p)
        
        for params in p2calc:
        
            ftype = params["ftype"]
            if ftype not in CorrelationFunction.allowed_types:
                raise Exception("Unknown Correlation Function Type")
    
            self.params.append(params)
                
            # We create CorrelationFunction and FTT it
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
                ndata = numpy.real(ftvals.data)

            self._add_me(self.axis,ndata)


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

    def __init__(self, axis, params, values=None):
        super().__init__()

        if not isinstance(axis, FrequencyAxis):
            faxis = axis.get_FrequencyAxis()
            self.axis = faxis
        else:
            self.axis = axis 
            
        # handle params
        self.params = []  # this will always be a list of components
        p2calc = []
        try:
            # if this passes, we assume params is a dictionary
            params.keys()
            self._is_composed = False
            p2calc.append(params)
            
        except:
            # othewise we assume it is a list of dictionaries 
            self._is_composed = True
            for p in params:
                p2calc.append(p)
        
        for params in p2calc:
            
            ftype = params["ftype"]
            if ftype not in CorrelationFunction.allowed_types:
                raise Exception("Unknown Correlation Function Type: "+ftype)


            self.params.append(params)
            
            # We create CorrelationFunction and FTT it
            if params["ftype"] == "Value-defined":
                if values is None:
                    raise Exception()
                else:
                    cfce = CorrelationFunction(axis, params, values=values)
            else:
                cfce = CorrelationFunction(axis, params)
                
            cfce.data = numpy.real(cfce.data)
    
            # data have to be protected from change of units
            with energy_units("int"):
                ftvals = cfce.get_Fourier_transform()
                ndata = numpy.real(ftvals.data)

            self._add_me(self.axis,ndata)


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
