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


class CorrelationFunction(DFunction, UnitsManaged):
    """Provides typical Bath correlation functions

    Parameters
    ----------

    timeAxis : TimeAxis
        TimeAxis object specifying the time interval on which the
        correlation function is defined.


    Types of correlation function provided
    --------------------------------------
    OverdampedBrownian-HighTemperature :
        OverdampedBrownian oscillator in high temperature limit

    OverdampedBrownian :
        General overdampedBrownian oscillator

    Examples
    --------

    >>> from quantarhei import TimeAxis
    >>> params = dict(ftype="OverdampedBrownian", cortime=100, reorg=20, \
                      T=300)
    >>> time = TimeAxis(0.0,1000,0.1)
    >>> with energy_units("1/cm"): \
            cf = CorrelationFunction(time,params)



    """


    allowed_types = ("OverdampedBrownian-HighTemperature",
                     "OverdampedBrownian", "Value-defined")

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

            self._make_overdumped_brownian_ht(params)

        elif self.ftype == "OverdampedBrownian":

            self._make_overdumped_brownian(params)

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

    def _make_overdumped_brownian(self, params):

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

        cfce = (lamb/(ctime*numpy.tan(1.0/(2.0*kBT*ctime))))\
            *numpy.exp(-time/ctime) \
            - 1.0j*(lamb/ctime)*numpy.exp(-time/ctime)

        cfce += (4.0*lamb*kBT/ctime) \
            *self._matsubara(kBT, ctime, nmatsu)

        self._make_me(self.axis, cfce)
        self.lamb = lamb
        self.temperature = temperature

        self.cutoff_time = 5.0*ctime


    def _make_overdumped_brownian_ht(self, params):
        temperature = params["T"]
        ctime = params["cortime"]
        # use the units in which params was defined
        lamb = self.manager.iu_energy(params["reorg"],
                                      units=self.energy_units)
        kBT = kB_intK*temperature
        time = self.axis.data

        cfce = 2.0*lamb*kBT*(numpy.exp(-time/ctime)
                             - 1.0j*(lamb/ctime)*numpy.exp(-time/ctime))

        self._make_me(self.axis, cfce)

        self.lamb = lamb
        self.temperature = temperature
        self.cutoff_time = 5.0*ctime

    def _make_value_defined(self, params, values):
        lamb = self.manager.iu_energy(params["reorg"],
                                      units=self.energy_units)
        self.lamb = lamb
        if "cutoff-time" in params.keys():
            tcut = params["cutoff-time"]
        else:
            tcut = self.axis.tmax
        self.cutoff_time = tcut

        if len(values) == self.axis.length:
            self.data = values
        else:
            raise Exception("Incompatible values")

    def get_temperature(self):
        """Returns the temperature of the correlation function

        """
        return self.temperature


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

        # params are saved in user defined units
        with energy_units(self.energy_units):
            spectd = SpectralDensity(frequencies, self.params)

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

    """

    def __init__(self, axis, params, values=None):
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

        self.params = params

        if values is None:
            # data have to be protected from change of units
            with energy_units("int"):
                ftvals = cfce.get_Fourier_transform()
                self.data = ftvals.data
            self.axis = ftvals.axis
        else:
            # This is not protected from change of units!!!!
            self.data = values
            self.axis = cfce.axis.get_FrequencyAxis()
            
            
        


class OddFTCorrelationFunction(DFunction, UnitsManaged):
    """Odd part of the Fourier transform of the correlation function

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

        cfce.data = 1j*numpy.imag(cfce.data)

        self.params = params

        # data have to be protected from change of units
        with energy_units("int"):
            ftvals = cfce.get_Fourier_transform()
            self.data = ftvals.data

        self.axis = ftvals.axis



class EvenFTCorrelationFunction(DFunction, UnitsManaged):
    """Even part of the Fourier transform of the correlation function

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
            self.data = ftvals.data

        self.axis = ftvals.axis


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
