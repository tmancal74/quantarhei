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
from .correlationfunctions import CorrelationFunction
from .correlationfunctions import FTCorrelationFunction

class SpectralDensity(DFunction, UnitsManaged):
    """This class represents the so-called spectral density

    Parameters
    ----------

    axis : TimeAxis, FrequencyAxis
        ValueAxis object specifying the frequency range directly or through
        Fourier transform frequencies corresponding to a TimeAxis

    params : dictionary
        Parameters of the spectral density



    Examples
    --------

    `SpectralDensity` object can be ctreated with the same parameters as
    `CorrelationFunction`. The temperature can be set, but it is not
    a compulsory parameter.

    >>> from quantarhei import TimeAxis
    >>> params = dict(ftype="OverdampedBrownian", cortime=100, reorg=20, T=300)
    >>> time = TimeAxis(0.0,1000,0.1)
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

    """

    def __init__(self, axis, params):
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

            ctime = params["cortime"]
            lamb = self.manager.iu_energy(params["reorg"],
                                          units=self.energy_units)

            # Masubara frequencies are needed only in time domain
            #try:
            #    N = params["matsubara"]
            #except:
            #    N = 10

            # Temperature is also not needed here, but it is stored
            #kBT = kB_intK*T

            # protect calculation from units management
            with energy_units("int"):
                omega = self.axis.data
                cfce = (2.0*lamb/ctime)*omega/(omega**2 + (1.0/ctime)**2)

            self._make_me(self.axis, cfce)

            # this is in internal units
            self.lamb = lamb


        else:
            raise Exception("Unknown correlation function type of"+
                            "type domain combination.")


    def get_temperature(self):
        """Returns the temperature of the correlation function

        """
        if "T" in self.params.keys():
            return self.temperature
        else:
            raise Exception("SpectralDensity was not set with temperature")


    def copy(self):
        """Creates a copy of the current correlation function"""
        return SpectralDensity(self.axis, self.params)


    def get_CorrelationFunction(self, temperature=None):
        """Returns correlation function corresponding to the spectral density

        """
        params = self.params.copy()
        if temperature is not None:
            params["T"] = temperature

        time = self.axis.get_TimeAxis()
        return CorrelationFunction(time, params)
        
        
    def get_FTCorrelationFunction(self,temperature=None):
        """Returns Fourier transformed correlation function
        
        Fourier transformed correlation function is calculated from the
        analytical formula connecting spectral density and FT correlation
        function.
        
        Parameters
        ----------
        
        temparature : optional
            Temperature which can be missing among the spectral density 
            parameters
            
            
        """
        params = self.params.copy()
        if temperature is not None:
            params["T"] = temperature
            
        ind_of_zero, diff = self.axis.locate(0.0)
        atol = 1.0e-7
        with energy_units("int"): #self.energy_units):
            if numpy.abs(diff) < atol:
                vals = (1.0 + (1.0/numpy.tanh(self.axis.data)))*self.data
            else:
                vals = self.data # (1.0 + (1.0/numpy.tanh(self.axis.data)))*self.data
                
        with energy_units(self.energy_units):
            ftc = FTCorrelationFunction(self.axis, params, values=vals)
 
        return ftc
        
        