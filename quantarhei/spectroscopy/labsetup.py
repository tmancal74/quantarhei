# -*- coding: utf-8 -*-
"""
    Laboratory set-up for non-linear spectroscopy
    
    

    Class Details
    -------------

"""

import numpy


from ..utils import Integer
from ..utils.vectors import X
from ..core.time import TimeAxis
from ..core.frequency import FrequencyAxis
from ..core.dfunction import DFunction

class LabSetup:
    """
    
    Class representing laboratory setup for non-linear spectroscopic
    experiments. It holds information about pulse shapes and polarizations.
    
    Pulses can be set in time- and/or frequency-domain. Consistency between
    the domains is not checked nor enforced. Consistent conversion between
    domains is provided by convenience routines [TO BE IMPLEMENTED]
    
    """
    
    number_of_pulses = Integer("number_of_pulses")
    
    def __init__(self, nopulses = 3):
        
        self.number_of_pulses = nopulses
    
        self.M4 = numpy.array([[4.0, -1.0, -1.0],
                               [-1.0, 4.0, -1.0],
                               [-1.0,-1.0,  4.0]])/30.0
    
        self.timeaxis = None
        self.freqaxis = None
        
        self.axis_type = None
        
        self.pulse_t = [None]*nopulses
        self.pulse_f = [None]*nopulses
        
        self.has_freqdomain = False
        self.has_timedomain = False
                        

    def set_pulse_shapes(self, axis, params):
        """Sets the pulse properties
        
        
        Pulse shapes or spectra are set in this routine. If axis is of 
        `TimeAxis` type, the parameters are understood as time domain, 
        if axis is of `FrequencyAxis` type, they are understood as frequency
        domain.
        
        
        Parameters
        ----------
        
        axis : TimeAxis or FrequencyAxis
            Quantarhei time axis object, which specifies the values for which 
            pulse properties are defined. If `TimeAxis` is specified, the 
            parameters are understood as time domain, if `FrequencyAxis`
            is specified, they are understood as frequency domain.
            
        params : dict
            Dictionary of pulse parameters. The parameters are the following:
            `ptype` is the pulse type with possible values `Gaussian` and 
            `numeric`. Time domain pulses are specified with their center
            at t = 0.
            
            `Gaussian` pulse has further parameters `amplitude`, `FWHM`,
            and `frequency` with obvious meanings. `FWHM` is speficied in `fs`,
            `frequency` is specified in energy units, while `amplitude`
            is in units of [energy]/[transition dipole moment]. 
            
            `numeric` pulse is specified by a second parameters `function` 
            which should be of DFunction type and specifies line shape around
            zero frequency. 

        Examples
        --------
        
        >>> import quantarhei as qr
        >>> lab = LabSetup()
        >>> 
        >>> # Time axis around 0
        >>> time = qr.TimeAxis(-500.0, 1000, 1.0)
        >>> # FrequencyAxis around 0 
        >>> freq = qr.FrequencyAxis(-2500, 1000, 5.0)
        
        `numeric` pulse shape if time domain
        
        >>> pulse1 = dict(ptype="numeric")
        >>> params = (pulse1, pulse1, pulse1)
        >>> lab.set_pulse_shapes(time, params)
        
        
        Gaussian pulse shape in time domain
        
        >>> pulse2 = dict(ptype="Gaussian", FWHM=10, amplitude=1.0)
        >>> params = (pulse2, pulse2, pulse2)
        >>> lab.set_pulse_shapes(time, params)
        
        Situations in which Exceptions are thrown
        
        >>> pulse3 = dict(ptype="other", FWHM=10, amplitude=1.0)
        >>> params = (pulse3, pulse3, pulse3)
        >>> lab.set_pulse_shapes(time, params)
        Traceback (most recent call last):
            ...
        Exception: Unknown pulse type
        
        >>> params = (pulse2, pulse2)
        >>> lab.set_pulse_shapes(time, params)
        Traceback (most recent call last):
            ...
        Exception: set_pulses requires 3 parameter sets
        
        
        """
        
        if isinstance(axis, TimeAxis):
            self.timeaxis = axis
            self.axis_type = "time"
            
        elif isinstance(axis, FrequencyAxis):
            self.freqaxis = axis
            self.axis_type = "frequency"
            
        else:
            raise Exception("Wrong axis paramater")

        if len(params) == self.number_of_pulses:
        
            k_p = 0
            for par in params:
            
                if par["ptype"] == "Gaussian":
                    
                    if self.axis_type == "time":
                        pass
                    
                    elif self.axis_type == "frequency":
                        pass
        
                elif par["ptype"] == "numeric":
            
                    fce = par["function"]

                    if self.axis_type == "time":
                        
                        #
                        # Create a new DFunction based on the submitted time 
                        # axis
                        #
                        data = numpy.zeros(self.timeaxis.length)
                        for t_p in self.timeaxis.data:
                            data[t_p] = fce.at(t_p)
                            
                        self.pulse_t[k_p] = DFunction(self.timeaxis, data)
                    
                    elif self.axis_type == "frequency":
                        
                        data = numpy.zeros(self.freqaxis.length)
                        for t_p in self.freqaxis.data:
                            data[t_p] = fce.at(t_p)
                            
                        self.pulse_f[k_p] = DFunction(self.freqaxis, data)
                    
        
                else:
                    raise Exception("Unknown pulse type")
                    
                k_p += 1
                
        else:
            text = "set_pulses requires "+str(self.number_of_pulses) \
                    +" parameter sets"
            raise Exception(text)
            
        

    
    def set_polarizations(self, pulse_polarizations=(X, X, X), 
                         detection_polarization=X):
        """Sets polarizations of the experimental pulses
        
        
        Parameters
        ----------
        
        pulse_polarization : tuple like
            Contains three vectors of polarization of the three pulses
            of the experiment. Currently we assume three pulse experiment
            per default.
            
        detection_polarization : array
            Vector of detection polarization
            
            
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> lab = LabSetup()
        >>> lab.set_polarizations(pulse_polarizations=(qr.utils.vectors.X,
        ...                                            qr.utils.vectors.Y,
        ...                                            qr.utils.vectors.Z))
        >>> print(lab.e[0,:])
        [ 1.  0.  0.]
        >>> print(lab.e[3,:])
        [ 1.  0.  0.]
        >>> print(lab.e[2,:])
        [ 0.  0.  1.]
        
        
        >>> lab.set_polarizations(pulse_polarizations=(qr.utils.vectors.X,
        ...                                            qr.utils.vectors.Y)) 
        Traceback (most recent call last):
            ...
        Exception: pulse_polarizations requires 3 values
        
        """
        if len(pulse_polarizations) == self.number_of_pulses:

            self.e = numpy.zeros((4,3))
            for i in range(3):
                self.e[i,:] = pulse_polarizations[i]
            self.e[3,:] = detection_polarization
            
            e = self.e
            
            F4e = numpy.zeros(3)
            F4e[0] = numpy.dot(e[3,:],e[2,:])*numpy.dot(e[1,:],e[0,:])
            F4e[1] = numpy.dot(e[3,:],e[1,:])*numpy.dot(e[2,:],e[0,:])
            F4e[2] = numpy.dot(e[3,:],e[0,:])*numpy.dot(e[2,:],e[1,:])
            
            self.F4eM4 = numpy.dot(F4e,self.M4)
            
            
        else:
            text = "pulse_polarizations requires "+ \
                    str(self.number_of_pulses)+" values"
            raise Exception(text)
            
        self.detection_polarization = detection_polarization

     
    def get_pulse_polarizations(self):
        """Returns polarizations of the laser pulses
        
        
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> lab = LabSetup()
        >>> lab.set_polarizations(pulse_polarizations=(qr.utils.vectors.X,
        ...                                            qr.utils.vectors.Y,
        ...                                            qr.utils.vectors.Z))
        >>> pols = lab.get_pulse_polarizations()
        >>> print(len(pols))
        3
        
        """
        pols = []
        for i in range(self.number_of_pulses):
            pols.append(self.e[i,:])
            
        return pols

        
    def get_detection_polarization(self):
        """Returns detection polarizations
        
        
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> lab = LabSetup()
        >>> lab.set_polarizations(pulse_polarizations=(qr.utils.vectors.X,
        ...                                            qr.utils.vectors.Y,
        ...                                            qr.utils.vectors.Z))
        >>> detpol = lab.get_detection_polarization()
        >>> print(detpol)
        [ 1.  0.  0.]
        
        """
        return self.e[3,:]


    def convert_to_time(self):
        pass

    
    def convert_to_frequency(self):
        pass

    
    def get_pulse_envelop(self, k, t):
        """Returns a numpy array with the pulse time-domain envelope
        
        """
        return self.pulse_t[k].at(t)
    
    
    def get_pulse_spectrum(self, k, omega):
        """Returns a numpy array with the pulse frequency-domain spectrum
        
        """
        return self.pulse_f[k].at(omega)
    
 
class labsetup(LabSetup):
    pass
