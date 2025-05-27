# -*- coding: utf-8 -*-
"""
    Laboratory set-up for non-linear spectroscopy
    
    This class controls calculations of non-linear optical spectra, and
    other experiments in which laboratory setting needs to be controlled.
    Examples are pulse polarization setting, pulse shapes and spectra 
    in non-linear spectroscopy.
            

    Class Details
    -------------

"""

import numpy
from functools import partial

from ..utils import Integer
from ..utils.vectors import X
from ..core.time import TimeAxis
from ..core.frequency import FrequencyAxis
from ..core.dfunction import DFunction
from .. import REAL, COMPLEX
from .. import Manager

class LabSetup:
    """Laboratory set-up for non-linear spectroscopy
    
    Class representing laboratory setup for non-linear spectroscopic
    experiments. It holds information about pulse shapes and polarizations.
    
    Pulses can be set in time- and/or frequency-domain. **Consistency between
    the domains is not checked nor enforced**. Consistent conversion between
    domains is provided by convenience routines [TO BE IMPLEMENTED]


    Parameters
    ----------

    nopulses : int
        Number of pulses in the experiment. Default is 3.



    """
    
    number_of_pulses = Integer("number_of_pulses")
    
    def __init__(self, nopulses = 3):
        
        # number of pulses in the set-up
        self.number_of_pulses = nopulses
        self.pulse_effects = "rescale_dip"  # Pulse shape effects accounted
                                        # for by rescaling transition dipoles
                                        # When the pulses are not defined no 
                                        # rescaling is used.
    
        # orientational averaging matrix for four-wave-mixing
        self.M4 = numpy.array([[4.0, -1.0, -1.0],
                               [-1.0, 4.0, -1.0],
                               [-1.0,-1.0,  4.0]])/30.0
    
        # auxiliary matrix for orientational averaging
        self.F4eM4 = None
        
        # pulse polarizations
        self.e = None
    
        self.timeaxis = None
        self.freqaxis = None
                        

        self.has_polarizations = False
        self.has_freqdomain = False
        self.has_timedomain = False

        # time or frequency
        self.axis_type = None

        # pulses in time- and frequency domain
        self.pulse_t = [None]*nopulses
        self.pulse_f = [None]*nopulses
        
        self._field_set = False

        self.dscaling = None
        
        #
        # Pulse characteristics
        #
        self.omega = numpy.zeros(nopulses, dtype=REAL)
        self.saved_omega = None
        self.pulse_centers = numpy.zeros(nopulses, dtype=REAL)
        self._centers_set = False
        self.phases = numpy.zeros(nopulses, dtype=REAL)
        self.delay_phases = numpy.zeros(nopulses, dtype=REAL)
        self.e = numpy.zeros((nopulses,3), dtype=REAL)
        
        self.saved_params = None
        
    
    def reset_pulse_shape(self):
        """Recalculates the pulse shapes 
        
        """

        if self.saved_params is not None:
            self.set_pulse_shapes(self.timeaxis, self.saved_params)
        #else:
        #    raise Exception("Pulse shapes must be set first.")
            
        

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
            
        params : dictionary
            Dictionary of pulse parameters. The parameters are the following:
            `ptype` is the pulse type with possible values `Gaussian` and 
            `numeric`. Time domain pulses are specified with their center
            at t = 0.
            
            **Gaussian** pulse has further parameters `amplitude`, `FWHM`,
            and `frequency` with obvious meanings. `FWHM` is speficied in `fs`,
            `frequency` is specified in energy units, while `amplitude`
            is in units of [energy]/[transition dipole moment]. The formula
            for the lineshape is
            
            .. math::
                 
                \\rm{shape}(\\omega) = 
                \\frac{2}{\\Delta}\\sqrt{\\frac{\\ln(2)}{\\pi}}
                \\exp\\left\\{-\\frac{4\\ln(2)\\omega^2}{\\Delta^2}\\right\\}

            The same formulae are used for time- and frequency domain
            definitions. For time domain, :math:`t` should be used in stead of 
            :math:`\omega`.

            **numeric** pulse is specified by a second parameters `function` 
            which should be of DFunction type and specifies line shape around
            zero frequency. 


        Examples
        --------
        
        >>> import quantarhei as qr
        >>> import matplotlib.pyplot as plt
        >>> lab = LabSetup()
        ... 
        >>> # Time axis around 0
        >>> time = qr.TimeAxis(-500.0, 1000, 1.0, atype="complete")

        Gaussian pulse shape in time domain
        
        >>> pulse2 = dict(ptype="Gaussian", FWHM=150, amplitude=1.0)
        >>> params = (pulse2, pulse2, pulse2)
        >>> lab.set_pulse_arrival_times([0.0, 0.0, 0.0])
        >>> lab.set_pulse_phases([0.0, 0.0, 0.0])  # these settings are compulsory
        >>> lab.set_pulse_shapes(time, params)
        
        Testing the pulse shape
                
        >>> dfc = lab.get_pulse_envelop(1, time.data) # doctest: +SKIP
        >>> pl = plt.plot(time.data, dfc)             # doctest: +SKIP
        >>> plt.show()                                # doctest: +SKIP
        

        .. plot::
            
            import quantarhei as qr
            import matplotlib.pyplot as plt
            
            lab = qr.LabSetup()
            time = qr.TimeAxis(-500.0, 1000, 1.0, atype="complete")
            pulse2 = dict(ptype="Gaussian", FWHM=150.0, amplitude=1.0)
            params = (pulse2, pulse2, pulse2)
            lab.set_pulse_arrival_times([0.0, 0.0, 0.0])
            lab.set_pulse_phases([0.0, 0.0, 0.0])
            lab.set_pulse_shapes(time, params)
            
            dfc = lab.get_pulse_envelop(1, time.data)
            pl = plt.plot(time.data, dfc)
            plt.show()

        
        `numeric` pulse shape in time domain
        
        >>> # We take the DFunction for creation of `numeric`ly defined
        >>> # pulse shape from the previous example
        >>> pls = lab.pulse_t[2]
        >>> # new lab object
        >>> lab2 = LabSetup()
        
        >>> pulse1 = dict(ptype="numeric", function=pls)
        >>> params = (pulse1, pulse1, pulse1)
        >>> lab2.set_pulse_arrival_times([0.0, 0.0, 0.0])
        >>> lab2.set_pulse_phases([0.0, 0.0, 0.0])
        >>> lab2.set_pulse_shapes(time, params)
        
        Testing the pulse shape
        
        >>> dfc = lab2.get_pulse_envelop(1, time.data) # doctest: +SKIP
        >>> pl = plt.plot(time.data, dfc)              # doctest: +SKIP
        >>> plt.show() # we skip output here           # doctest: +SKIP
        
        
        Gaussian pulse shape in frequency domain
        
        >>> lab = LabSetup()
        >>> # FrequencyAxis around 0 
        >>> freq = qr.FrequencyAxis(-2500, 1000, 5.0)
        ... 
        >>> pulse2 = dict(ptype="Gaussian", FWHM=800, amplitude=1.0)
        >>> params = (pulse2, pulse2, pulse2)
        >>> lab.set_pulse_arrival_times([0.0, 0.0, 0.0])
        >>> lab.set_pulse_phases([0.0, 0.0, 0.0])
        >>> lab.set_pulse_shapes(freq, params)
        
        Testing the pulse shape
        
        >>> # getting differnt frequency axis
        >>> freq2 = qr.FrequencyAxis(-1003, 100, 20.0)
        >>> # and reading spectrum at two different sets of points
        >>> dfc1 = lab.get_pulse_spectrum(1, freq.data)  
        >>> dfc2 = lab.get_pulse_spectrum(1, freq2.data) 
        >>> pl1 = plt.plot(freq.data, dfc1)             # doctest: +SKIP
        >>> pl2 = plt.plot(freq2.data, fdc2)            # doctest: +SKIP
        >>> plt.show()                                  # doctest: +SKIP
        
        We plot in two different sets of points.
        
        .. plot::
            
            import quantarhei as qr
            import matplotlib.pyplot as plt
            
            lab = qr.LabSetup()
            freq = qr.FrequencyAxis(-2500, 1000, 5.0)
            pulse2 = dict(ptype="Gaussian", FWHM=800.0, amplitude=1.0)
            params = (pulse2, pulse2, pulse2)
            lab.set_pulse_arrival_times([0.0, 0.0, 0.0])
            ab.set_pulse_phases([0.0, 0.0, 0.0])
            lab.set_pulse_shapes(freq, params)
            
            freq2 = qr.FrequencyAxis(-1000, 100, 20.0)
           
            dfc1 = lab.get_pulse_spectrum(1, freq.data)
            dfc2 = lab.get_pulse_spectrum(1, freq2.data)
            pl1 = plt.plot(freq.data, dfc1) 
            pl2 = plt.plot(freq2.data, dfc2)
            plt.show()


        `numeric` pulse shape in frequency domain
        
        >>> # We take the DFunction for creation of `numeric`ly defined
        >>> # pulse shape from the previous example
        >>> pls = lab.pulse_f[2]
        >>> # new lab object
        >>> lab2 = LabSetup()
        
        >>> pulse1 = dict(ptype="numeric", function=pls)
        >>> params = (pulse1, pulse1, pulse1)
        >>> lab2.set_pulse_arrival_times([0.0, 0.0, 0.0])
        >>> lab2.set_pulse_phases([0.0, 0.0, 0.0])
        >>> lab2.set_pulse_shapes(freq, params)
        
        Testing the pulse shape
        
        >>> dfc = lab2.get_pulse_envelop(1, freq.data) # doctest: +SKIP
        >>> pl = plt.plot(freq.data, dfc)              # doctest: +SKIP
        >>> plt.show() # we skip output here           # doctest: +SKIP

   
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
        

        >>> params = (pulse2, pulse2)
        >>> lab.set_pulse_shapes(time.data, params)
        Traceback (most recent call last):
            ...
        Exception: Wrong axis paramater
        
        >>> time = qr.TimeAxis(0.0, 1000, 1.0)
        >>> lab.set_pulse_shapes(time, params)
        Traceback (most recent call last):
            ...
        Exception: TimeAxis has to be of 'complete' type use atype='complete' as a parameter of TimeAxis

        
        """
        
        if isinstance(axis, TimeAxis):
            if axis.atype == "complete":
                self.timeaxis = axis
                self.axis_type = "time"
            else:
                raise Exception("TimeAxis has to be of 'complete' type"+
                                " use atype='complete' as a parameter"+
                                " of TimeAxis")
            
        elif isinstance(axis, FrequencyAxis):
            self.freqaxis = axis
            self.axis_type = "frequency"
            
        else:
            raise Exception("Wrong axis paramater")


        if not self._centers_set:
            #print("Use 'set_pulse_arrival_times' function before this function.")
            #raise Exception("Pulse arrival times have to specified before "+
            #                "the pulse shape is set.")
            self.set_pulse_arrival_times([0.0 for ii in range(self.number_of_pulses)])


        if len(params) == self.number_of_pulses:
        
            self.saved_params = params    
        
            k_p = 0
            for par in params:
            
                if par["ptype"] == "Gaussian":
                    
                    if self.axis_type == "time":
                        #
                        # Time domain Gaussian pulse around 0.0 as a DFunction
                        #
                        tma = self.timeaxis
                        fwhm = par["FWHM"]
                        amp = par["amplitude"]
                        
                        tc = self.pulse_centers[k_p]
                        
                        # normalized Gaussian mupliplied by amplitude
                        lfc = 4.0*numpy.log(2.0)
                        pi = numpy.pi
                        val = (2.0/fwhm)*numpy.sqrt(numpy.log(2.0)/pi) \
                        *amp*numpy.exp(-lfc*((tma.data-tc)/fwhm)**2)
                        
                        self.pulse_t[k_p] = DFunction(tma, val)
                    
                    elif self.axis_type == "frequency":
                        #
                        # Frequency domain Gaussian pulse around 0.0 
                        # as a DFunction
                        #
                        fra = self.freqaxis
                        fwhm = par["FWHM"]
                        amp = par["amplitude"]
                        try:
                            freq = par["frequency"]
                        except:
                            freq = 0.0
                        
                        # normalized Gaussian mupliplied by amplitude
                        val = (2.0/fwhm)*numpy.sqrt(numpy.log(2.0)/3.14159) \
                        *amp*numpy.exp(-4.0*numpy.log(2.0)\
                                       *((fra.data-freq)/fwhm)**2)
                        
                        self.pulse_f[k_p] = DFunction(fra, val) 
                        
                elif par["ptype"] == "numeric":
            
                    fce = par["function"]

                    if self.axis_type == "time":
                        
                        #
                        # Create a new DFunction based on the submitted time 
                        # axis
                        #
                        data = numpy.zeros(self.timeaxis.length)
                        i_p = 0
                        for t_p in self.timeaxis.data:
                            data[i_p] = fce.at(t_p)
                            i_p += 1
                            
                        self.pulse_t[k_p] = DFunction(self.timeaxis, data)
                    
                    elif self.axis_type == "frequency":
                        
                        data = numpy.zeros(self.freqaxis.length)
                        i_p = 0
                        for t_p in self.freqaxis.data:
                            data[i_p] = fce.at(t_p)
                            i_p += 1
                            
                        self.pulse_f[k_p] = DFunction(self.freqaxis, data)
                    
        
                else:
                    raise Exception("Unknown pulse type")
                    
                k_p += 1
             
            if self.axis_type == "time":
                self.has_timedomain = True
            elif self.axis_type == "frequency":
                self.has_freqdomain = True
                
            self._field_set = True
                
        else:
            text = "set_pulses requires "+str(self.number_of_pulses) \
                    +" parameter sets"
            raise Exception(text)
    
    

    def set_pulse_polarizations(self, pulse_polarizations=(X, X, X), 
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
        >>> lab.set_pulse_polarizations(pulse_polarizations=(qr.utils.vectors.X,
        ...                                            qr.utils.vectors.Y,
        ...                                            qr.utils.vectors.Z))
        >>> print(lab.e[0,:])
        [ 1.  0.  0.]
        >>> print(lab.e[3,:])
        [ 1.  0.  0.]
        >>> print(lab.e[2,:])
        [ 0.  0.  1.]
        
        
        >>> lab.set_pulse_polarizations(pulse_polarizations=(qr.utils.vectors.X,
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
        >>> lab.set_pulse_polarizations(pulse_polarizations=(qr.utils.vectors.X,
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
        >>> lab.set_pulse_polarizations(pulse_polarizations=(qr.utils.vectors.X,
        ...                                            qr.utils.vectors.Y,
        ...                                            qr.utils.vectors.Z))
        >>> detpol = lab.get_detection_polarization()
        >>> print(detpol)
        [ 1.  0.  0.]
        
        """
        return self.e[3,:]


    def convert_to_time(self):
        """Converts pulse information from frequency domain to time domain
        
        
        
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> import matplotlib.pyplot as plt
        >>> lab = LabSetup()
        >>> freq = qr.FrequencyAxis(-100, 200, 1.0) # atype="complete" is default
        >>> pulse = dict(ptype="Gaussian", FWHM=20, amplitude=1.0)
        >>> params = (pulse, pulse, pulse)
        >>> lab.set_pulse_arrival_times([0.0, 0.0, 0.0])
        >>> lab.set_pulse_phases([0.0, 0.0, 0.0])
        >>> lab.set_pulse_shapes(freq, params)
        >>> lab.convert_to_time()
        >>> # plot the original and the FT pulses
        >>> pls_1f = lab.pulse_f[1]                          # doctest: +SKIP
        >>> p1 = plt.plot(pls_1f.axis.data, pls_1f.data)     # doctest: +SKIP
        >>> pls_1t = lab.pulse_t[1]                          # doctest: +SKIP
        >>> p2 = plt.plot(pls_1t.axis.data, pls_1t.data)     # doctest: +SKIP 
        >>> plt.show()                                       # doctest: +SKIP
        
        .. plot::

            import quantarhei as qr
            import matplotlib.pyplot as plt
            lab = qr.LabSetup()
            freq = qr.FrequencyAxis(-100,200,1.0)
            pulse = dict(ptype="Gaussian", FWHM=5, amplitude=1.0)
            params = (pulse, pulse, pulse)
            lab.set_pulse_arrival_times([0.0, 0.0, 0.0])
            lab.set_pulse_phases([0.0, 0.0, 0.0])
            lab.set_pulse_shapes(freq, params)
            lab.convert_to_time()
            
            pls_1f = lab.pulse_f[1]
            plt.plot(pls_1f.axis.data, pls_1f.data)
            pls_1t = lab.pulse_t[1]
            plt.plot(pls_1t.axis.data, pls_1t.data)
            plt.show()


        Now we compare back and forth Fourier transform with the original

        >>> import quantarhei as qr
        >>> import numpy
        >>> lab = LabSetup()
        >>> freq = qr.FrequencyAxis(-100,200,1.0) # atype="complete" is default
        >>> pulse = dict(ptype="Gaussian", FWHM=20, amplitude=1.0)
        >>> params = (pulse, pulse, pulse)
        >>> lab.set_pulse_arrival_times([0.0, 0.0, 0.0])
        >>> lab.set_pulse_phases([0.0, 0.0, 0.0])        
        >>> lab.set_pulse_shapes(freq, params)
        >>> freq_vals_1 = lab.get_pulse_spectrum(2, freq.data)
        >>> lab.convert_to_time()
        
        Here we override the original frequency domain definition
        
        >>> lab.convert_to_frequency()
        >>> freq_vals_2 = lab.get_pulse_spectrum(2, freq.data)
        >>> numpy.allclose(freq_vals_2, freq_vals_1)
        True
        
        and now the other way round
        
        >>> import quantarhei as qr
        >>> import numpy
        >>> lab = LabSetup()
        >>> time = qr.TimeAxis(-100,200,1.0, atype="complete")
        >>> pulse = dict(ptype="Gaussian", FWHM=20, amplitude=1.0)
        >>> params = (pulse, pulse, pulse)
        >>> lab.set_pulse_arrival_times([0.0, 0.0, 0.0])
        >>> lab.set_pulse_phases([0.0, 0.0, 0.0])        
        >>> lab.set_pulse_shapes(time, params)  
        >>> time_vals_1 = lab.get_pulse_envelop(2, time.data)
        >>> lab.convert_to_frequency()        
 
        Here we override the original time domain definition
        
        >>> lab.convert_to_time()
        >>> time_vals_2 = lab.get_pulse_envelop(2, freq.data)
        >>> numpy.allclose(time_vals_2, time_vals_1)
        True
        
        
        Situation in which excetions are thrown
        
        >>> lab = LabSetup()
        >>> lab.convert_to_time()
        Traceback (most recent call last):
            ...
        Exception: Cannot convert to time domain: frequency domain not set
        
        
        """
        if self.has_freqdomain:
            
            freq = self.freqaxis
            time = freq.get_TimeAxis()
            
            k_p = 0
            for pulse in self.pulse_f:
                ft_pulse = pulse.get_Fourier_transform()
                # we replace the DFunction's axis attribute with the one
                # calculated above; in time domain the pulses also share
                # the same TimeAxis object
                ft_pulse.axis = time
                self.pulse_t[k_p] = ft_pulse
                k_p += 1
            
            self.timeaxis = time
            self.has_timedomain = True
        
        else:
            raise Exception("Cannot convert to time domain: "+ 
                            "frequency domain not set")

    
    def convert_to_frequency(self):
        """Converts pulse information from time domain to frequency domain
        
        
        
        Examples
        --------

        >>> import quantarhei as qr
        >>> lab = LabSetup()
        >>> time = qr.TimeAxis(-100,200,1.0, atype="complete")
        >>> pulse = dict(ptype="Gaussian", FWHM=20, amplitude=1.0)
        >>> params = (pulse, pulse, pulse)
        >>> lab.set_pulse_arrival_times([0.0, 0.0, 0.0])
        >>> lab.set_pulse_phases([0.0, 0.0, 0.0])
        >>> lab.set_pulse_shapes(time, params)
        >>> lab.convert_to_frequency()
        >>> # plot the original and the FT pulses
        >>> pls_1f = lab.pulse_f[1]                      # doctest: +SKIP
        >>> plt.plot(pls_1f.axis.data, pls_1f.data)      # doctest: +SKIP
        >>> pls_1t = lab.pulse_t[1]                      # doctest: +SKIP
        >>> plt.plot(pls_1t.axis.data, pls_1t.data)      # doctest: +SKIP
        >>> plt.show()                                   # doctest: +SKIP
        
        .. plot::

            import quantarhei as qr
            import matplotlib.pyplot as plt
            lab = qr.LabSetup()
            time = qr.TimeAxis(-100,200,1.0, atype="complete")
            pulse = dict(ptype="Gaussian", FWHM=5, amplitude=1.0)
            params = (pulse, pulse, pulse)
            lab.set_pulse_arrival_times([0.0, 0.0, 0.0])
            lab.set_pulse_phases([0.0, 0.0, 0.0])
            lab.set_pulse_shapes(time, params)
            lab.convert_to_frequency()
            
            pls_1f = lab.pulse_f[1]
            plt.plot(pls_1f.axis.data, pls_1f.data)
            pls_1t = lab.pulse_t[1]
            plt.plot(pls_1t.axis.data, pls_1t.data)
            plt.show()
            
        
        Situation in which excetions are thrown
        
        >>> lab = LabSetup()
        >>> lab.convert_to_frequency()
        Traceback (most recent call last):
            ...
        Exception: Cannot convert to frequency domain: time domain not set
        
        
        """
        if self.has_timedomain:
            
            time = self.timeaxis
            freq = time.get_FrequencyAxis()
            
            k_p = 0
            for pulse in self.pulse_t:
                ft_pulse = pulse.get_Fourier_transform()
                # we replace the DFunction's axis attribute with the one
                # calculated above; in time domain the pulses also share
                # the same TimeAxis object
                ft_pulse.axis = freq
                self.pulse_f[k_p] = ft_pulse
                k_p += 1
            
            self.freqaxis = freq
            self.has_freqdomain = True
                
        else:
            raise Exception("Cannot convert to frequency domain: "+ 
                            "time domain not set")
        
    
    def get_pulse_envelop(self, k, t):
        """Returns a numpy array with the pulse time-domain envelope
        
        
        Parameters
        ----------
        
        k : int
            Index of the pulse to be returned
            
        t : array like
            Array of time points at which the pulse is returned
            
            
        Examples
        --------

        >>> import quantarhei as qr
        >>> lab = LabSetup()
        >>> time = qr.TimeAxis(-100, 200, 1.0, atype="complete")
        >>> pulse2 = dict(ptype="Gaussian", FWHM=30.0, amplitude=1.0)
        >>> params = (pulse2, pulse2, pulse2)
        >>> lab.set_pulse_arrival_times([0.0, 0.0, 0.0])
        >>> lab.set_pulse_phases([0.0, 0.0, 0.0])
        >>> lab.set_pulse_shapes(time, params)
        >>> dfc = lab.get_pulse_envelop(1, [-50.0, -30.0, 2.0, 30.0])
        >>> print(dfc)
        [  1.41569209e-05   1.95716100e-03   3.09310662e-02   1.95716100e-03]
        
        .. plot::
            :include-source:
            
            import quantarhei as qr
            import matplotlib.pyplot as plt
            
            lab = qr.LabSetup()
            time = qr.TimeAxis(-500.0, 1000, 1.0, atype="complete")
            pulse2 = dict(ptype="Gaussian", FWHM=150.0, amplitude=1.0)
            params = (pulse2, pulse2, pulse2)
            lab.set_pulse_arrival_times([0.0, 0.0, 0.0])
            lab.set_pulse_phases([0.0, 0.0, 0.0])
            lab.set_pulse_shapes(time, params)
            
            pls = lab.pulse_t[2]
            lab2 = qr.LabSetup()
            
            pulse1 = dict(ptype="numeric", function=pls)
            params = (pulse1, pulse1, pulse1)
            lab2.set_pulse_arrival_times([0.0, 0.0, 0.0])
            lab2.set_pulse_phases([0.0, 0.0, 0.0])
            lab2.set_pulse_shapes(time, params)
            
            dfc = lab2.get_pulse_envelop(1, time.data)
            pl = plt.plot(time.data, dfc)
            plt.show()

        
        """
        return self.pulse_t[k].at(t)
    
    
    def get_pulse_spectrum(self, k, omega):
        """Returns a numpy array with the pulse frequency-domain spectrum

        Parameters
        ----------
        
        k : int
            Index of the pulse to be returned
            
        omega : array like
            Array of frequency points at which the pulse is returned
            
            
        Examples
        --------

        >>> import quantarhei as qr
        >>> lab = LabSetup()
        >>> freq = qr.FrequencyAxis(-2500, 1000, 5.0)
        >>> pulse2 = dict(ptype="Gaussian", FWHM=800.0, amplitude=1.0)
        >>> params = (pulse2, pulse2, pulse2)
        >>> lab.set_pulse_arrival_times([0.0, 0.0, 0.0])
        >>> lab.set_pulse_phases([0.0, 0.0, 0.0])
        >>> lab.set_pulse_shapes(freq, params)
        >>> dfc = lab.get_pulse_spectrum(1, [600.0, 700.0, 800.0, 900.0])
        >>> print(dfc)
        [  2.46865554e-04   1.40563844e-04   7.33935684e-05   3.51409609e-05]
        
        Here is a complete example with setting, getting and plotting spectrum:

        .. plot::
            :include-source:
                
            import quantarhei as qr
            import matplotlib.pyplot as plt
            
            lab = qr.LabSetup()
            freq = qr.FrequencyAxis(-2500, 1000, 5.0)
            pulse2 = dict(ptype="Gaussian", FWHM=800.0, amplitude=1.0)
            params = (pulse2, pulse2, pulse2)
            lab.set_pulse_arrival_times([0.0, 0.0, 0.0])
            lab.set_pulse_phases([0.0, 0.0, 0.0])
            lab.set_pulse_shapes(freq, params)
            
            pls = lab.pulse_f[2]
            lab2 = qr.LabSetup()
            
            pulse1 = dict(ptype="numeric", function=pls)
            params = (pulse1, pulse1, pulse1)
            lab2.set_pulse_arrival_times([0.0, 0.0, 0.0])
            lab2.set_pulse_phases([0.0, 0.0, 0.0])
            lab2.set_pulse_shapes(freq, params)
            
            dfc = lab2.get_pulse_spectrum(1, freq.data)
            pl = plt.plot(freq.data, dfc)
            plt.show()

        
        """
        return self.pulse_f[k].at(omega)
    
    
    def set_pulse_frequencies(self, omegas):
        """Sets pulse frequencies
        
        
        Parameters
        ----------
        
        omegas : array of floats
            Frequencies of pulses
        
        Examples
        --------
        
        >>> lab = LabSetup()
        >>> lab.set_pulse_frequencies([1.0, 2.0, 1.0])
        >>> print(lab.omega)
        [ 1.  2.  1.]
        
        Situation which throws an exception

        >>> lab = LabSetup()
        >>> lab.set_pulse_frequencies([1.0, 2.0, 1.0, 6.0])       
        Traceback (most recent call last):
            ...
        Exception: Wrong number of frequencies: 3 required
        
        """
        
        # FIXME: energy unit control has to be in place
        if len(omegas) == self.number_of_pulses:

            self.omega = Manager().convert_energy_2_internal_u(numpy.array(omegas, dtype=REAL))
            
        else:
            raise Exception("Wrong number of frequencies: "+
                            str(self.number_of_pulses)+" required")


    def get_pulse_frequency(self, k):
        """Returns frequency of the pulse with index k
        
        Parameters
        ----------
        
        k : int
            Pulse index
           
        Examples
        --------
        
        >>> lab = LabSetup()
        >>> lab.set_pulse_frequencies([1.0, 2.0, 1.0])
        >>> print(lab.get_pulse_frequency(1))
        2.0
            
        """
        return self.omega[k]
    
    
    def set_pulse_arrival_times(self, times):
        """Sets the arrival time (i.e. centers) of the pulses
 
        Parameters
        ----------
        
        times : array of floats
            Arrival times (centers) of the pulses
        
        Examples
        --------
        
        >>> lab = LabSetup()
        >>> lab.set_pulse_arrival_times([1.0, 20.0, 100.0])
        >>> print(lab.pulse_centers)
        [1.0, 20.0, 100.0]
        
        Situation which throws an exception

        >>> lab = LabSetup()
        >>> lab.set_pulse_arrival_times([1.0, 2.0, 1.0, 6.0])       
        Traceback (most recent call last):
            ...
        Exception: Wrong number of arrival times: 3 required
 
        """
        if len(times) == self.number_of_pulses:
            
            self.pulse_centers = times
            self._centers_set = True
            self.reset_pulse_shape()
            
        else:
            raise Exception("Wrong number of arrival times: "+
                            str(self.number_of_pulses)+" required")        
        

    def get_pulse_arrival_times(self):
        """Returns frequency of the pulse with index k
        
           
        Examples
        --------
        
        >>> lab = LabSetup()
        >>> lab.set_pulse_arrival_times([1.0, 20.0, 100.0])
        >>> print(lab.get_pulse_arrival_times())
        [1.0, 20.0, 100.0]
            
        """
        return self.pulse_centers


    def get_pulse_arrival_time(self, k):
        """Returns frequency of the pulse with index k
        
        Parameters
        ----------
        
        k : int
            Pulse index
           
        Examples
        --------
        
        >>> lab = LabSetup()
        >>> lab.set_pulse_arrival_times([1.0, 20.0, 100.0])
        >>> print(lab.get_pulse_arrival_time(1))
        20.0
            
        """
        return self.pulse_centers[k]


    def set_pulse_phases(self, phases):
        """Sets the phases of the individual pulses 
 
        Parameters
        ----------
        
        phases : array of floats
            Phases of the pulses
        
        Examples
        --------
        
        >>> lab = LabSetup()
        >>> lab.set_pulse_phases([1.0, 3.14, -1.0])
        >>> print(lab.phases)
        [1.0, 3.14, -1.0]
        
        Situation which throws an exception

        >>> lab = LabSetup()
        >>> lab.set_pulse_phases([1.0, 2.0, 1.0, 6.0])       
        Traceback (most recent call last):
            ...
        Exception: Wrong number of phases: 3 required
 
        """
        if len(phases) == self.number_of_pulses:
            
            self.phases = phases
            
        else:
            raise Exception("Wrong number of phases: "+
                            str(self.number_of_pulses)+" required") 


    def get_pulse_phases(self):
        """Returns frequency of the pulse with index k
        
           
        Examples
        --------
        
        >>> lab = LabSetup()
        >>> lab.set_pulse_phases([1.0, 3.14, -1.0])
        >>> print(lab.get_pulse_phases())
        [1.0, 3.14, -1.0]
            
        """
        return self.phases


    def get_pulse_phase(self, k):
        """Returns frequency of the pulse with index k
        
        Parameters
        ----------
        
        k : int
            Pulse index
           
        Examples
        --------
        
        >>> lab = LabSetup()
        >>> lab.set_pulse_phases([1.0, 3.14, -1.0])
        >>> print(lab.get_pulse_phase(1))
        3.14
            
        """
        return self.phases[k]
    
    
    # def get_field(self, k, rwa=0.0):
    #     """Returns an EField object corresponding to the k-th field
        
    #     """

    #     ef = EField(self,omega=self.omega[k],polar=self.e[k,:],
    #                 ftype="Data", data=self.pulse_t[k].data)
    #     ef.subtract_frequency(rwa)
        
    #     return ef
    
    
    # def get_fields(self, rwa=0.0):
    #     """Retruns a list of EField objects for the lab's pulses
        
    #     """
        
    #     fields = [None]*self.number_of_pulses
    #     for kk in range(self.number_of_pulses):
            
    #         ff = self.get_field(kk, rwa=rwa)
    #         fields[kk] = ff
            
    #     return fields
    
    def get_labfield(self, k):
        return LabField(self, k)

        

    def get_labfields(self):
        """Retruns a list of EField objects for the lab's pulses
        
        """
        
        fields = [None]*self.number_of_pulses
        
        for kk in range(self.number_of_pulses):
            
            ff = self.get_labfield(kk)
            fields[kk] = ff
            
        return fields
    
    
    def set_rwa(self, om):
        
        self.saved_omega = numpy.zeros((self.number_of_pulses), dtype=REAL)
        
        self.saved_omega[:] = self.omega[:]
        self.omega[:] -= om
        
    def restore_rwa(self):
        
        if self.saved_omega is None:
            raise Exception("RWA has to be set first")
            
        self.omega[:] = self.saved_omega[:] 
        
        
    def get_field(self, kk=None):
        """Returns the total field of the lab or a single field
        
        """
        #N = t.shape[0]
        
        #fld = numpy.zeros(N, dtype=COMPLEX)
        if kk is None:

            flds = self.get_labfields()
            
            kk = 0
            for fl in flds:
                if kk == 0:
                    fld = fl.get_field()
                else:
                    fld += fl.get_field()
                kk += 1
                
            return fld
        
        else:

            fl = self.get_labfield(kk)
            return fl.get_field()


    def get_field_derivative(self):
        """ Returns the time derivative of the total field

        """

        flds = self.get_labfields()
        
        kk = 0
        for fl in flds:

            if kk == 0:
                fld_d = fl.get_field_derivative()
            else:
                fld_d += fl.get_field_derivative()
            
            kk += 1
            
        return fld_d        
        

 
class labsetup(LabSetup):
    """labsetup is just a different name for the class LabSetup.
    
    All details about usage of the labsetup can be found in the documentation
    of the LabSetup
    
    """
    pass



def _labattr(name, target, flag=None):
    """Pointer to a field attribute of the LabSep object
    
    """
    
    storage_name = target
    access_flag = flag
    
    @property
    def prop(self):
        at = getattr(self.labsetup,storage_name)
        return at[self.index]
    
    @prop.setter
    def prop(self,value):  
        at = getattr(self.labsetup,storage_name)
        if access_flag is not None:
            setattr(self, access_flag, True)
        at[self.index] = value
        
    return prop

def _labarray(name, target):
    """Pointer to a field array attribute of the LabSep object
    
    """
    
    storage_name = target
    
    @property
    def prop(self):
        at = getattr(self.labsetup,storage_name)
        return at[self.index,:]
    
    @prop.setter
    def prop(self,value):  
        at = getattr(self.labsetup,storage_name)  
        at[self.index,:] = value
        
    return prop


def _fieldprop(name, flag, sign):
    """Property returning field values over time
    
    """
    cmplx_sign = sign
    
    @property 
    def prop(self):
        if getattr(self.labsetup, flag):
            if cmplx_sign == 1:
                return self.get_field()
            elif cmplx_sign == -1:
                return numpy.conj(self.get_field())
            elif cmplx_sign == 0:
                fld = self.get_field()
                return (fld + numpy.conj(fld))/2.0
            else:
                raise Exception("Only signs of -1, 0 and 1 are allowed.")
        else:
            raise Exception("The property '"+name+"' is not initialited.")
    
    @prop.setter
    def prop(self, value):
        raise Exception("The property '"+name+"' is protected"+
                        " and cannot be set.")
        
    return prop        
    
        
def _get_example_lab():
    """Returns a LabSetup instance for doctests
    
    """
    lab = LabSetup(nopulses=3)
    
    return lab
    

class LabField():
    """Class representing electric field of a laser pulse defined in LabSetup


    Objects of this class are linked to their "mother" object of the LabSetup
    type. Their properties can be changed locally, with an effect on the
    LabSetup, or globally from the LabSetup with an effect on the EField
    objects.


    Examples
    --------
    
    Only the number of pulses has to be specified when LabSetup is created.
    
    >>> lab = LabSetup(nopulses=3)
    
    We can ask for a LabField object right away, even before field parameters
    are set.
    
    >>> lf = LabField(lab, 1)

    This object has all parameters "empty"
    
    >>> lf.pol
    array([ 0.,  0.,  0.])
    
    >>> lf.om
    0.0
    
    >>> lf.tc
    0.0
    
    >> lf.phi
    0.0
    
    The 'field' property, however, refuses to return values
    
    >>> print(lf.field)
    Traceback (most recent call last):
        ...
    Exception: The property 'field' is not initialited.
    
    Nor it can be set
    
    >>> lf.field = 10.0
    Traceback (most recent call last):
        ...
    Exception: The property 'field' is protected and cannot be set.
    
    The LabField properties will be initialited through the LabSetup object.
    The only rule to follow is that arrival times of the pulses have to be
    specified before the pulse shape.
        
    >>> lab.set_pulse_arrival_times([0.0, 0.0, 100.0])
    
    >>> time = TimeAxis(-500.0, 1000, 1.0, atype="complete")    
    >>> pulse2 = dict(ptype="Gaussian", FWHM=150, amplitude=1.0)  
    >>> params = (pulse2, pulse2, pulse2)    
    >>> lab.set_pulse_shapes(time, params) 

    Everything else can be set before we ask for the field's time dependence.
    The LabField object can be created even 

    >>> lab.set_pulse_polarizations(pulse_polarizations=(X,X,X),
    ...                             detection_polarization=X)
    >>> lab.set_pulse_frequencies([1.0, 1.0, 1.0])
    >>> lab.set_pulse_phases([0.0, 1.0, 0.0])    

    
    >>> lf = LabField(lab, 2)
    >>> print(lf.get_phase() == lab.phases[2])
    True
    
    >>> lf.set_phase(3.14)
    >>> print(lf.get_phase() == lab.phases[2])
    True
    
    >>> lab.phases[2] = 6.28
    >>> print(lf.get_phase() == lab.phases[2])
    True
    
    >>> print(lf.get_center() == lab.pulse_centers[2])
    True
    
    >>> lf.set_center(12.0)
    >>> print(lab.pulse_centers[2])
    12.0
    
    >>> print(lf.get_frequency() == lab.omega[2])
    True
    
    >>> lf.set_frequency(12.0)
    >>> print(lab.omega[2])
    12.0
    
    >>> lf.get_polarization()
    array([ 1.,  0.,  0.])
    
    >>> lab.e[2,:] = [0.0, 1.0, 0.0]
    >>> lf.get_polarization()
    array([ 0.,  1.,  0.])
    
    >>> lf.set_polarization([0.0, 0.0, 1.0])
    >>> lab.e[2,:]
    array([ 0.,  0.,  1.])
    
    
    # we also have some quick access attributes
    
    >>> print(lf.phi)
    6.28
    
    >>> lf.phi = 1.2
    >>> lf.phi
    1.2

    >>> print(lf.phi == lab.phases[2])
    True

    >>> print(lf.tc)
    12.0
    
    >>> lf._center_changed
    False
    
    >>> lf.tc = 10.0
    >>> lf.tc
    10.0
    
    >> lf._center_changed
    True

    >>> print(lf.tc == lab.pulse_centers[2])
    True


    >>> print(lf.om)
    12.0
    
    >>> lf.om = 10.0
    >>> lf.om
    10.0

    >>> print(lf.om == lab.omega[2])
    True


    >>> lf.pol
    array([ 0.,  0.,  1.])
    
    >>> lab.e[2,:] = [0.0, 1.0, 0.0]
    >>> lf.pol
    array([ 0.,  1.,  0.])
    
    >>> lf.pol = [0.0, 0.0, 1.0]
    >>> lab.e[2,:]
    array([ 0.,  0.,  1.])

 
    Most importantly, we can access the field values 
    
    >>> fld = lf.field
    >>> fld.shape
    (1000,)
    
    and this property cannot be directly changed.
    >>> lf.field = 10.0
    Traceback (most recent call last):
        ...
    Exception: The property 'field' is protected and cannot be set.

    """
    
    
    phi = _labattr("phi","phases")
    delay_phi = _labattr("delay_phi","delay_phases")
    tc = _labattr("tc", "pulse_centers", flag="_center_changed")
    om = _labattr("om", "omega")
    pol = _labarray("pol", "e")
    field_p = _fieldprop("field_p","_field_set", 1)
    field_m = _fieldprop("field_p","_field_set", -1)
    field = _fieldprop("field","_field_set", 0)
    
    def __init__(self, labsetup, k):
        
        self.labsetup = labsetup
        self.index = k
        self._center_changed = False

        # delay phase 
        self.set_delay_phase(self.tc)

    
    def get_phase(self):
        """Returns the phase of the pulse
        
        """
        return self.labsetup.phases[self.index]
    
    
    def get_delay_phase(self):
        """Returns the phase caused by the pulse delay
        
        """
        return self.labsetup.delay_phases[self.index]


    def get_total_phase(self):
        
        return (self.labsetup.delay_phases[self.index]
                +self.labsetup.phases[self.index])
    
    
    def set_phase(self, val):
        """Sets the phase of the pulse
        
        Parameters
        ----------
        
        val : float
        
        
        """
        self.labsetup.phases[self.index] = val
        # FIXME: The pulse field has to be updated


    def set_delay_phase(self, val):
        """Returns the phase caused by the pulse delay
        
        """
        raise Exception("Setting delay phases independently is not allowed")
    

    def get_center(self):
        """Returns the pulse center time
        
        """
        return self.labsetup.pulse_centers[self.index]
    

    def set_center(self, val):
        """Sets the pulse center time
        
        val: float
            The center of the pulse
            
        """
        self.labsetup.pulse_centers[self.index] = val
      
        # calculate phase shift associated with the delay
        self.set_delay_phase(val)
        #om = self.labsetup.omega[self.index]
        #phi = val*om
        #self.labsetup.delay_phases[self.index] = phi 
        
        # reset the pulse shapes
        self.labsetup.reset_pulse_shape()

        self._center_changed = False

        
    def set_delay_phase(self, val):
        """Calculate and set delay phases
        
        """
        om = self.labsetup.omega[self.index]
        phi = val*om

        #print("Setting delay phase of", phi, "val=",val, "om=", om)
        self.labsetup.delay_phases[self.index] = phi 

        
    def get_frequency(self):
        return self.labsetup.omega[self.index]

    def set_frequency(self, val):
        self.labsetup.omega[self.index] = val
        
    def get_polarization(self):
        return self.labsetup.e[self.index,:]
    
    def set_polarization(self, pol):
        self.labsetup.e[self.index,:] = pol

    def get_fwhm(self):
        return self.labsetup.saved_params[self.index]["FWHM"]
        
        
    def get_field(self, time=None, sign=1):
        """Returns the electric field of the pulses
        
        """
        if self._center_changed:
            # recalculate pulses
            self.labsetup.reset_pulse_shape()
            # FIXME: might require reseting the phase too!!!
        
        if time is None:
            tt = self.labsetup.timeaxis.data
            env = self.labsetup.pulse_t[self.index].data
            om = self.om
            phi = self.phi
            delay_phi = self.delay_phi
            #print(phi, delay_phi)
            fld = env*numpy.exp(-1j*sign*om*tt)*numpy.exp(1j*sign*phi) \
                     *numpy.exp(1j*sign*delay_phi)
            return fld
        
        else:
            return self.labsetup.pulse_t[self.index].at(time)


    def get_time_axis(self):
        return self.labsetup.timeaxis 
    
    
    def set_rwa(self, om):
        
        self.labsetup.set_rwa(om)
        

    def restore_rwa(self):
        
        self.labsetup.restore_rwa()
        

    def as_spline_function(self):
        """Returns the spline represention of this field
        
        """
        df = DFunction(self.labsetup.timeaxis, self.get_field())
        return df.as_spline_function()
        
            
    def get_pulse_envelop(self, tt):
        """Returns the envelop values
        
        
        """

        #tma = self.timeaxis
        if self.labsetup.saved_params[self.index]["ptype"] == "Gaussian":
            fwhm = self.labsetup.saved_params[self.index]["FWHM"]
            amp = self.labsetup.saved_params[self.index]["amplitude"]
            
            #tc = self.pulse_centers[k_p]
            
            # normalized Gaussian mupliplied by amplitude
            lfc = 4.0*numpy.log(2.0)
            pi = numpy.pi
            val = (2.0/fwhm)*numpy.sqrt(numpy.log(2.0)/pi) \
            *amp*numpy.exp(-lfc*(tt/fwhm)**2) 
            
            return val
        
        else:
            
            raise Exception()
            
            
    def get_pulse_envelop_function(self):
        """Return a function to be called later
        
        
        """

        if self.labsetup.saved_params[self.index]["ptype"] == "Gaussian":
            fwhm = self.labsetup.saved_params[self.index]["FWHM"]
            amp = self.labsetup.saved_params[self.index]["amplitude"]
            lfc = 4.0*numpy.log(2.0)
            pi = numpy.pi

            def env(tt):

                val = (2.0/fwhm)*numpy.sqrt(numpy.log(2.0)/pi) \
                    *amp*numpy.exp(-lfc*(tt/fwhm)**2) 
            
                return val

            return env
        
        else:

            raise Exception()  



