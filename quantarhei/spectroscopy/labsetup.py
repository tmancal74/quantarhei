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


from ..utils import Integer
from ..utils.vectors import X
from ..core.time import TimeAxis
from ..core.frequency import FrequencyAxis
from ..core.dfunction import DFunction

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
        
        self.number_of_pulses = nopulses
    
        self.M4 = numpy.array([[4.0, -1.0, -1.0],
                               [-1.0, 4.0, -1.0],
                               [-1.0,-1.0,  4.0]])/30.0
    
        self.timeaxis = None
        self.freqaxis = None
                        
        self.F4eM4 = None
        self.e = None
        self.has_polarizations = False
        
        self.has_freqdomain = False
        self.has_timedomain = False

        # time or frequency
        self.axis_type = None

        self.pulse_t = [None]*nopulses
        self.pulse_f = [None]*nopulses
        
        self.omega = None
        
                        

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
                \\frac{2}{\\Delta}\\sqrt{\\frac{4\\ln(2)}{\\pi}}
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

        if len(params) == self.number_of_pulses:
        
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
                        
                        # normalized Gaussian mupliplied by amplitude
                        val = (2.0/fwhm)*numpy.sqrt(numpy.log(2.0)/3.14159) \
                        *amp*numpy.exp(-4.0*numpy.log(2.0)*(tma.data/fwhm)**2)
                        
                        self.pulse_t[k_p] = DFunction(tma, val)
                    
                    elif self.axis_type == "frequency":
                        #
                        # Frequency domain Gaussian pulse around 0.0 
                        # as a DFunction
                        #
                        fra = self.freqaxis
                        fwhm = par["FWHM"]
                        amp = par["amplitude"]
                        
                        # normalized Gaussian mupliplied by amplitude
                        val = (2.0/fwhm)*numpy.sqrt(numpy.log(2.0)/3.14159) \
                        *amp*numpy.exp(-4.0*numpy.log(2.0)*(fra.data/fwhm)**2)
                        
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
        """Converts pulse information from frequency domain to time domain
        
        
        
        Examples
        --------
        
        >>> import quantarhei as qr
        >>> import matplotlib.pyplot as plt
        >>> lab = LabSetup()
        >>> freq = qr.FrequencyAxis(-100, 200, 1.0) # atype="complete" is default
        >>> pulse = dict(ptype="Gaussian", FWHM=20, amplitude=1.0)
        >>> params = (pulse, pulse, pulse)
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
        >>> lab.set_pulse_shapes(time, params)
        >>> dfc = lab.get_pulse_envelop(1, [-50.0, -30.0, 2.0, 30.0])
        >>> print(dfc)
        [  1.41569269e-05   1.95716182e-03   3.09310793e-02   1.95716182e-03]
        
        .. plot::
            :include-source:
            
            import quantarhei as qr
            import matplotlib.pyplot as plt
            
            lab = qr.LabSetup()
            time = qr.TimeAxis(-500.0, 1000, 1.0, atype="complete")
            pulse2 = dict(ptype="Gaussian", FWHM=150.0, amplitude=1.0)
            params = (pulse2, pulse2, pulse2)
            lab.set_pulse_shapes(time, params)
            
            pls = lab.pulse_t[2]
            lab2 = qr.LabSetup()
            
            pulse1 = dict(ptype="numeric", function=pls)
            params = (pulse1, pulse1, pulse1)
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
            lab.set_pulse_shapes(freq, params)
            
            pls = lab.pulse_f[2]
            lab2 = qr.LabSetup()
            
            pulse1 = dict(ptype="numeric", function=pls)
            params = (pulse1, pulse1, pulse1)
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
        [1.0, 2.0, 1.0]
        
        Situation which throws an exception

        >>> lab = LabSetup()
        >>> lab.set_pulse_frequencies([1.0, 2.0, 1.0, 6.0])       
        Traceback (most recent call last):
            ...
        Exception: Wrong number of frequencies: 3 required
        
        """
        
        # FIXME: energe unit control has to be in place
        if len(omegas) == self.number_of_pulses:
            
            self.omega = omegas
            
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
    
 
class labsetup(LabSetup):
    pass



class EField():
    """Class representing electric field of a laser pulse
    
    
    """
    
    def __init__(self, time, omega=0.0, polar=[1.0, 0.0, 0.0], 
                 ftype="Gaussian", params=None):
        
        self.time = time
        # FIXME: energy conversion
        self.omega = omega
        self.pol = polar
        
        if ftype=="Gaussian":
            
            self.ftype=ftype
            self.Emax = params["Emax"]
            # FIXME: energy conversion
            self.fwhm = params["fwhm"]
            self.tc = params["tc"]
            
            self.envelop = self.Emax*numpy.sqrt(numpy.log(2.0)/numpy.pi)\
                     *numpy.exp(-numpy.log(2.0)*((self.time.data
                                -self.tc)/self.fwhm)**2) \
                     /self.fwhm


    def subtract_omega(self, om):
        """
        
        """
        self.saved_omega = self.omega
        self.omega -= om


    def restore_omega(self):
        """
        
        """
        self.omega = self.saved_omega


    def field_i(self, sign, i):
        """Field at index i
        
        """
        if sign is None:
            return self.envelop[i]*numpy.cos(self.omega*self.time.data[i])
        
        if sign == 1:
            return 0.5*self.envelop[i]* \
                   numpy.exp(1j*self.omega*self.time.data[i])        
        elif sign == -1:
            return 0.5*numpy.conj(self.envelop[i])* \
                   numpy.exp(-1j*self.omega*self.time.data[i]) 
        else:
            raise Exception("Unknown field component")
    
    
    def field(self, sign=None):
        """Field in an array
        
        """
        if sign is None:
            return self.envelop*numpy.cos(self.omega*self.time.data)
        
        if sign == 1:
            return 0.5*self.envelop*numpy.exp(1j*self.omega*self.time.data)
        elif sign == -1:
            return 0.5*numpy.conj(self.envelop) \
                   *numpy.exp(-1j*self.omega*self.time.data)
        else:
            raise Exception("Unknown field component")
                   



