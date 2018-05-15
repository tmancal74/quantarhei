# -*- coding: utf-8 -*-
"""




"""

import numpy


from ..utils import Integer
from ..utils.vectors import X

class labsetup:
    """Laboratory set-up for non-linear spectroscopy
    
    
    
    """
    
    number_of_pulses = Integer("number_of_pulses")
    
    def __init__(self, nopulses = 3):
        
        self.number_of_pulses = nopulses
    
        self.M4 = numpy.array([[4.0, -1.0, -1.0],
                               [-1.0, 4.0, -1.0],
                               [-1.0,-1.0,  4.0]])/30.0
                        

    def set_pulses(self, params):
        pass

    
    def set_polarizations(self, pulse_polarizations=(X, X, X), 
                         detection_polarization=X):
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
            text = "pulse_polarizations requires "+str(self.order)+" values"
            raise Exception(text)
            
        self.detection_polarization = detection_polarization

     
    def get_pulse_polarizations(self):
        return (self.e[0,:],self.e[1,:],self.e[2,:])

        
    def get_detection_polarization(self):
        return self.e[3,:]


    def set_pulse_spectra(self):
        pass
    
    
    def set_pulse_envelop(self, k, t):
        pass
    
    
    def get_pulse_spectrum(self, k, omega):
        pass
    
    