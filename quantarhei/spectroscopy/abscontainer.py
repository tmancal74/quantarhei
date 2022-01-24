# -*- coding: utf-8 -*-

from ..core.saveable import Saveable


class AbsSpectrumContainer(Saveable):
    
    def __init__(self, axis=None):

        self.axis = axis
        self.count = 0
        self.spectra = {}
        
    def set_axis(self, axis):
        self.axis = axis
        
    def set_spectrum(self, spect, tag=None):
        """Stores absorption spectrum 
        
        Checks compatibility of its frequency axis
        
        """
        frq = spect.axis
        
        if self.axis is None:
            self.axis = frq
            
        if self.axis.is_equal_to(frq):
            if tag is None:
                tag1 = str(self.count)
            else:
                tag1 = str(tag)
            self.spectra[tag1] = spect
            self.count += 1
        else:
            raise Exception("Incompatible time axis (equal axis required)")
            
            
    def get_spectrum(self, tag):
        """Returns spectrum corresponing to time t2
        
        Checks if the time t2 is present in the t2axis
        
        """  
        if not isinstance(tag, str):
            tag = str(tag)
            
        if tag in self.spectra.keys():
            return self.spectra[tag]     
        else:
            raise Exception("Unknown spectrum")

        
    def get_spectra(self):
        """Returns a list or tuple of the calculated spectra
        
        """
        
        ven = [value for (key, value) in sorted(self.spectra.items())]
        return ven        
