# -*- coding: utf-8 -*-

import numpy

from .twod2 import TwoDSpectrum
from .twod2 import TwoDSpectrumContainer
from .twod2 import TwoDSpectrumCalculator

import matplotlib.pyplot as plt

class PumpProbeSpectrum(TwoDSpectrum):
    """Class representing pump-probe spectra
    
    
    
    """
    
    def __init__(self):
        self.t2 = -1.0
        
    def set_axis(self, axis):
        self.xaxis = axis
        
    def set_data(self, data):
        self.data = data

    def get_PumpProbeSpectrum(self):
        """Returns self
        
        This method is here to override the one inherented from TwoDSpectrum
        
        """
        return self
    
    def plot(self, fig=None):
        """Plots 2D spectrum
        
        
        """
        if fig is None:
            fig, ax = plt.subplots(1,1)
        else:
            fig.clear()
            fig.add_subplot(1,1,1)
            ax = fig.axes[0]
            
        plt.plot(self.xaxis.data,self.data)

        
    

class PumpProbeSpectrumContainer(TwoDSpectrumContainer):
    """Container for a set of pump-probe spectra
    
    """
    def __init__(self, t2axis=None):
        
        self.t2axis = t2axis
        self.spectra = {}
        
    def plot(self):
        
        plt.clf()
        spctr = self.get_spectra()
        for sp in spctr:
            plt.plot(sp.xaxis.data,sp.data)
            
            


class PumpProbeSpectrumCalculator(TwoDSpectrumCalculator):
    
    def __init__(self):
        pass
    
    def calculate_from_2D(self, twod):
        """Calculates pump-probe spectrum from 2D spectrum
        
        Calculates pump-probe spectrum from 2D spectrum
        using projection theorem
        
        Parameters
        ----------
        
        twod: TwoDSpectrum
            2D spectrum from which pump-probe will be calculated
            
        """
        
        pp = PumpProbeSpectrum()
        
        # waiting time
        t2 = twod.get_t2()
        pp.set_t2(t2)
        
        # time step for integration
        xaxis = twod.xaxis
        dw = xaxis.step
        
        # real part of the total 2D spectrum
        tddata = numpy.real(twod.reph2D + twod.nonr2D)
        
        # integration over omega_1 axis
        ppdata = -numpy.sum(tddata,axis=1)*dw
        
        # setting pump-probe data
        pp.set_data(ppdata)
        pp.set_axis(twod.yaxis)
        
        return pp
        