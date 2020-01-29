# -*- coding: utf-8 -*-

import numpy

from .twod import TwoDSpectrum
from .twodcontainer import TwoDResponseContainer
from .twodcalculator import TwoDResponseCalculator
from .mocktwodcalculator import MockTwoDResponseCalculator
from ..core.dfunction import DFunction

import matplotlib.pyplot as plt


class PumpProbeSpectrum(DFunction):
    """Class representing pump-probe spectra
    
    
    
    """
    
    #data = None
    
    def __init__(self):
        self.t2 = -1.0
        
    def set_axis(self, axis):
        self.xaxis = axis
        self.axis = self.xaxis
        
    def set_data(self, data):
        self.data = data


    def get_PumpProbeSpectrum(self):
        """Returns self
        
        This method is here to override the one inherented from TwoDSpectrum
        
        """
        return self


    def set_t2(self, t2):
        """Sets the t2 (waiting time) of the spectrum
        
        
        """
        self.t2 = t2


    def get_t2(self):
        """Returns the t2 (waiting time) of the spectrum
        
        """
        return self.t2



class PumpProbeSpectrumContainer(TwoDResponseContainer):
    """Container for a set of pump-probe spectra
    
    """
    def __init__(self, t2axis=None):
        
        self.t2axis = t2axis
        self.spectra = {}
        self.itype = None
        
    def plot(self):
        
        plt.clf()
        spctr = self.get_spectra()
        for sp in spctr:
            plt.plot(sp.xaxis.data,sp.data)
            
    def set_spectrum(self, spec, tag):
        self.spectra[tag] = spec
        
    def amax(self):
        mxs = []
        for s in self.get_spectra():       
            spect = numpy.real(s.data)   
            mx = numpy.amax(spect)
            mxs.append(mx)
        return numpy.amax(numpy.array(mxs))   
     
    def amin(self):
        mxs = []
        for s in self.get_spectra():       
            spect = numpy.real(s.data)   
            mx = numpy.amin(spect)
            mxs.append(mx)
        return numpy.amin(numpy.array(mxs))   



    def make_movie(self, filename, axis=None,
                   cmap=None, vmax=None, vmin=None,
                   frate=20, dpi=100, start=None, end=None,
                   show_states=None, progressbar=False):
        
        import matplotlib.pyplot as plt
        import matplotlib.animation as manimation
        
        FFMpegWriter = manimation.writers["ffmpeg"]

        metadata = dict(title="Test Movie", artist='Matplotlib',
                comment='Movie support!')
        writer = FFMpegWriter(fps=frate, metadata=metadata)
        
        fig = plt.figure() 
        
        spctr = self.get_spectra()
        l = len(spctr)
        last_t2 = spctr[l-1].get_t2()
        first_t2 = spctr[0].get_t2()
        
        if vmax is None:
            mx = self.amax()
        else:
            mx = vmax
            
        if vmin is None:
            mn = self.amin()
        else:
            mn = vmin
            
        mxabs = max(numpy.abs(mx), numpy.abs(mn))
        mx = mx+0.05*mxabs
        mn = mn-0.05*mxabs
        
        if start is None:
            start = first_t2
        if end is None:
            end = last_t2
        
                
        with writer.saving(fig, filename, dpi):  
            k = 0
            # Initial call to print 0% progress
            sp2write = self.get_spectra(start=start, end=end)
            l = len(sp2write)
            if progressbar:
                self._printProgressBar(0, l, prefix = 'Progress:',
                                       suffix = 'Complete', length = 50)
            for sp in self.get_spectra(start=start, end=end):
                # FIXME: this does not work as it should yet
                sp.plot(show=False, fig=fig, axis=axis,
                        label="T="+str(sp.get_t2())+"fs") #, vmax=mx, vmin=mn,
                          #)
                writer.grab_frame()
                if progressbar:
                    self._printProgressBar(k + 1, l, prefix = 'Progress:',
                                           suffix = 'Complete', length = 50)
                
                k += 1


class PumpProbeSpectrumCalculator(TwoDResponseCalculator):
    
    def __init__(self, t1axis, t2axis, t3axis):
        pass
    


def calculate_from_2D(twod):
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
    #tddata = numpy.real(twod.d__data)
    tddata = numpy.real(twod.data)
    
    # integration over omega_1 axis
    ppdata = -numpy.sum(tddata,axis=1)*dw
    
    # setting pump-probe data
    pp.set_data(ppdata)
    pp.set_axis(twod.yaxis)
    
    return pp
    


class MockPumpProbeSpectrumCalculator(MockTwoDResponseCalculator):
    """Effective line shape pump-probe spectrum calculator
    
    
    """
    
    
    def calculate_all_system(self, sys, H, eUt, lab):
        """Calculates all spectra corresponding to a specified t2axis
        
        """
        
        temporary_fix = True
        
        if temporary_fix:
            # calculation via 2D spectrum
            tcont = super().calculate_all_system(sys, H, eUt, lab)
            return tcont.get_PumpProbeSpectrumContainer()
      
    
    def calculate_one_system(self, t2, sys, H, eUt, lab):
        """Calculates one spectru corresponding to a specified t2 time
        
        """
        
        temporary_fix = True
        
        if temporary_fix:
            # calculation via 2D spectrum
            tcont = super().calculate_one_system(t2, sys, H, eUt, lab)
            return tcont.get_PumpProbeSpectrum()


        