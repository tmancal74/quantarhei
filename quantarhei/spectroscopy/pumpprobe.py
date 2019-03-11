# -*- coding: utf-8 -*-

import numpy

from .twod2 import TwoDSpectrum
from .twodcontainer import TwoDSpectrumContainer
from .twodcalculator import TwoDSpectrumCalculator
from .mocktwodcalculator import MockTwoDSpectrumCalculator
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


#    def plot(self, fig=None, axis=None, title=None, text=None, label=None,
#             xlabel=None,
#             ylabel=None, label_font=None, text_font=None,
#             text_loc=[0.05,0.1], fontsize="20",
#             vmax=None, vmin=None, color=None, show=False):
#        """Plots 2D spectrum
#        
#        
#        """
#        if axis is not None:
#            axs = list(axis)
#            
#            if len(axs) == 2:
#                axs.append(numpy.amin(self.data))
#                axs.append(numpy.amax(self.data))
#                
#            if len(axs) != 4:
#                raise Exception("Wrong axis specification")
#        else:
#            axs = [self.xaxis.min, self.xaxis.max,
#                   numpy.amin(self.data),numpy.amax(self.data)]
#            
#        if fig is None:
#            fig, ax = plt.subplots(1,1)
#        else:
#            fig.clear()
#            fig.add_subplot(1,1,1)
#            ax = fig.axes[0]
#            
#        if color is not None:
#            plt.plot(self.xaxis.data,self.data, color)
#        else:
#            plt.plot(self.xaxis.data,self.data)
#
#        if vmax is not None:
#            axs[3] = vmax
#            
#        if vmin is not None:
#            axs[2] = vmin
#            
#        if axs is not None:
#            plt.axis(axs)
#
#        if title is not None:
#            plt.title(title)
#
#        if text is not None:
#            if text_font is not None:
#                plt.text(text[0],text[1],
#                     text[2], fontdict=text_font)
#            else:
#                plt.text(text[0], text[1], text[2])
#
#        if label_font is not None:
#            font = label_font
#        else:
#            font={'size':20}        
#
#        #
#        # Label
#        #
#        levo = axs[0]
#        prvo = axs[1]
#        dole = axs[2]
#        hore = axs[3]
#        pos = text_loc
#        if label is not None:
#            label = label    
#            ax.text((prvo-levo)*pos[0]+levo,
#                (hore-dole)*pos[1]+dole,
#                label,
#                fontsize=str(fontsize))
#
#        if xlabel is None:
#            xl = ""
#        if ylabel is None:
#            yl = ""    
#
#        if xlabel is not None:
#            xl = xlabel
#        if ylabel is not None:
#            yl = ylabel
#
#        if xl is not None:
#            plt.xlabel(xl, **font)
#        if xl is not None:
#            plt.ylabel(yl, **font)
#
#        if show:
#            plt.show()


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
                sp.plot(fig=fig, axis=axis, vmax=mx, vmin=mn,
                        label="T="+str(sp.get_t2())+"fs")
                writer.grab_frame()
                if progressbar:
                    self._printProgressBar(k + 1, l, prefix = 'Progress:',
                                           suffix = 'Complete', length = 50)
                
                k += 1


class PumpProbeSpectrumCalculator(TwoDSpectrumCalculator):
    
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
    tddata = numpy.real(twod.d__data)
    
    # integration over omega_1 axis
    ppdata = -numpy.sum(tddata,axis=1)*dw
    
    # setting pump-probe data
    pp.set_data(ppdata)
    pp.set_axis(twod.yaxis)
    
    return pp
    


class MockPumpProbeSpectrumCalculator(MockTwoDSpectrumCalculator):
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


        