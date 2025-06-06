# -*- coding: utf-8 -*-
"""Two-dimensional Fourier Transform Spectrum and its calculation


#
# Includes past Fluorescence detected 2D code
#


"""
#import h5py
import matplotlib.pyplot as plt  
import numpy

from ..core.time import TimeAxis
from ..core.frequency import FrequencyAxis
from ..builders.aggregates import Aggregate
from ..builders.molecules import Molecule
from ..core.managers import eigenbasis_of
from ..core.managers import energy_units
from ..qm.propagators.poppropagator import PopulationPropagator 
from ..core.units import convert
from ..core.saveable import Saveable

from ..utils import derived_type

import quantarhei as qr

import time

try:

    import aceto.nr3td as nr3td            
    from aceto.lab_settings import lab_settings
    from aceto.band_system import band_system            
    _have_aceto = True
    
except:
    #
    # FIXME: There should be an optional warning and a fall back onto
    # quantarhei.implementations.aceto module
    #
    #raise Exception("Aceto not available")
    #from ..implementations.aceto import nr3td            
    _have_aceto = False 


class DFunction2:
    """Descrete two-dimensional function 
    
    Note: This will be later moved into the core package of Quantarhei
    It handles basic saving and loading operatoins of a two-dimensional
    discrete function
    """
    def __init__(self, x=None, y=None, z=None):
        pass
    
    def save(self, filename):
        pass
    
    def load(self, filename):
        pass
    
    def plot(self,**kwargs):
        pass    


class TwoDSpectrumBase(DFunction2):
    """Basic class of a two-dimensional spectrum
    
    
    """
    
    # spectral types
    stypes = ["rephasing","non-rephasing","nonrephasing", "total"]
    
    # spectral parts
    sparts = ["real","imaginary"]
    
    # to keep Liouville pathways separate?
    keep_pathways = False
    
    # to keep stypes separate?
    keep_stypes = True
    
    def __init__(self):
        super().__init__()
        
        self.data = None
        self.xaxis = None
        self.yaxis = None
            
        self.reph2D = None
        self.nonr2D = None


    def set_axis_1(self, axis):
        """Sets the x-axis of te spectrum (omega_1 axis)
        
        """
        self.xaxis = axis
        
    def set_axis_3(self, axis):
        """Sets the y-axis of te spectrum (omega_3 axis)
        
        """
        self.yaxis = axis
        
    def set_data(self, data, dtype="Tot"):
        if dtype == "Tot":
            self.data = data
            
        elif dtype == "Reph":
            
            self.reph2D = data
            
        elif dtype == "Nonr":
            
            self.nonr2D = data
            
        else:
            
            raise Exception("Unknow type of data: "+dtype)

    def add_data(self, data, dtype="Tot"):
        if dtype == "Tot":
            
            if self.data is None:
                self.data = numpy.zeros(data.shape, dtype=data.dtype)
            self.data += data
            
        elif dtype == "Reph":
            if self.reph2D is None:
                self.reph2D = numpy.zeros(data.shape, dtype=data.dtype)
            self.reph2D += data                
            
        elif dtype == "Nonr":
            if self.nonr2D is None:
                self.nonr2D = numpy.zeros(data.shape, dtype=data.dtype)                
            self.nonr2D += data
                
            
        else:
            
            raise Exception("Unknow type of data: "+dtype)

        
    def save(self, filename):
        super().save(filename)
    
    def load(self, filename):
        super().load(filename)
        
    
    
class TwoDSpectrum(TwoDSpectrumBase, Saveable):
    """This class represents a single 2D spectrum
    
    Methods
    -------
    
    plot(fig=None, window=None, stype="total", spart="real",
         vmax=None, vmin_ratio=0.5, 
         colorbar=True, colorbar_loc="right",
         cmap=None, Npos_contours=10,
         show_states=None,
         text_loc=[0.05,0.9], fontsize="20", label=None)
    
    
    """
    
    def __init__(self, keep_pathways=False, keep_stypes=True):
        self.keep_pathways = keep_pathways
        self.keep_stypes = keep_stypes
        self.t2 = -1.0
        super().__init__()
    
    def set_t2(self, t2):
        """Sets the t2 (waiting time) of the spectrum
        
        
        """
        self.t2 = t2
        
    def get_t2(self):
        """Returns the t2 (waiting time) of the spectrum
        
        """
        return self.t2
    
    
    def get_value_at(self, x, y):
        """Returns value of the spectrum at a given coordinate
        
        """
        (ix, dist) = self.xaxis.locate(x)
        (iy, dist) = self.yaxis.locate(y)    
        
        return numpy.real(self.reph2D[ix,iy]+self.nonr2D[ix,iy])
    
    def get_max_value(self):
        
        return numpy.amax(numpy.real(self.reph2D+self.nonr2D))
            
    
    def devide_by(self, val):
        """Devides the total spectrum by a value
        
        """
        self.reph2D = self.reph2D/val
        self.nonr2D = self.nonr2D/val  


    def get_PumpProbeSpectrum(self):
        """Returns a PumpProbeSpectrum corresponding to the 2D spectrum
        
        """
        from .pumpprobe import PumpProbeSpectrumCalculator
        ppc = PumpProbeSpectrumCalculator()
        return ppc.calculate_from_2D(self)
    
    
        
    
    def plot(self, fig=None, window=None, stype="total", spart="real",
             vmax=None, vmin_ratio=0.5, 
             colorbar=True, colorbar_loc="right",
             cmap=None, Npos_contours=10,
             show_states=None,
             text_loc=[0.05,0.9], fontsize="20", label=None, zero_contour=True,
             plot_units=r'cm$^{-1}$'):
        """Plots the 2D spectrum
        
        Parameters
        ----------
        
        fig : matplotlib.figure
            Figure into which plotting will be done. This is used e.g. when
            making a movie using moview writter (may be obsolete)
            
        window : list
            Specifies the plotted window in current energy units. When axes
            are x and y, the window is specified as window=[x_min,x_max,y_min,y_max]
            
        stype : {"total", "rephasing", "non-rephasing"}
            type of the spectrum 
            
            
        """
        
        if stype == "total":
            if (self.reph2D is not None) and (self.nonr2D is not None):
                spect2D = self.reph2D + self.nonr2D 
            elif self.reph2D is not None:
                spect2D = self.reph2D 
            elif self.nonr2D is not None:
                spect2D = self.nonr2D
                
            
        elif stype == "rephasing":
            spect2D = self.reph2D
        elif stype == "non-rephasing":
            spect2D = self.nonr2D            
        else:
            raise Exception("Undefined spectrum type"+stype)
        
        if spart == "real":
            spect2D = numpy.real(spect2D)
        elif spart == "imaginary":
            spect2D = numpy.imag(spect2D)
        elif spart == "abs":
            spect2D = numpy.abs(spect2D)
        else:
            raise Exception("Undefined part of the spectrum: "+spart)
         
            
        if window is not None: 
            axis = window
            w1_min = axis[0]
            w1_max = axis[1]
            w3_min = axis[2]
            w3_max = axis[3]

            (i1_min, dist) = self.xaxis.locate(w1_min)
            (i1_max, dist) = self.xaxis.locate(w1_max)

            (i3_min, dist) = self.yaxis.locate(w3_min)
            (i3_max, dist) = self.yaxis.locate(w3_max)   
            
        else:
            i1_min = 0
            i1_max = self.xaxis.length
            i3_min = 0
            i3_max = self.yaxis.length
            
    
        #
        # Plotting with given units on axes
        #
  
        realout = spect2D[i3_min:i3_max,i1_min:i1_max]
    
        if fig is None:
            fig, ax = plt.subplots(1,1)
        else:
            fig.clear()
            fig.add_subplot(1,1,1)
            ax = fig.axes[0]
            
        if cmap is None:
            cmap = plt.cm.rainbow
            
        if vmax is None:
            vmax = numpy.amax(realout)

        vmin = numpy.amin(realout)
        if vmin < -vmax*vmin_ratio:
            vmax = -vmin
        else:
            vmin = -vmax*vmin_ratio
        
        Npos = Npos_contours
        poslevels = [i*vmax/Npos for i in range(1, Npos)]
        neglevels = [-i*vmax/Npos for i in range(Npos,1,-1)]
        
        levo = self.xaxis.data[i1_min]
        prvo = self.xaxis.data[i1_max-1]
        dole = self.yaxis.data[i3_min]
        hore = self.yaxis.data[i3_max-1]
        
        cm = plt.imshow(realout, extent=[self.xaxis.data[i1_min],
                                    self.xaxis.data[i1_max-1],
                                    self.yaxis.data[i3_min],
                                    self.yaxis.data[i3_max-1]],
                   origin='lower', vmax=vmax, vmin=vmin,
                   interpolation='bilinear', cmap=cmap)  
        
        pos = text_loc
        
        # text
        if label is not None:
            label = label    
            ax.text((prvo-levo)*pos[0]+levo,
                (hore-dole)*pos[1]+dole,
                label,
                fontsize=str(fontsize))
        
        #axis labels
        plt.xlabel(r'$\omega_{1}$ ('+plot_units+')')
        plt.ylabel(r'$\omega_{3}$ ('+plot_units+')')
        # positive contours
        plt.contour(self.xaxis.data[i1_min:i1_max],
                     self.yaxis.data[i3_min:i3_max],
                     realout, levels=poslevels, colors="k",
                     linewidth=1)
        
        if zero_contour:
            # zero contour
            plt.contour(self.xaxis.data[i1_min:i1_max],
                         self.yaxis.data[i3_min:i3_max],
                         realout, levels=[0],colors="b",
                         linewidth=1)
        
        # negatove contours
        plt.contour(self.xaxis.data[i1_min:i1_max],
                     self.yaxis.data[i3_min:i3_max],
                     realout, levels=neglevels,colors="k",
                     linewidth=1)  
        
        
        if colorbar:
            plt.clim(vmin=vmin,vmax=vmax)
            fig.colorbar(cm)
            
        if show_states is not None:
            for en in show_states:  
                plt.plot([en,en],[dole,hore],'--k',linewidth=1.0)
                plt.plot([levo,prvo],[en,en],'--k',linewidth=1.0)
            

    def trim_to(self, window=None):
        """Trims the 2D spectrum to a specified region
        
        """
        if window is not None:
            axis = window
            w1_min = axis[0]
            w1_max = axis[1]
            w3_min = axis[2]
            w3_max = axis[3]

            (i1_min, dist) = self.xaxis.locate(w1_min)
            (i1_max, dist) = self.xaxis.locate(w1_max)

            (i3_min, dist) = self.yaxis.locate(w3_min)
            (i3_max, dist) = self.yaxis.locate(w3_max)    
            
            # create minimal off-set
            i1_min -=1
            i1_max +=1
            i3_min -=1
            i3_max +=1
            
            # reconstruct xaxis
            start_1 = self.xaxis.data[i1_min]
            length_1 = i1_max - i1_min
            step_1 = self.xaxis.step
            xaxis = FrequencyAxis(start_1,length_1,step_1, 
                                  atype=self.xaxis.atype,
                                  time_start=self.xaxis.time_start)
            self.xaxis = xaxis
            
            # reconstruct yaxis
            start_3 = self.yaxis.data[i3_min]
            length_3 = i3_max - i3_min
            step_3 = self.yaxis.step
            yaxis = FrequencyAxis(start_3,length_3,step_3, 
                                  atype=self.yaxis.atype,
                                  time_start=self.yaxis.time_start)                
            self.yaxis = yaxis            
            

            # reconstruct data
            if self.keep_stypes:
                reph2D = self.reph2D[i1_min:i1_max,i3_min:i3_max]
                self.reph2D = reph2D   
                nonr2D = self.nonr2D[i1_min:i1_max,i3_min:i3_max]
                self.nonr2D = nonr2D 
                
            else:
                data = self.data[i1_min:i1_max,i3_min:i3_max]
                self.data = data
                
        else:
            # some automatic trimming in the future
            pass
                
    def show(self):
        
        plt.show()
        
    def savefig(self, filename, dpi):
        
        plt.savefig(filename,dpi=dpi)
        
    def _create_root_group(self, start, name):
        return start.create_group(name)
    
    def _save_attributes(self,rt):
        rt.attrs.create("t2", self.t2)
        keeps = []
        if self.keep_pathways:
            keeps.append(1)
        else:
            keeps.append(0)
        if self.keep_stypes:
            keeps.append(1)
        else:
            keeps.append(0)
            
        rt.attrs.create("keeps",keeps)

    def _load_attributes(self,rt):
        self.t2 = rt.attrs["t2"]
        keeps = rt.attrs["keeps"]
        self.keep_pathways = (keeps[0] == 1)
        self.keep_stypes = (keeps[1] == 1)
                    
    def _save_data(self,rt):
        if self.keep_stypes:
            rt.create_dataset("reph2D",data=self.reph2D)
            rt.create_dataset("nonr2D",data=self.nonr2D)
        else:
            rt.create_dataset("data",data=self.data)

    def _load_data(self,rt):
        if self.keep_stypes:
            self.reph2D = numpy.array(rt["reph2D"])
            self.nonr2D = numpy.array(rt["nonr2D"])
        else:
            self.data = numpy.array(rt["data"]) 
            
    def _save_axis(self, rt, name, ax):
        axdir = rt.create_group(name)
        axdir.attrs.create("start",ax.start)
        axdir.attrs.create("length",ax.length)
        axdir.attrs.create("step",ax.step)
        #FIXME: atype and time_start

    def _load_axis(self, rt, name):
        axdir = rt[name]
        start = axdir.attrs["start"]
        length = axdir.attrs["length"]
        step = axdir.attrs["step"]
        return FrequencyAxis(start, length, step)   
    
#    def save(self, filename, units="int"):
#        """Saves the whole object into file
#        
#        
#        """
#        with h5py.File(filename,"w") as f:
#            rt = self._create_root_group(f,"spectrum")
#            self._save_attributes(rt)
#            self._save_data(rt)
#            with energy_units(units):
#                self._save_axis(rt, "xaxis", self.xaxis)
#                self._save_axis(rt, "yaxis", self.yaxis)
#            
#            
#            
#    def load(self, filename, units="int"):
#        """Loads the whole object from a file
#        
#        
#        """
#        with h5py.File(filename,"r") as f:
#            rt = f["spectrum"]
#            self._load_attributes(rt)
#            self._load_data(rt)
#            with energy_units(units):
#                self.xaxis = self._load_axis(rt, "xaxis")
#                self.yaxis = self._load_axis(rt, "yaxis")    
#                
#        

        
class TwoDSpectrumContainer(Saveable):
    """Class holding a set of TwoDSpectra
    
    Parameters
    ----------
    
    t2axis: TimeAxis
       object holding waiting times at which spectra are calculated
       
    keep_pathways: bool
       if set True, the container will keep all types of Liouville pathways
       stored separately
       
    keep_stypes: bool
       if se t True, the container will keep rephasing and non-rephasing 
       spectra stored separately
       
       
    """
    
    def __init__(self, t2axis=None, keep_pathways=False, keep_stypes=True):
        
        self.t2axis = t2axis
        self.keep_pathways = keep_pathways
        self.keep_stypes = keep_stypes
        
        
        if self.keep_pathways:
            raise Exception("Container keeping pathways not available yet")
            
        self.spectra = {}
        
        
    def set_spectrum(self, spect):
        """Stores spectrum for time t2
        
        Checks if the time t2 is present in the t2axis
        
        """
        t2 = spect.get_t2()
        if t2 in self.t2axis.data:
            self.spectra[t2] = spect
        else:
            raise Exception("Waiting time not compatible with the t2 axis")
            
            
    def get_spectrum(self, t2):
        """Returns spectrum corresponing to time t2
        
        Checks if the time t2 is present in the t2axis
        
        """        
        if t2 in self.t2axis.data:
            return self.spectra[t2]     
        else:
            raise Exception("Waiting time not compatible with the t2 axis")

    def length(self):
        return len(self.spectra.keys())
        
    def get_spectra(self, start=None, end=None):
        """Returns a list or tuple of the calculated spectra
        
        """
        
        ven = [value for (key, value) in sorted(self.spectra.items())]
        
        if (start is None) and (end is None): 
            return ven
        else:
            ven2 = []
            vkeys = [key for (key, value) in sorted(self.spectra.items())]
            for k in vkeys:
                if k >= start and k <= end:
                    ven2.append(self.spectra[k])
            return ven2
        
    def get_PumpProbeSpectrumContainer(self, skip=0):
        """Converts this container into PumpProbeSpectrumContainer
        
        """
        
        from .pumpprobe import PumpProbeSpectrumContainer
        
        k = 0
        ppc = []
        for sp in self.get_spectra():
            if k == 0:
                pp = sp.get_PumpProbeSpectrum()
                ppc.append(pp)
            k += 1
            if k > skip:
                k = 0
            
        length = len(ppc)
        start = ppc[0].get_t2()
        step = ppc[1].get_t2()-start

        naxis = TimeAxis(start,length,step)        
        ppcont = PumpProbeSpectrumContainer(t2axis=naxis)

        for sp in ppc:
            ppcont.set_spectrum(sp)
            
        return ppcont
    
    def get_point_evolution(self, x, y, times):
        """Tracks an evolution of a single point on the 2D spectrum
        
        """
        
        vals = numpy.zeros(times.length)
        k = 0
        for t2 in times.data:
            
            sp = self.get_spectrum(t2)
            vals[k] = sp.get_value_at(x,y)
            k +=1
            
        return vals
            
    def _create_root_group(self, start, name):
        return start.create_group(name)

    def _save_axis(self, rt, name, ax):
        axdir = rt.create_group(name)
        axdir.attrs.create("start",ax.start)
        axdir.attrs.create("length",ax.length)
        axdir.attrs.create("step",ax.step)

    def _load_axis(self, rt, name):
        axdir = rt[name]
        start = axdir.attrs["start"]
        length = axdir.attrs["length"]
        step = axdir.attrs["step"]
        return TimeAxis(start, length, step) 
        
#    def save(self, filename):
#        """Saves the whole object into file
#        
#        
#        """
#        with energy_units("int"):
#            with h5py.File(filename,"w") as f:
#                self._save_axis(f,"t2axis",self.t2axis)
#                rt = self._create_root_group(f, "spectra")            
#                for sp in self.get_spectra():
#                    t2 = sp.get_t2
#                    rgname = "spectrum_"+str(t2)
#                    srt = sp._create_root_group(rt,rgname)
#                    sp._save_attributes(srt)
#                    sp._save_data(srt)
#                    sp._save_axis(srt,"xaxis",sp.xaxis,)
#                    sp._save_axis(srt,"yaxis",sp.yaxis)
#            
#                
#            
#            
#            
#    
#    def load(self, filename):
#        """Loads the whole object from a file
#        
#        
#        """
#        with energy_units("int"):
#            with h5py.File(filename,"r") as f:
#                self.t2axis = self._load_axis(f, "t2axis")
#                rt = f["spectra"]
#                for key in rt.keys():
#                    sp = TwoDSpectrum()
#                    srt = rt[key]
#                    sp._load_attributes(srt)
#                    sp._load_data(srt)
#                    sp.xaxis = sp._load_axis(srt,"xaxis")
#                    sp.yaxis = sp._load_axis(srt,"yaxis")
#                    
#                    self.set_spectrum(sp)
                
    def trimall_to(self,window=None):
        """Trims all spectra in the container
        
        """
        if window is not None:
            axes = window
            for s in self.get_spectra():
                s.trim_to(window=axes)
           
    def amax(self, spart="real"):
        """Returns maximum amplitude of the spectra in the container
        
        """
        mxs = []
        for s in self.get_spectra():       
            spect2D = numpy.real(s.reph2D) + numpy.real(s.nonr2D)   
            mx = numpy.amax(spect2D)
            mxs.append(mx)
        return numpy.amax(numpy.array(mxs))

    # Print iterations progress
    def _printProgressBar(self, iteration, total, 
                          prefix = '', suffix = '', 
                          decimals = 1, length = 100,
                          fill='*'):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
            
        Based on: 
        https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
        """
#                          fill = '█'):
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
        # Print New Line on Complete
        if iteration == total: 
            print()
    
             
    def make_movie(self, filename, window=None,
                   stype="total", spart="real", 
                   cmap=None, Npos_contours=10,
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
        
        mx = self.amax()
        
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
                sp.plot(fig=fig, window=window, cmap=cmap, vmax=mx, 
                        Npos_contours=Npos_contours,
                        stype=stype,spart=spart,
                        show_states=show_states,
                        label="T="+str(sp.get_t2())+"fs")
                writer.grab_frame()
                if progressbar:
                    self._printProgressBar(k + 1, l, prefix = 'Progress:',
                                           suffix = 'Complete', length = 50)
                
                k += 1
#                if k == 20:
#                    return
            
            
class TwoDSpectrumCalculator:
    """Calculator of the 2D spectrum
    
    
    Enables setting up parameters of 2D spectrum calculation for later
    evaluation. The method `calculate` returns TwoDSpectrumContainer
    with a 2D spectrum.
    
    Parameters
    ----------
    twodtype: either 2DES for traditional 2D or F-2DES for FL-detected 2D
    for F-2DES: either gamma_factor=2 (for annihilation yield) or
    timegate[0,1000] for detection window, k_rel for relaxation and k_an for 
    annihilation rate (all in ps)
    
    """

    t1axis = derived_type("t1axis",TimeAxis)
    t2axis = derived_type("t2axis",TimeAxis)
    t3axis = derived_type("t3axis",TimeAxis)
    
    system = derived_type("system",[Molecule,Aggregate])
    
    def __init__(self, t1axis, t2axis, t3axis,
                 system=None,
                 dynamics="secular",
                 relaxation_tensor=None,
                 rate_matrix=None,
                 effective_hamiltonian=None,
                 twodtype="2DES",
                 gamma_factor=None,
                 Population_factors=None,
                 ee_states=None, f_states=None,
                 no_transfer=False):
            
            
        self.t1axis = t1axis
        self.t2axis = t2axis
        self.t3axis = t3axis
        self.twodtype = twodtype
        self.gamma_factor = gamma_factor
        self.Population_factors=Population_factors
        self.ee_states=ee_states
        self.f_states=f_states
        self.no_transfer=no_transfer
        
        #FIXME: check the compatibility of the axes 
        
        if system is not None:
            self.system = system
        
        #FIXME: properties to be protected
        self.dynamics = dynamics
        
        # unprotected properties
        self.data = None
        
        self._relaxation_tensor = None
        self._rate_matrix = None
        self._relaxation_hamiltonian = None
        self._has_relaxation_tensor = False
        if relaxation_tensor is not None:
            self._relaxation_tensor = relaxation_tensor
            self._has_relaxation_tensor = True
        if effective_hamiltonian is not None:
            self._relaxation_hamiltonian = effective_hamiltonian
        if rate_matrix is not None:
            self._rate_matrix = rate_matrix
            self._has_rate_matrix = True
            
        #self._have_aceto = False
        
        # after bootstrap information
        self.sys = None
        self.lab = None
        self.t1s = None
        self.t3s = None
        self.rmin = None
        self.rwa = None
        self.oa1 = None
        self.oa3 = None
        self.Uee = None
        self.Uc0 = None
        
        self.tc = 0
        
       
    def _vprint(self, string):
        """Prints a string if the self.verbose attribute is True
        
        """
        if self.verbose:
            print(string)
            
    def bootstrap(self,rwa=0.0, lab=None, verbose=False):
        """Sets up the environment for 2D calculation
        
        """


        self.verbose = verbose
    
    
        if True:
            
            # calculate 2D spectrum using aceto library

            ###############################################################################
            #
            # Create band_system from quantarhei classes
            #
            ###############################################################################
            
            if isinstance(self.system, Aggregate):
            
                pass
            
            else:
                
                raise Exception("Molecule 2D not implememted")
                
            agg = self.system
            
            #
            # hamiltonian and transition dipole moment operators
            #
            H = agg.get_Hamiltonian()
            D = agg.get_TransitionDipoleMoment()
            
           
            #
            # Construct band_system object
            #
            Nb = 3
            Ns = numpy.zeros(Nb, dtype=numpy.int)
            Ns[0] = 1
            Ns[1] = agg.nmono
            Ns[2] = Ns[1]*(Ns[1]-1)/2
            self.sys = band_system(Nb, Ns)
            self.Ns=Ns
            
            
            #
            # Set energies
            #
            en = numpy.zeros(self.sys.Ne, dtype=numpy.float64)
            #if True:
            with eigenbasis_of(H):
                for i in range(self.sys.Ne):
                    en[i] = H.data[i,i]
                self.sys.set_energies(en)
            
                #
                # Set transition dipole moments
                #
                dge_wr = D.data[0:Ns[0],Ns[0]:Ns[0]+Ns[1],:]
                def_wr = D.data[Ns[0]:Ns[0]+Ns[1],
                                (Ns[0]+Ns[1]):(Ns[0]+Ns[1]+Ns[2]),:]
            
                dge = numpy.zeros((3,Ns[0],Ns[1]), dtype=numpy.float64)
                deff = numpy.zeros((3,Ns[1],Ns[2]), dtype=numpy.float64)
                
                for i in range(3):
                    dge[i,:,:] = dge_wr[:,:,i]
                    deff[i,:,:] = def_wr[:,:,i]
                self.sys.set_dipoles(0,1,dge)
                self.sys.set_dipoles(1,2,deff)
            
            
            #
            # Relaxation rates
            #
            KK = agg.get_RedfieldRateMatrix()
            
            # relaxation rate in single exciton band
            Kr = KK.data[Ns[0]:Ns[0]+Ns[1],Ns[0]:Ns[0]+Ns[1]] #*10.0
            #print(1.0/Kr)
            
            self.sys.init_dephasing_rates()
            self.sys.set_relaxation_rates(1,Kr)
            
            
            #
            # Lineshape functions
            #
            sbi = agg.get_SystemBathInteraction()
            cfm = sbi.CC
            cfm.create_double_integral()
            
            
            #
            # Transformation matrices
            #
            SS = H.diagonalize()
            SS1 = SS[1:Ns[1]+1,1:Ns[1]+1]
            SS2 = SS[Ns[1]+1:,Ns[1]+1:]
            H.undiagonalize()
            
            self.sys.set_gofts(cfm._gofts)    # line shape functions
            self.sys.set_sitep(cfm.cpointer)  # pointer to sites
            self.sys.set_transcoef(1,SS1)  # matrix of transformation coefficients  
            self.sys.set_transcoef(2,SS2)  # matrix of transformation coefficients  

            #
            # Finding population evolution matrix
            #
            prop = PopulationPropagator(self.t1axis, Kr)
      #      Uee, Uc0 = prop.get_PropagationMatrix(self.t2axis,
      #                                            corrections=True)
            self.Uee, cor = prop.get_PropagationMatrix(self.t2axis,
                                                  corrections=3)

            # FIXME: Order of transfer is set by hand here - needs to be moved
            # to some reasonable place
            
            #Ucor = Uee
            self.Uc0 = cor[0]
            
            #for ko in range(No+1):
            #    print("Subtracting ", ko)
            #    Ucor[:,:,tc] -= cor[ko]

            #
            # define lab settings
            #
            if lab is None:
                self.lab = lab_settings(lab_settings.FOUR_WAVE_MIXING)
                X = numpy.array([1.0, 0.0, 0.0], dtype=numpy.float64)
                self.lab.set_laser_polarizations(X,X,X,X)
            else:
                self.lab = lab
            
            #
            # Other parameters
            #
            #dt = self.t1axis.step
            self.rmin = 0.0001
            self.t1s = self.t1axis.data 
            self.t3s = self.t3axis.data
            self.rwa = rwa



            atype = self.t1axis.atype
            self.t1axis.atype = 'complete'
            self.oa1 = self.t1axis.get_FrequencyAxis() 
            self.oa1.data += self.rwa
            self.oa1.start += self.rwa
            self.t1axis.atype = atype
            
            atype = self.t3axis.atype
            self.t3axis.atype = 'complete'
            self.oa3 = self.t3axis.get_FrequencyAxis() 
            self.oa3.data += self.rwa
            self.oa3.start += self.rwa
            self.t3axis.atype = atype
            
        else:
            
            raise Exception("So far, no 2D outside aceto")
            
        self.tc = 0
            
    def calculate_next(self):

        sone = self.calculate_one(self.tc)
        self.tc += 1
        return sone
    
        
    def calculate_one(self, tc):

        tt2 = self.t2axis.data[tc]        
        Nr1 = self.t1axis.length
        Nr3 = self.t3axis.length        
        #
        # Initialize response storage
        #
        resp_r = numpy.zeros((Nr1, Nr3), 
                             dtype=numpy.complex128, order='F')
        resp_n = numpy.zeros((Nr1, Nr3), 
                             dtype=numpy.complex128, order='F')
        #additional response storage for the ESAs
        respf_r = numpy.zeros((Nr1, Nr3), 
                             dtype=numpy.complex128, order='F')
        respf_n = numpy.zeros((Nr1, Nr3), 
                             dtype=numpy.complex128, order='F')
        #additional response storage for the ESAs with state indices
        respf_r_i = numpy.zeros((self.Ns[2],Nr1, Nr3), 
                             dtype=numpy.complex128, order='F')
        respf_n_i = numpy.zeros((self.Ns[2],Nr1, Nr3), 
                             dtype=numpy.complex128, order='F')

        # FIXME: on which axis we should be looking for it2 ??? 
        (it2, err) = self.t1axis.locate(tt2) 
        self._vprint("t2 = "+str(tt2)+"fs (it2 = "+str(it2)+")")
        #tc = it2
    
        #
        # calcute response
        #
        self._vprint("calculating response: ")

        t1 = time.time()
        
        self._vprint(" - ground state bleach")
        # GSB
        nr3td.nr3_r3g(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_r) 
        nr3td.nr3_r4g(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_n)
    
        self._vprint(" - stimulated emission")
        # SE
        nr3td.nr3_r1g(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_n)
        nr3td.nr3_r2g(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_r)
        
        self._vprint(" - excited state absorption")
        # ESA
#        nr3td.nr3_r1fs(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_r)
#        nr3td.nr3_r2fs(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_n)
        nr3td.nr3_r1fs(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, respf_r_i, plist=True)
        nr3td.nr3_r2fs(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, respf_n_i, plist=True)
        
        # Transfer
        #with current version of ACETO works only for 2DES#
        #because the transfer pathways don't return list of pathways#
        #easy to fix in Aceto if needed#
        
        if self.twodtype == "2DES":
            Utr = self.Uee[:,:,self.tc]-self.Uc0[:,:,self.tc] #-Uc1[:,:,tc]-Uc2[:,:,tc]
            self.sys.set_population_propagation_matrix(Utr) 
            
            if not self.no_transfer:        
                self._vprint(" - stimulated emission with transfer")    
                # SE
                nr3td.nr3_r1g_trans(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_n)
                nr3td.nr3_r2g_trans(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_r)
                
        #                # This contributes only when No > 0
        #                nr3td.nr3_r2g_trN(lab, sys, No, it2, t1s, t3s, rwa, rmin, resp_r)
        #                
            
                self._vprint(" - excited state absorption with transfer") 
                
                # ESA
        #        nr3td.nr3_r1fs_trans(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_r)
        #        nr3td.nr3_r2fs_trans(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, resp_n)
                nr3td.nr3_r1fs_trans(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, respf_r)
                nr3td.nr3_r2fs_trans(self.lab, self.sys, it2, self.t1s, self.t3s, self.rwa, self.rmin, respf_n)
        

        
        #Normal ESA is summed over all states
        for i in range(0,self.Ns[2]):
            respf_n+=respf_n_i[i,:,:]
            respf_r+=respf_r_i[i,:,:]
        
        if self.twodtype == "2DES":
            resp_n=resp_n+respf_n
            resp_r=resp_r+respf_r
        
        elif self.twodtype == "F-2DES":
            if self.gamma_factor is not None:
                resp_n=resp_n+(self.gamma_factor-1)*respf_n
                resp_r=resp_r+(self.gamma_factor-1)*respf_r   
                self._vprint("Using gamma of: "+str(self.gamma_factor))
            elif self.Population_factors is not None:
                Pe=self.Population_factors[0]
                Pee=self.Population_factors[1]
                Pf=self.Population_factors[2]
                resp_n=(resp_n-respf_n)*Pe
                resp_r=(resp_r-respf_r)*Pe
                #distinguishing the ee and f states - to be made more general, by construction
                if self.ee_states is None:
                    self.ee_states=range(0,1)
                if self.f_states is None:
                    self.f_states=range(1,self.Ns[2])
                for i in self.ee_states:
                        resp_n+=Pee*respf_n_i[i,:,:]
                        resp_r+=Pee*respf_r_i[i,:,:]
                for i in self.f_states:
                        resp_n+=Pf*respf_n_i[i,:,:]
                        resp_r+=Pf*respf_r_i[i,:,:]
                        
            else:
                raise Exception("Not enough parameters for F-2DES")
            
        else:
            raise Exception("Unknown type of 2DES")
            
        t2 = time.time()
        self._vprint("... calculated in "+str(t2-t1)+" sec")


        #
        # Calculate corresponding 2D spectrum
        #
        
        ftresp = numpy.fft.fft(resp_r,axis=1)
        ftresp = numpy.fft.ifft(ftresp,axis=0)
        reph2D = numpy.fft.fftshift(ftresp)
        
        ftresp = numpy.fft.ifft(resp_n,axis=1)
        ftresp = numpy.fft.ifft(ftresp,axis=0)*ftresp.shape[1]
        nonr2D = numpy.fft.fftshift(ftresp)


        onetwod = TwoDSpectrum()
        onetwod.set_axis_1(self.oa1)
        onetwod.set_axis_3(self.oa3)
        onetwod.set_data(reph2D, dtype="Reph")
        onetwod.set_data(nonr2D, dtype="Nonr")
        
        onetwod.set_t2(self.t2axis.data[tc])
        
        
        return onetwod
                
                
            
    def calculate(self):
        """Returns 2D spectrum
        
        Calculates and returns TwoDSpectrumContainer containing 2D spectrum
        based on the parameters specified in this object.
        
        
        """            
                       
        if _have_aceto:

            twods = TwoDSpectrumContainer(self.t2axis)
            
            teetoos = self.t2axis.data
            for tt2 in teetoos:

                onetwod = self.calculate_next()
                twods.set_spectrum(onetwod)   
            
            return twods
        
        else:
            
            # fall back on quantarhei's own implementation
        
            ret = TwoDSpectrumContainer()
            
        
        return ret
    
    
class MockTwoDSpectrumCalculator(TwoDSpectrumCalculator):
    """Calculator of the 2D spectrum from LiouvillePathway objects
    
    
    This class is used to represent LiouvillePatjway objects. Lineshape is
    Gaussian 
    
    """

    def __init__(self, t1axis, t3axis):
        t2axis = TimeAxis()
        super().__init__(t1axis, t2axis, t3axis)
        self.widthx = convert(300, "1/cm", "int")
        self.widthy = convert(300, "1/cm", "int")
        self.dephx = convert(300, "1/cm", "int")
        self.dephy = convert(300, "1/cm", "int")        

        
    def bootstrap(self,rwa=0.0, pathways=None, verbose=False, 
                  shape="Gaussian", all_positive=False):
        
        self.shape = shape
        self.all_positive = all_positive
        
        self.verbose = verbose
        self.rwa = rwa
        self.pathways = pathways

        atype = self.t1axis.atype
        self.t1axis.atype = 'complete'
        self.oa1 = self.t1axis.get_FrequencyAxis() 
        self.oa1.data += self.rwa
        self.oa1.start += self.rwa
        self.t1axis.atype = atype
        
        atype = self.t3axis.atype
        self.t3axis.atype = 'complete'
        self.oa3 = self.t3axis.get_FrequencyAxis() 
        self.oa3.data += self.rwa
        self.oa3.start += self.rwa
        self.t3axis.atype = atype        
        

    def set_width(self, val):
        self.widthx = val
        self.widthy = val
        
    def set_deph(self, val):
        self.dephx = val
        self.dephy = val

        
    def calculate(self):
        """Calculate the 2D spectrum for all pathways
        
        """
        
        onetwod = TwoDSpectrum()
        onetwod.set_axis_1(self.oa1)
        onetwod.set_axis_3(self.oa3)
        
        for pwy in self.pathways:
            
            data = self.calculate_pathway(pwy, shape=self.shape)
            
            if pwy.pathway_type == "R":
                onetwod.add_data(data, dtype="Reph")
            elif pwy.pathway_type == "NR":
                onetwod.add_data(data, dtype="Nonr")
            else:
                raise Exception("Unknown pathway type")

        onetwod.set_t2(0.0)    
            
        return onetwod


    def calculate_pathway(self, pathway, shape="Gaussian"):
        """Calculate the shape of a Liouville pathway
        
        """
 
        noe = 1+pathway.order+pathway.relax_order 
        
        cen1 = pathway.frequency[0]
        cen3 = pathway.frequency[noe-2]
        if self.all_positive:
            pref = numpy.abs(pathway.pref)
        else:
            pref = pathway.pref
            
        N1 = self.oa1.length
        N3 = self.oa3.length
        
        if pathway.widths[1] < 0.0:
            widthx = self.widthx
        else:
            widthx = pathway.widths[1]
            
        if pathway.widths[3] < 0.0:
            widthy = self.widthy
        else:
            widthy = pathway.widths[3]
            
        if pathway.dephs[1] < 0.0:
            dephx = self.dephx
        else:
            dephx = pathway.dephs[1]
            
        if pathway.widths[3] < 0.0:
            dephy = self.dephy
        else:
            dephy = pathway.dephs[3]
        
        #print(shape, widthx, widthy)
        
        if pathway.pathway_type == "R":

            reph2D = numpy.zeros((N1, N3), dtype=qr.COMPLEX)
            
            if shape == "Gaussian":
                oo3 = self.oa3.data[:]
                for i1 in range(N1):
                    o1 = -self.oa1.data[i1]                    
                        
                    reph2D[:, i1] = \
                    pref*numpy.exp(-((o1-cen1)/widthx)**2)\
                        *numpy.exp(-((oo3-cen3)/widthy)**2)
            
            elif shape == "Lorentzian":
                oo3 = self.oa3.data[:]
                for i1 in range(N1):
                    o1 = -self.oa1.data[i1]                    
                        
                    reph2D[:, i1] = \
                    pref*(dephx/((o1-cen1)**2 + dephx**2))\
                        *(dephy/((oo3-cen3)**2 + dephy**2))
                        
            else:
                raise Exception("Unknown line shape: "+shape)   
            
            return reph2D
            
        elif pathway.pathway_type == "NR":
           
            nonr2D = numpy.zeros((N1, N3), dtype=qr.COMPLEX)
            
            if shape == "Gaussian":
                oo3 = self.oa3.data[:]
                for i1 in range(N1):
                    o1 = self.oa1.data[i1]                    
                    
                    nonr2D[:, i1] = \
                    pref*numpy.exp(-((o1-cen1)/widthx)**2)\
                        *numpy.exp(-((oo3-cen3)/widthy)**2)
                        
            elif shape == "Lorentzian":
                oo3 = self.oa3.data[:]
                for i1 in range(N1):
                    o1 = self.oa1.data[i1]                    
                        
                    nonr2D[:, i1] = \
                    pref*(dephx/((o1-cen1)**2 + dephx**2))\
                        *(dephy/((oo3-cen3)**2 + dephy**2))

            else:
                raise Exception("Unknown line shape: "+shape)
            
            return nonr2D
        
            
        
        