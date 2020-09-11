# -*- coding: utf-8 -*-
"""Class holding a set of TwoDSpectra



    Class Details
    -------------

"""
import numbers

#import h5py
#import matplotlib.pyplot as plt  
import numpy

from ..core.time import TimeAxis
from ..core.valueaxis import ValueAxis
from ..core.frequency import FrequencyAxis
from ..core.dfunction import DFunction
#from .twod2 import TwoDResponse
from .twod import TwoDSpectrum

from ..core.managers import Manager, energy_units

#from ..core.managers import energy_units
from .. import COMPLEX
from .. import REAL

from ..core.saveable import Saveable

from .. import part_REAL, part_IMAGINARY, part_COMPLEX, part_ABS
from .. import signal_TOTL #, signal_REPH, signal_NONR


class TwoDResponseContainer(Saveable):
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
        
        self.axis = None
        
        self.itype = None
        self.index = 0
        self.tags = []
        
        if self.keep_pathways:
            raise Exception("Container keeping pathways not available yet")
            
        self.spectra = {}
        self._which = None
        
        if t2axis is not None:
            self.use_indexing_type(itype=t2axis)
        
        
    def use_indexing_type(self, itype):
        """Sets the type of indices used to identify spectra
        
        Parameters
        ----------
        
        itype : string, ValueAxis, TimeAxis, FrequencyAxis
            Type of indexig. If itype is a string, it should have values
            either 'integer' or 'string' in which case the specra will be 
            stored by integer index or by string (as in dictionary). If 
            itype is a ValueAxis (TimeAxis, FrequencyAxis), spectra will
            be indexed by the values in the axis object.
        
        """
        
        if isinstance(itype, str):
            if itype == "integer":
                self.itype = "integer"
            elif itype == "string":
                self.itype = "string"
            else:
                raise Exception("Unknown indexing type")
        elif isinstance(itype, ValueAxis):
            if isinstance(itype, TimeAxis):
                self.itype = "TimeAxis"
                self.axis = itype
                # This axis must FFT in the "standard" way 
                # -- it must of the "complete" type 
                self.axis.atype = "complete" 
            elif isinstance(itype, FrequencyAxis):
                self.itype = "FrequencyAxis"
                self.axis = itype
            else:
                self.itype = "ValueAxis"
                self.axis = itype
        else:
            raise Exception("Unknown indexing type")
        

    def set_spectrum(self, spect, tag=None):
        """Stores spectrum with a tag (time, index, etc.)
        
        Stores the spectrum according to present indexing scheme
        
        Parameters
        ----------
        
        spect : TwoDSpectrum
            Object holding the spectrum; when not tag is specified for a
            spectrum which as its t2 time set, the tag is set to t2 time.
            
        tag : {int, string, ValuesAxis, TimeAxis, FrequencyAxis}
            Tag which will be used for retrieval of the spectrum from 
            the container.
        
        """
        
        
        if self.itype == "integer":
            
            if tag is None:
                self.spectra[self.index] = spect
                self.index += 1
                return self.index
            else:
                if isinstance(tag, numbers.Integral):
                    self.spectra[tag] = spect
                else:
                    raise Exception("The spectrum has to be tagged by an integer")
                return tag
        
        elif self.itype in ["ValueAxis", "TimeAxis", "FrequencyAxis"]:
            
            if tag is None:
                # we will read the spectrum intrinsic t2 time and set it as tag
                if spect.t2 >= 0.0:
                    tag = spect.t2
                    
            if tag is not None:
                if tag in self.axis.data:
                    self.spectra[tag] = spect
                    self.tags.append(tag)
                    self.index += 1
                else:
                    raise Exception("Tag not compatible with the ValueAxis")
            else:
                raise Exception("No tag specified, and spectrum"
                                +" does not have t2 time set"
                                +" - cannot store spectrum")
            return self.index
        
        elif self.itype == "string":
            if tag is not None:
                stag = str(tag)
                self.spectra[stag] = spect
                self.tags.append(stag)
                self.index += 1
            else:
                raise Exception("No tag specified - cannot store spectrum")
            return self.index

        else:
            
            raise Exception("Unknown type of indexing")    


    def _lousy_equal(self, x1, x2, dx, frac=0.25):
        """Equals up to fraction of dx
        
        This function returns True if x1 is closer to x2 than `frac` of 
        a specified interval. In addition it saves the value of x2 to which
        x1 is equal in the attribute _which of the present class.
        
        
        """
        if abs(x1-x2) < dx*frac: 
            self._which = x2
            return True
        
        self._which = None
        return False


    def get_spectrum_by_index(self, indx):
        """Returns spectrum by integet index
        
        The integer index is assigned to all spectra in the order they were
        saved to the container. They can be retrieved in this order
        
        Parameters
        ----------
        
        indx : int
            Index of the spectrum to be retrieved
            
        """
        
        if self.itype == "integer":
            
            return self.get_spectrum(indx)
        
        else:

            return self.spectra[self.tags[indx]]


    def get_response(self, tag):
        """Same as get_spectrum, but the name more sense for a response container
        
        """
        return self.get_spectrum(tag)


    def get_spectrum(self, tag):
        """Returns spectrum corresponing to time t2
        
        Checks if the time t2 is present in the t2axis
        
        Parameters
        ----------
        
        t2 : float
            Waiting time for which spectrum should be returned
            
            
        """        
        if self.itype in ["integer"]:
            
            return self.spectra[tag]
            
        elif self.itype in ["string"]:
            
            return self.spectra[str(tag)]

        elif self.itype in ["ValueAxis", "TimeAxis", "FrequencyAxis"]:
            
            with energy_units("int"):
                if any(self._lousy_equal(tag, li, self.axis.step) 
                   for li in self.axis.data):
    
                    try:
                        return self.spectra[self._which]     
                    except KeyError:
                        print(self.spectra)
                        raise Exception()      
                else:
                    raise Exception("Tag not compatible with the ValueAxis")
            
        else:
            
            raise Exception("Unknown type of indexing")

    
    def set_data_flag(self, flag):
        """Sets data flag for all spectra in the container
        
        
        """
        
        for tag in self.spectra:
            
            sp = self.spectra[tag]
            sp.set_data_flag(flag)


    def get_TwoDSpectrumContainer(self, stype=signal_TOTL):
        """Returns a container with specific spectra
        
        """
        if self.itype in ["ValueAxis", "TimeAxis", "FrequencyAxis"]:
            axis = self.axis.deepcopy()
        
            cont = TwoDSpectrumContainer(axis)
        
            for val in self.axis.data:
                sp = self.get_spectrum(val)
                nsp = sp.get_TwoDSpectrum(dtype=stype)
                cont.set_spectrum(nsp, tag=val)
                
            return cont
        
        else:
            
            raise Exception("")
            

    def get_nearest(self, val):
        
        if self.itype == "FrequencyAxis":
            #print(Manager().current_units["frequency"])
            #print(Manager().current_units["energy"])
            nval = Manager().convert_energy_2_internal_u(val)
            # get tags and convert them to numbers
            ntags = numpy.zeros(len(self.tags), dtype=REAL)
            k = 0
            for stag in self.tags:
                #print(stag)
                ntag = float(stag)
                ntags[k] = ntag
                k += 1
            dtags = numpy.abs(ntags - nval)
            imin = numpy.argmin(dtags)
            #print("Returning spectrum at: ", 
            #      Manager().convert_energy_2_current_u(self.tags[imin]))
            return self.spectra[self.tags[imin]], imin
            

    def length(self):
        """Returns the length of the container
        
        
        """
        return len(self.spectra.keys())


    def get_spectra(self, start=None, end=None):
        """Returns a list of the calculated spectra
        
        Returns all spectra or an interval of spectra when `start` and `end`
        are specified
        
        Parameters
        ----------
        
        start : int
            Index of the first spectrum to be returned
            
        end : int
            Index of the last spectrum to be returned

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
        ttc = []
        ii = 0
        for sp in self.get_spectra():
            if k == 0:
                pp = sp.get_PumpProbeSpectrum()
                ppc.append(pp)
                ttc.append(self.axis.data[ii])
            k += 1
            if k > skip:
                k = 0
            ii += 1
            
        length = len(ppc)
        start = ppc[0].get_t2()
        step = ppc[1].get_t2()-start

        naxis = TimeAxis(start,length,step)        
        ppcont = PumpProbeSpectrumContainer(t2axis=naxis)
        ppcont.itype = self.itype

        ii = 0
        for sp in ppc:
            tt = ttc[ii]
            ppcont.set_spectrum(sp, tt)
            ii += 1
            
        return ppcont


    def get_integrated_area_evolution(self, times, area, dpart=part_REAL):
        """Returns the integrated area of the 2D spectra as a function of their index

        """
        vals = numpy.zeros(times.length, dtype=COMPLEX)
        k = 0
  
        # this only acts on Frequency axis
        with energy_units("int"):
            tms = times.data      
  
        for t2 in tms:
            
            sp = self.get_spectrum(t2)
            vals[k] = sp.get_area_integral(area, dpart=part_REAL)
            k +=1
            
        return DFunction(times, vals)        

    
    def get_point_evolution(self, x, y, times):
        """Tracks an evolution of a single point on the 2D spectrum
        
        
        Parameters
        ----------
        
        x : float
            x coordinate in the 2D spectrum (usually omega_1 axis)

        y : float
            y coordinate in the 2D spectrum (usually omega_3 axis)
            
        times : ValueAxis
            Times (usually waiting t_2 times) in which spectra are taken
            
        """

        vals = numpy.zeros(times.length, dtype=COMPLEX)
        k = 0
  
        # this only acts on Frequency axis
        with energy_units("int"):
            tms = times.data      
  
        for t2 in tms:
            
            sp = self.get_spectrum(t2)
            vals[k] = sp.get_value_at(x, y)
            k +=1
            
        return DFunction(times, vals)

    
    def global_fit_exponential(self, guess=None):
        """Global fit of the data with exponentials
        
        
        """
        from scipy.optimize import least_squares
        from functools import partial
        
        if guess is None:
            guess = [1.0, 1.0/100.0, 0.0]
 
        _exp_2D_fcion = partial(_exp_2D_data0, times=self.axis.data, cont=self)
            
        params = least_squares(_exp_2D_fcion, guess)           
        
        return params

    
    
    def fft(self, ffttype="complex-positive", window=None, offset=0.0,
            dtype=None, dpart=part_COMPLEX, tag=None):
        """Fourier transform in t2 time
        
        This method performs FFT on the container data determined by the
        value of the `dtype` argument. The new container is created and 
        the storage resolution of its components is set `off`. This means
        that the container and its spectra have no idea about what data they
        store. Even when plotting the spectra, one has to set plotting of
        the `total` spectrum.
        
        Parameters
        ----------
        
        ffttype : string
            Specifies the type Fourier transform we perform
            
        window : DFunction
            Windowing function for the data. Default is None
        
        """
        
        if dtype is None:
            raise Exception("Type of the data for FFT has to be specified")
            
        if self.itype not in ["ValueAxis", "TimeAxis", "FrequencyAxis"]:
            raise Exception("FFT cannot be performed for"+
                            " this type of indexing")

        # even when no window function is supplied, we create one with
        # all elements equal to one
        if window is None:
            winfce = DFunction(self.axis, 
                               numpy.ones(self.axis.length, dtype=REAL))
        else:
            winfce = window
             
        if isinstance(self.axis, TimeAxis):
            # restrict the time axis by the off-set
            tlist = []
            for tt in self.axis.data:
                if tt >= offset:
                    tlist.append(tt)
            if len(tlist) > 1:
                dt = tlist[1]-tlist[0]
                Nt = len(tlist)
                t0 = tlist[0]
                # effective time axis for fft
                eff_axis = TimeAxis(t0, Nt, dt, atype="complete") 
            else:
                raise Exception("Offset too large")
        else:
            eff_axis = self.axis
            
        # put all data into one array
        
        #raise Exception()
        self.set_data_flag(dtype)
        
        tags = eff_axis.data #self.axis.data
        Nos = self.length()

        if len(tags) <= Nos:
            tag1 = eff_axis.data[0] #self.axis.data[0]
            sp1 = self.get_spectrum(tag1)
            
            N1, N2 = sp1.d__data.shape
            data = numpy.zeros((N1, N2, len(tags)), dtype=sp1.d__data.dtype)


            for k_n in range(len(tags)):
                tag = eff_axis.data[k_n] #self.axis.data[k_n]
                
                spect = self.get_spectrum(tag)
                spect.set_data_flag(dtype)
                if dpart == part_COMPLEX:
                    data[:,:,k_n] = spect.d__data
                elif dpart == part_REAL:
                    data[:,:,k_n] = numpy.real(spect.d__data)
                elif dpart == part_IMAGINARY:
                    data[:,:,k_n] = numpy.imag(spect.d__data)
                elif dpart == part_ABS:
                    data[:,:,k_n] = numpy.abs(spect.d__data)
                    
            
        else:
            raise Exception("Number of spectra not consistent"+
                            " with ValueAxis object")
            
        #
        # FFT of the axis
        #
        axis = eff_axis # = self.axis            
        
        if isinstance(axis, TimeAxis):
            axis.shift_to_zero()
            new_axis = axis.get_FrequencyAxis()
            
        elif isinstance(axis, FrequencyAxis):
            new_axis = axis.get_TimeAxis()
            
        else: 
            # this must be ValueAxis

            ftaxdata = (2.0*numpy.pi)*numpy.fft.fftfreq(axis.length,
                                                        axis.step)
            ftaxdata = numpy.fft.fftshift(ftaxdata)
            
            start = ftaxdata[0]
            length = len(ftaxdata)
            step = ftaxdata[1]-ftaxdata[0]

            new_axis = ValueAxis(start, length, step)            

        
        #
        # FFT of the data
        #
        
        # window function
        ftdata = numpy.zeros(data.shape, dtype=data.dtype)
        Nwin = len(winfce.data)
        Ndat = data.shape[2]
        for i_n in range(data.shape[0]):
            for j_n in range(data.shape[1]):
                ftdata[i_n,j_n,:] = data[i_n,j_n,:]*winfce.data[Nwin-Ndat:Nwin]
                
        ftdata = numpy.fft.ifft(ftdata, axis=2)
        ftdata = numpy.fft.fftshift(ftdata, axes=2)
        
        # save it to a new container
        new_container = TwoDSpectrumContainer()
        new_container.use_indexing_type(new_axis)

        for k_n in range(Ndat):
            tag = new_axis.data[k_n]
            spect = TwoDSpectrum()
            spect.set_axis_1(sp1.xaxis)
            spect.set_axis_3(sp1.yaxis)
            
            spect.set_data(ftdata[:, :, k_n], dtype=signal_TOTL)

            new_container.set_spectrum(spect, tag=tag)
        
        return new_container

          
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

        
    def trimall_to(self, window=None):
        """Trims all spectra in the container

        Parameters
        ----------
        
        window: list of floats
            Window, specified by four float number, to which all spectra
            in the container should be trimmed
            
        """
        if window is not None:
            axes = window
            for s in self.get_spectra():
                s.trim_to(window=axes)

    # FIXME: this is still legacy version        
    def amax(self, spart=part_REAL):
        """Returns maximum amplitude of the spectra in the container
        
        """
        

        mxs = []
        for s in self.get_spectra():       
            #spect2D = numpy.real(s.d__data)
            if spart == part_REAL:
                spect2D = numpy.real(s.data)
            elif spart == part_IMAGINARY:
                spect2D = numpy.imag(s.data)
            elif spart == part_ABS:
                spect2D = numpy.abs(s.data)
            else:
                raise Exception("Unknow part of the spectrum:", spart)
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
                   stype=signal_TOTL, spart=part_REAL, 
                   cmap=None, 
                   Npos_contours=10,
                   vmax=None,
                   xlabel=None,
                   ylabel=None,
                   axis_label_font=None,
                   start=None, end=None,
                   frate=20, dpi=100, 
                   show_states=None, 
                   show_states_func=None,
                   label=None,
                   label_func=None,
                   text_loc=None,
                   progressbar=False, 
                   use_t2=True, 
                   title="Quantarhei movie",
                   comment="Created with Quantarhei"):
        """Creates a movie out of the spectra in the container
        
        
        Parameters
        ----------
        
        Npos_contours : int
            Nomber of positive value contours in the plot
            
        
        
        """
        
        
        import matplotlib.pyplot as plt
        import matplotlib.animation as manimation
        
        FFMpegWriter = manimation.writers["ffmpeg"]

        metadata = dict(title=title, artist='Quantarhei',
                comment=comment)
        writer = FFMpegWriter(fps=frate, metadata=metadata)
        
        fig = plt.figure() 
        
        spctr = self.get_spectra()
        l = len(spctr)
        
        if use_t2 and (spctr[0].get_t2() < 0.0):
            #print("Warning: switching off usage of t2"
            #      +" information obtained from the spectrum object (t2 < 0)")
            use_t2 = False
        
        if use_t2:
            last_t2 = spctr[l-1].get_t2()
            first_t2 = spctr[0].get_t2()
        
            if start is None:
                start = first_t2
            if end is None:
                end = last_t2

        if vmax is None:
            mx = self.amax(spart=spart)
        else:
            mx = vmax        
                
        with writer.saving(fig, filename, dpi):  
            k = 0
            # Initial call to print 0% progress
            if use_t2:
                sp2write = self.get_spectra(start=start, end=end)
            else:
                sp2write = self.get_spectra()
            l = len(sp2write)
            
            if progressbar:
                self._printProgressBar(0, l, prefix = 'Progress:',
                                       suffix = 'Complete', length = 50)
                                
            for sp in sp2write: #self.get_spectra(start=start, end=end):
                
                if label_func is not None:
                    (label, text_loc) = label_func(sp)
                if show_states_func is not None:
                    show_states = show_states_func(sp)

                sp.plot(fig=fig, window=window, cmap=cmap, vmax=mx, 
                        Npos_contours=Npos_contours,
                        stype=stype,spart=spart,
                        show_states=show_states,
                        xlabel=xlabel, ylabel=ylabel,
                        axis_label_font=axis_label_font,
                        label=label, text_loc=text_loc) #"T="+str(sp.get_t2())+"fs")
                writer.grab_frame()
                if progressbar:
                    self._printProgressBar(k + 1, l, prefix = 'Progress:',
                                           suffix = 'Complete', length = 50)
                
                k += 1
#                if k == 20:
#                    return


def _exp_2D_data0(params, times=None, cont=None):
    """Returns a residue between time dependent matrix data and a matrix
       multipled by a sum of exponentials
    
    """
    
    if times is None:
        raise Exception("Times have to be supplied")
    if cont is None:
        raise Exception("Spectra container has to be supplied")
        
    np = len(params)
    
    data0 = cont.get_spectrum(0.0)
    
    residues = numpy.zeros((times.shape[0], data0.data.shape[0],
                            data0.data.shape[1]), dtype=data0.dtype)
    
    N = times.shape[0]*data0.data.shape[0]*data0.data.shape[0]
    
    nexp = int((np - 1)/2)
    
    ii = 0
    for t in times:

        ret = 0.0
        kp = 0

        for kk in range(nexp):
            ret += params[kp]*numpy.exp(-params[kp+1]*t)
            kp += 2
        ret += params[kp]
        
        spect = cont.get_spectrum(t)
        
        residues[ii,:,:] = numpy.abs(spect.data - ret*data0.data)
        
        ii += 1
        
    return residues.reshape(N)






class TwoDSpectrumContainer(TwoDResponseContainer):
    
    def __init__(self, t2axis=None, dtype=signal_TOTL):
        
        self.t2axis = t2axis

        self.axis = None
        
        self.itype = None
        self.index = 0
        self.tags = []
        self.dtype = dtype
            
        self.spectra = {}
        self._which = None
        
        if t2axis is not None:
            self.use_indexing_type(itype=t2axis)
            

    def set_data_flag(self, flag):
        """Sets data flag for all spectra in the container
        
        
        """
        
        if flag != self.dtype:
            raise Exception("Cannot change spectra type")

            
    def get_TwoDSpectrumContainer(self, stype=signal_TOTL):
        """Returns a container with specific spectra
        
        """
        if stype == self.dtype:
            return self
        else:
            raise Exception("Cannot change spectra type in this container")

           
    # FIXME: This needs to be reimplemented
    def get_PumpProbeSpectrumContainer(self, skip=0):
        """Converts this container into PumpProbeSpectrumContainer
        
        """
        
        if self.dtype == signal_TOTL:
        
            from .pumpprobe import PumpProbeSpectrumContainer
            
            k = 0
            ppc = []
            ttc = []
            ii = 0
            for sp in self.get_spectra():
                if k == 0:
                    pp = sp.get_PumpProbeSpectrum()
                    ppc.append(pp)
                    ttc.append(self.axis.data[ii])
                k += 1
                if k > skip:
                    k = 0
                ii += 1
                
            length = len(ppc)
            start = ppc[0].get_t2()
            step = ppc[1].get_t2()-start
    
            naxis = TimeAxis(start,length,step)        
            ppcont = PumpProbeSpectrumContainer(t2axis=naxis)
    
            ii = 0
            for sp in ppc:
                tt = ttc[ii]
                ppcont.set_spectrum(sp, tt)
                ii += 1
                
            return ppcont  
        
        else:
            
            raise Exception("Cannot calculate Pump-probe from 2D"+
                            " spectra of type"+self.dtype)
            
    def normalize2(self, norm=1.0, each=False, dpart=part_REAL):
        """Normalize the whole container of spectra
        
        Normalization of the whole container so that the maximum
        value of the spectrum is equal to the `norm`.
        
        Parameters
        ----------
        
        norm : float
            Value to which we normalize the spectra
            
        each: bool
            If False, we normalize the global maximum of the container, i.e.
            a maximum accross whole spectra. if True, each spectrum is
            normalized individually against its own maximum.
            
        dpart: string
            Part of the spectrum from which the maximum is calculated, it
            can be part_REAL, part_IMAGINARY or part_ABS. The values of these
            constants are defined in the highest level of namespace in
            Quantarhei.
        
        """
        
        nsp = len(self.spectra)
        mxs = numpy.zeros(nsp, dtype=REAL)
        ii = 0
        for tag in self.spectra.keys():
            sp = self.get_spectrum(tag)
            mxs[ii] = sp.get_max_value(dpart=dpart)
            ii += 1
            
        mx = numpy.amax(mxs)
            
        nmax = [mx]
        for tag in self.spectra.keys():
            sp = self.get_spectrum(tag)
            if each:
                nmax = [sp.get_max_value(dpart=dpart)]
            sp.normalize2(norm, dpart=dpart, nmax=nmax)
            
        return nmax

            
    def fft(self, ffttype="complex-positive", window=None, offset=0.0,
            dtype=None, dpart=part_COMPLEX, tag=None):
        """Fourier transform in t2 time
        
        This method performs FFT on the container data determined by the
        value of the `dtype` argument. The new container is created and 
        the storage resolution of its components is set `off`. This means
        that the container and its spectra have no idea about what data they
        store. Even when plotting the spectra, one has to set plotting of
        the `total` spectrum.
        
        Parameters
        ----------
        
        ffttype : string
            Specifies the type Fourier transform we perform
            
        window : DFunction
            Windowing function for the data. Default is None
        
        """
        
        if dtype is not None:
            if dtype != self.dtype:
                raise Exception("Cannot change spectra type"+
                                " in TwoDSpectrumContainer")
                
        if self.itype not in ["ValueAxis", "TimeAxis", "FrequencyAxis"]:
            raise Exception("FFT cannot be performed for"+
                            " this type of indexing")

        # even when no window function is supplied, we create one with
        # all elements equal to one
        if window is None:
            winfce = DFunction(self.axis, 
                               numpy.ones(self.axis.length, dtype=REAL))
        else:
            winfce = window
             
        if isinstance(self.axis, TimeAxis):
            # restrict the time axis by the off-set
            tlist = []
            for tt in self.axis.data:
                if tt >= offset:
                    tlist.append(tt)
            if len(tlist) > 1:
                dt = tlist[1]-tlist[0]
                Nt = len(tlist)
                t0 = tlist[0]
                # effective time axis for fft
                eff_axis = TimeAxis(t0, Nt, dt, atype="complete") 
            else:
                raise Exception("Offset too large")
        else:
            eff_axis = self.axis
            
        # put all data into one array
        
        #raise Exception()
        
        tags = eff_axis.data #self.axis.data
        Nos = self.length()

        if len(tags) <= Nos:
            tag1 = eff_axis.data[0] #self.axis.data[0]
            sp1 = self.get_spectrum(tag1)
            
            N1, N2 = sp1.data.shape
            data = numpy.zeros((N1, N2, len(tags)), dtype=sp1.data.dtype)


            for k_n in range(len(tags)):
                tag = eff_axis.data[k_n] #self.axis.data[k_n]
                
                spect = self.get_spectrum(tag)

                if dpart == part_COMPLEX:
                    data[:,:,k_n] = spect.data
                elif dpart == part_REAL:
                    data[:,:,k_n] = numpy.real(spect.data)
                elif dpart == part_IMAGINARY:
                    data[:,:,k_n] = numpy.imag(spect.data)
                elif dpart == part_ABS:
                    data[:,:,k_n] = numpy.abs(spect.data)
                    
            
        else:
            raise Exception("Number of spectra not consistent"+
                            " with ValueAxis object")
            
        #
        # FFT of the axis
        #
        axis = eff_axis # = self.axis            
        
        if isinstance(axis, TimeAxis):
            axis.shift_to_zero()
            new_axis = axis.get_FrequencyAxis()
            
        elif isinstance(axis, FrequencyAxis):
            new_axis = axis.get_TimeAxis()
            
        else: 
            # this must be ValueAxis

            ftaxdata = (2.0*numpy.pi)*numpy.fft.fftfreq(axis.length,
                                                        axis.step)
            ftaxdata = numpy.fft.fftshift(ftaxdata)
            
            start = ftaxdata[0]
            length = len(ftaxdata)
            step = ftaxdata[1]-ftaxdata[0]

            new_axis = ValueAxis(start, length, step)            

        
        #
        # FFT of the data
        #
        
        # window function
        ftdata = numpy.zeros(data.shape, dtype=data.dtype)
        Nwin = len(winfce.data)
        Ndat = data.shape[2]
        for i_n in range(data.shape[0]):
            for j_n in range(data.shape[1]):
                ftdata[i_n,j_n,:] = data[i_n,j_n,:]*winfce.data[Nwin-Ndat:Nwin]
                
        ftdata = numpy.fft.ifft(ftdata, axis=2)
        ftdata = numpy.fft.fftshift(ftdata, axes=2)
        
        # save it to a new container
        new_container = TwoDSpectrumContainer()
        new_container.use_indexing_type(new_axis)

        for k_n in range(Ndat):
            tag = new_axis.data[k_n]
            spect = TwoDSpectrum()
            spect.set_axis_1(sp1.xaxis)
            spect.set_axis_3(sp1.yaxis)
            
            spect.set_data(ftdata[:, :, k_n], dtype=self.dtype)

            new_container.set_spectrum(spect, tag=tag)
        
        return new_container


    def unitedir(self, dname):
        
        clist = self.loaddir(dname)

        n = len(clist)

        cont = clist[1]
        
        for k in range(1, n):
            
            ctn = clist[k+1]
            
            for tag in ctn.spectra:
                cont.set_spectrum(ctn.spectra[tag], tag=tag)
                
        return cont

        