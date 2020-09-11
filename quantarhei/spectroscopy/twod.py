# -*- coding: utf-8 -*-
import warnings

import numpy
import matplotlib.pyplot as plt

from ..core.datasaveable import DataSaveable
from ..core.saveable import Saveable
from ..core.dfunction import DFunction
from ..core.valueaxis import ValueAxis
from ..core.frequency import FrequencyAxis

from .. import signal_TOTL
from .. import TWOD_SIGNALS
from .. import part_REAL, part_IMAGINARY, part_ABS

import quantarhei as qr

class TwoDSpectrum(DataSaveable, Saveable):
    
    dtypes = TWOD_SIGNALS
    
    def __init__(self):
        
        self.xaxis = None
        self.yaxis = None
        self.data = None
        self.dtype = None
        
        self.t2 = -1.0
        
        self.params = None
        
        
    def set_axis_1(self, axis):
        """Sets the x-axis of te spectrum (omega_1 axis)
        
        """
        self.xaxis = axis


    def set_axis_3(self, axis):
        """Sets the y-axis of te spectrum (omega_3 axis)
        
        """
        self.yaxis = axis


    def set_data_type(self, dtype=signal_TOTL): #"Tot"):
        """Sets the data type for this 2D spectrum
        
        Parameters
        ----------
        
        dtype : string
           Specifies the type of data stored in this TwoDSpectrum object 
        """
        
        if dtype in self.dtypes.values():
            self.dtype = dtype
        else:
            raise Exception("Unknown data type for TwoDSpectrum object")


    def get_spectrum_type(self):
        
        return self.dtype


    def set_data(self, data, dtype=signal_TOTL): #"Tot"):
        """Sets the data of the 2D spectrum
        
        Sets the object data depending on the specified type and stores the
        type
        
        
        Parameters
        ----------
        
        data : 2D array
            Data of the spectrum, float or complex
            
        dtype : string
            Type of the data stored. Three values are allowed: if quantarhei
            is import as `qr` that they are qr.signal_REPH
            for rephasing spectra, qr.signal_NONR for non-rephasing spectra,
            and qr.signal_TOTL for total spectrum, which is the sum of both
            
        """
        self.set_data_type(dtype)
        
        if (self.xaxis is None) or (self.yaxis is None):
            raise Exception("Axes of the 2D spectrum are not set")
            
        if ((self.xaxis.length == data.shape[0]) and
            (self.yaxis.length == data.shape[1])):
            
            self.data = data
            
        else:
            raise Exception("Axes not compatible with data\n"+
                            "xaxis.length = "+str(self.xaxis.length)+"\n"+
                            "yaxis.length = "+str(self.yaxis.length)+"\n"+
                            "data shape = "+str(data.shape))
            

    def add_data(self, data): 
        """Sets the data of the 2D spectrum

        
        """
        if self.data is None:
            raise Exception("Data is not initialized: use set_data method.")
        self.data += data 

       
    def set_t2(self, t2):
        """Sets the t2 (waiting time) of the spectrum
        
        
        """
        self.t2 = t2


    def get_t2(self):
        """Returns the t2 (waiting time) of the spectrum
        
        """
        return self.t2


    def log_params(self, params=None):
        """Store custom information about the spectrum
        
        
        This can be used e.g. for storage of parameters of the system for
        which the spectrum was calculated
        
        """
        self.params = params
        

    def get_log_params(self):
        return self.params
    
    
    
    def get_value_at(self, x, y):
        """Returns value of the spectrum at a given coordinate
        
        """

        if self.dtype is None:
            raise Exception("Data type not set")

        (ix, dist) = self.xaxis.locate(x)
        (iy, dist) = self.yaxis.locate(y)    

        return self.data[iy,ix]


    def get_cut_along_x(self, y0):
        """Returns a DFunction with the cut of the spectrum along the x axis
        
        """
        (iy, dist) = self.yaxis.locate(y0)
        
        ax = self.xaxis
        vals = numpy.zeros(ax.length, dtype=self.data.dtype)
        for ii in range(ax.length):
            vals[ii] = self.data[iy, ii]
    
        return DFunction(ax, vals)
    

    def get_cut_along_y(self, x0):
        """Returns a DFunction with the cut of the spectrum along the y axis
        
        """
        (ix, dist) = self.xaxis.locate(x0)
        
        ay = self.yaxis
        vals = numpy.zeros(ay.length, dtype=self.data.dtype)
        for ii in range(ay.length):
            vals[ii] = self.data[ii, ix]
    
        return DFunction(ay, vals)   
    

    def get_cut_along_line(self, point1, point2, which_step=None, step=None):
        """Returns a cut along a line specified by two points
        
        """
        vx1 = point1[0]
        vy1 = point1[1]
        vx2 = point2[0]
        vy2 = point2[0]
        
        (x1, dist) = self.xaxis.locate(vx1)
        (y1, dist) = self.yaxis.locate(vy1)
        (x2, dist) = self.xaxis.locate(vx2)
        (y2, dist) = self.yaxis.locate(vy2)
        
        length = numpy.sqrt((vx1-vx2)**2 + (vy1-vy2)**2)
        
        if which_step is None:
            wstep = "x"
        else:
            wstep = which_step
            
        if step is None:
            if wstep == "x":
                dx = self.xaxis.step
            elif wstep == "y":
                dx = self.yaxis.step
            else:
                raise Exception("Step along the cut line is not defined")
                
        else:
            dx = step
             
        Nstep = int(length/dx)
        
        axis = ValueAxis(0.0, Nstep+1, dx)
            
        vals = numpy.zeros(Nstep+1, dtype=self.data.dtype)
        ii = 0
        for val in axis.data:
            vx1 = self.xaxis.data[x1]
            vx2 = self.xaxis.data[x2]
            vy1 = self.yaxis.data[y1]
            vy2 = self.yaxis.data[y2]
            x = vx1 + val*(vx2-vx1)/(Nstep*dx)
            y = vy1 + val*(vy2-vy1)/(Nstep*dx)
            vals[ii] = self.get_value_at(x,y)
            ii += 1
            
        return DFunction(axis, vals)    
        
    
    def get_diagonal_cut(self):
        """Returns cut of the spectrum along the diagonal 
        
        
        """
        
        point1 = [self.xaxis.min, self.yaxis.min]
        point2 = [self.xaxis.max, self.yaxis.max]
        
        fce = self.get_cut_along_line(point1, point2, which_step="x")
        
        fce.axis.data += point1[0]
        
        return fce
    

    def get_anti_diagonal_cut(self, point):
        
        pass


    def get_max_value(self, dpart=part_REAL):
        """Maximum value of the real part of the spectrum
        
        
        """
        if dpart == part_REAL:
            return numpy.amax(numpy.real(self.data))
        elif dpart == part_IMAGINARY:
            return numpy.amax(numpy.imag(self.data))
        elif dpart == part_ABS:
            return numpy.amax(numpy.abs(self.data))
        else:
            raise Exception("Unknown data part")


    def get_min_value(self, dpart=part_REAL):
        """Minimum value of the real part of the spectrum
        
        
        """
        if dpart == part_REAL:
            return numpy.amin(numpy.real(self.data))
        elif dpart == part_IMAGINARY:
            return numpy.min(numpy.imag(self.data))
        elif dpart == part_ABS:
            return numpy.amin(numpy.abs(self.data))
        else:
            raise Exception("Unknown data part")


    def get_area_integral(self, area, dpart=part_REAL):
        """Returns an integral of a given area in the 2D spectrum
        
        """
        def integral_square(x1, x2, y1, y2, data, dx, dy):
            (n1, n2) = data.shape
            data.reshape(n1*n2)
            return numpy.sum(data)*dy*dy
        
        area_shape = area[0]
        x1 = area[1][0]
        x2 = area[1][1]
        y1 = area[1][2]
        y2 = area[1][3]
        
        dx = self.xaxis.step
        dy = self.yaxis.step
        
        (nx1, derr) = self.xaxis.locate(x1)
        (nx2, derr) = self.xaxis.locate(x2)
        (ny1, derr) = self.yaxis.locate(y1)
        (ny2, derr) = self.yaxis.locate(y2)
        
        if area_shape == "square":
            int_fce = integral_square
        else:
            raise Exception("Unknown area type: "+area_shape)   
            
        data = self.data[nx1:nx2, ny1:ny2]
        
        if dpart == part_REAL:
            return int_fce(x1, x2, y1, y2, numpy.real(data), dx, dy)
        elif dpart == part_IMAGINARY:
            return int_fce(x1, x2, y1, y2, numpy.imag(data), dx, dy)
        elif dpart == part_ABS:
            return int_fce(x1, x2, y1, y2, numpy.abs(data))
        else:
            raise Exception("Unknown data part")


    def get_area_max(self, area, dpart=part_REAL):
        """Returns a max value in a given area in the 2D spectrum
        
        """
        def find_in_square(x1, x2, y1, y2, data, dx, dy):
            return numpy.amax(data)
        
        area_shape = area[0]
        x1 = area[1][0]
        x2 = area[1][1]
        y1 = area[1][2]
        y2 = area[1][3]
        
        dx = self.xaxis.step
        dy = self.yaxis.step
        
        (nx1, derr) = self.xaxis.locate(x1)
        (nx2, derr) = self.xaxis.locate(x2)
        (ny1, derr) = self.yaxis.locate(y1)
        (ny2, derr) = self.yaxis.locate(y2)
        
        if area_shape == "square":
            int_fce = find_in_square
        else:
            raise Exception("Unknown area type: "+area_shape)   
            
        data = self.data[ny1:ny2, nx1:nx2]
        
        if dpart == part_REAL:
            return int_fce(x1, x2, y1, y2, numpy.real(data), dx, dy)
        elif dpart == part_IMAGINARY:
            return int_fce(x1, x2, y1, y2, numpy.imag(data), dx, dy)
        elif dpart == part_ABS:
            return int_fce(x1, x2, y1, y2, numpy.abs(data), dx, dy)
        else:
            raise Exception("Unknown data part")


    def normalize2(self, norm=1.0, dpart=part_REAL, nmax=None, use_max=False):
        """Normalizes the spectrum to the given maximum.
        
        Normalizes the spectrum to the maximum of a given part of the spectrum.
        Parts are real, imaginary or absolute values. Normalization is to 
        a supplied norm. If `nmax` is not specified, `use_max` is not
        taken into account. If `nmax` is a list-like object, the maximum
        of the current spectrum is appended to it. If `use_max` is set to False
        we use the first element of the list-like object `nmax` as the
        maximum of the present data.
        
        
        Parameters
        ----------
        
        norm : float
            Norm to which we normalize. Default value is 1.
            
        dpart : string
            Part of the spectrum we normalize against.
            
        namx : list-like or None
            If `nmax` is a list-like object, the maximum of the current 
            spectrum is appended to it.
            
        use_max : bool
            If `use_max` is set to False we use the first element of the 
            list-like object `nmax` as the  maximum of the present data.
        
        
        """           
        mx = self.get_max_value(dpart=dpart)
        if nmax is not None:
            nmax.append(mx)
            if not use_max:
                mx = nmax[0]
        if mx != 0.0:
            self.data = (self.data/mx)*norm
        
    
    def devide_by(self, val):
        """Devides the total spectrum by a value
        
        
        Parameters
        ----------
        
        val : float, int
            Value by which we devide the spectrum
        
        """
        
        self.data = self.data/val

    
    def _interpolate(self):
        """Interpolate the spectrum by splines
        
        
        """
        pass


    def zeropad(self, fac=1):
        if fac == 1:
            return

        Nxa = self.xaxis.length
        Nya = self.yaxis.length
        
        print("Start: ", self.xaxis.start)
        print("Step:  ", self.xaxis.step)
        print("end:   ", self.xaxis.step*self.xaxis.length + self.xaxis.start)
        
        nNxa = fac*Nxa
        nNya = fac*Nya
        
        print("New length x: ", nNxa)
        print("New length y: ", nNya)
        

    def shift_energy(self, dE, interpolation="linear"):
        """Shift the spectrum in both frequency axis by certain amount
        
        The shift is specified as the energy that has to be added to the
        2D spectrum to shifted to its expected position. As a result
        we have a new spectrum 
        
        S(E_3, E_1) = S(E_3-dE, E_1-dE)
        
        Because of the negative sign above, we have to change the sign of dE
        as a first thing in the code. This is for user convenience.

        """

        dE = -dE
        
        ndata = numpy.zeros(self.data.shape, dtype=self.data.dtype)
        ndata1 = numpy.zeros(self.data.shape, dtype=self.data.dtype)
        
        if interpolation == "linear":
            
            print("Shifting spectrum: linear interpolation")
            
            data = self.data
    
            N1 = self.data.shape[0]
            N2 = self.data.shape[1]
            
            Dx = self.xaxis.step
            
            A = numpy.abs(dE)
            n = int(numpy.floor(A/Dx))
            dx = A - n*Dx   # dx is positive
            
            s = int(numpy.sign(dE))
            
            # dimesion 1
            # make sure we stay within the defined array
            if s == 1:
                n1 = 0
                n2 = N1 - (n+1)
            else:
                n1 = n+1
                n2 = N1
                
            for i1 in range(n1, n2):
                
                ndata1[i1, :] = (data[i1+s*n,:]
                + (data[i1+s*(n+1),:] - data[i1+s*n,:])*(dx/Dx))
            
            # dimension 2
            # make sure we stay within the defined array
            if s == 1:
                n1 = 0
                n2 = N2 - (n+1)
            else:
                n1 = n+1
                n2 = N2
    
            for i2 in range(n1, n2):
                
                ndata[:, i2] = (ndata1[:, i2+s*n]
                + (ndata1[:, i2+s*(n+1)] - ndata1[:,i2+s*n])*(dx/Dx))
                
        elif interpolation == "spline":
            
            pass
        
        elif interpolation == "fft":
            
            print("Shifting spectrum: fft interpolation")
            
            # inverse FFT
            ndata = numpy.fft.fftshift(self.data)
            ndata1 = numpy.fft.fft2(ndata)
            ndata = numpy.fft.fftshift(ndata1)

            timex = numpy.fft.fftfreq(self.xaxis.length,
                                                     self.xaxis.step)
            timex = numpy.fft.fftshift(timex)
            timey = numpy.fft.fftfreq(self.yaxis.length,
                                                     self.yaxis.step)
            timey = numpy.fft.fftshift(timey)
            
            # multiply by exponentials
            etx = numpy.exp(1j*dE*timex)
            ety = numpy.exp(1j*dE*timey)
            
            for k in range(ndata.shape[0]):
                ndata1[k,:] = etx[k]*ndata[k,:]*ety[:]
            
            # back FFT
            ndata = numpy.fft.ifft2(ndata1)
            ndata = numpy.fft.fftshift(ndata)
            
        else:
            raise Exception("Unknown interpolation type")
            
        self.data[:,:] = ndata[:,:]
            

    # FIXME: implement this
    def get_PumpProbeSpectrum(self):
        """Returns a PumpProbeSpectrum corresponding to the 2D spectrum
        
        """
        #from .pumpprobe import PumpProbeSpectrumCalculator
        from . import pumpprobe as pp
        #from ..core.time import TimeAxis
        #fake_t = TimeAxis(0,1,1.0)
        #ppc = PumpProbeSpectrumCalculator(fake_t, fake_t, fake_t)
        #return ppc.calculate_from_2D(self)
        return pp.calculate_from_2D(self)


    def plot(self, fig=None, window=None, 
             stype=None, spart=part_REAL,
             vmax=None, vmin_ratio=0.5, 
             colorbar=True, colorbar_loc="right",
             cmap=None, Npos_contours=10,
             show_states=None,
             text_loc=[0.05,0.9], fontsize="20", label=None,
             show=False,
             show_diagonal=None,
             xlabel=None,
             ylabel=None,
             axis_label_font=None):
        """Plots the 2D spectrum
        
        Parameters
        ----------
        
        fig : matplotlib.figure
            Figure into which plotting will be done. This is used e.g. when
            making a movie using moview writter (may be obsolete).
            If fig is None, we create a new figure object.
            
        window : list
            Specifies the plotted window in current energy units. When axes
            are x and y, the window is specified as window=[x_min,x_max,
            y_min,y_max]
            
        spart : {part_REAL, part_IMAGINARY, part_ABS}
            part of the spectrum to be plotted, constants such as part_REAL are
            defined at the highest import level of quantarhei.
            
        vmax : float
            max of the plotting range in the z-direction. If vmax is None,
            maximum of the real part of the spectrum is used to determine
            the values of `vmax`
            
            
            
            
        """
        spect2D = self.data
        
        #
        # What part of the spectrum to plot
        #
        if spart == part_REAL:
            spect2D = numpy.real(spect2D)
        elif spart == part_IMAGINARY:
            spect2D = numpy.imag(spect2D)
        elif spart == part_ABS:
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
    
        #
        #  How to treat the figures
        #
        if fig is None:
            fig, ax = plt.subplots(1,1)
        else:
            fig.clear()
            fig.add_subplot(1,1,1)
            ax = fig.axes[0]
            
        #
        # Color map
        #
        if cmap is None:
            cmap = plt.cm.rainbow
            
            
        #
        # Actual plotting
        #
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

        #
        # Label
        #
        if label is not None:
            if isinstance(label, str) and len(text_loc) == 2:
                text_loc = [text_loc]
                label = [label]
                ln_t = 1
                ln_l = 1
            else:
                try:
                    ln_t = len(text_loc)
                    ln_l = len(label)
                except:
                    raise Exception("text_loc and label parameters must be"+
                                    " lists of the same lengths")
    
            if ln_t != ln_l:
                raise Exception("text_loc and label parameters have to have the"
                                +" same number of members")
                
            kk = 0
            for pos in text_loc:
                
                lbl = label[kk]
                if lbl is not None:
                    ax.text((prvo-levo)*pos[0]+levo,
                        (hore-dole)*pos[1]+dole,
                        lbl,
                        fontsize=str(fontsize))
                kk += 1
                
        
        #
        # Contours
        #
        
        # positive contours are always plotted
        warnings.filterwarnings("error")
        
        try:
            plt.contour(self.xaxis.data[i1_min:i1_max],
                     self.yaxis.data[i3_min:i3_max],
                     realout, levels=poslevels, colors="k")
                     #linewidth=1)
        except:
            pass
            #print("No positive contours found; not plotted")
              
        # other contours only if we do not plot absolute values
        if spart != "abs":
            # zero contour
            try:
                plt.contour(self.xaxis.data[i1_min:i1_max],
                         self.yaxis.data[i3_min:i3_max],
                         realout, levels=[0],colors="b")
                         #linewidth=1)
            except:
                pass
                #print("Zero contour not found; not plotting")
        
        
            # negative contours
            try:
                plt.contour(self.xaxis.data[i1_min:i1_max],
                         self.yaxis.data[i3_min:i3_max],
                         realout, levels=neglevels, colors="k")
                         #linewidth=1) 
            except:
                pass
                #print("Negative contour not found; not plotting")
        
        warnings.resetwarnings()
        
        #
        # Color bar presence
        #
        if colorbar:
            plt.clim(vmin=vmin,vmax=vmax)
            fig.colorbar(cm)
            
        #
        # Plot lines denoting positions of selected transitions
        #
        if show_states is not None:
            for en in show_states:
                try:
                    en1 = en[0]
                    co1 = en[1]
                except:
                    en1 = en
                    co1 = '--k'
                    
                if en1 >= levo and en1 <= prvo:
                    plt.plot([en1,en1],[dole,hore],co1,linewidth=1.0)
                if en1 >= dole and en1 <= hore:
                    plt.plot([levo,prvo],[en1,en1],co1,linewidth=1.0)
          
        #
        # show diagonal line
        #
        if show_diagonal is not None:
            plt.plot([levo,prvo],[dole,hore], '-k',linewidth=1.0)

        #
        # axis labels
        #
        if axis_label_font is not None:
            font = axis_label_font
        else:
            font={'size':20}

        if xlabel is None:
            xl = ""
        if ylabel is None:
            yl = ""
            
        if xlabel is not None:
            xl = r'$\omega$ [fs$^{-1}$]'

        if isinstance(self.xaxis, FrequencyAxis):
            units = self.xaxis.unit_repr_latex()
            xl = r'$\omega_{1}$ ['+units+']'
            yl = r'$\omega_{3}$ ['+units+']'
#        if isinstance(self.axis, TimeAxis):
#            xl = r'$t$ [fs]'
#            yl = r'$f(t)$'

        if xlabel is not None:
            xl = xlabel
        if ylabel is not None:
            yl = ylabel

        if xl is not None:
            plt.xlabel(xl, **font)
        if xl is not None:
            plt.ylabel(yl, **font)            
            
            
        #
        # Should the spectra be showed now?
        #
        if show:
            self.show()

            
    def show(self):
        """Show the plot of 2D spectrum
        
        By default, plots are not shown. It is waited until explicit show()
        is called
        
        """
        plt.show()
        

    def savefig(self, filename):
        """Saves the fige of the plot into a file
        
        """
        plt.savefig(filename)
            

    def trim_to(self, window=None):
        """Trims the 2D spectrum to a specified region
        
        Parameters
        ----------
        
        window : array
            Spectral window to which the present 2D spectrum 
            should be trimmed. The window is specified as
            an array [w1_min, w1_max, w3_min, w3_max], where 
            w1_min is the lower bound of the w1 axis (the x-axis),
            w1_max is the upper bound of the w1 axis and similarly
            for the w3 axis (y-axis) of the spectrum.
            
        
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
            atype = self.xaxis.atype
                
            xaxis = FrequencyAxis(start_1,length_1,step_1, 
                                  atype=atype,
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
            if self.data is not None:
                data = self.data[i1_min:i1_max,i3_min:i3_max]
                self.data = data
                        
        else:
            # some automatic trimming in the future
            pass
         