# -*- coding: utf-8 -*-
#mport ..utils as utils
from ..core.valueaxis import ValueAxis
from ..core.units import cm2int

import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from ..utils import Integer, Float
from ..utils.vectors import X


class twodspect:
    """ Two-Dimensional Fourier Transformed Optical Spectrum 
    
    
    
    """
    
    def __init__(self,xaxis=ValueAxis(0.0,10,0.1),yaxis=ValueAxis(0.0,10,0.1)):
        
        if (isinstance(xaxis,ValueAxis) and isinstance(yaxis,ValueAxis)):
            self.xaxis = xaxis
            self.yaxis = yaxis
        else:
            raise TypeError("Expecting ValueAxis")
        
        self.renderable_components = []
        
        
    def add_component(self,lineshape):
        
        # CHECK if lineshape is a RENDERABLE
        if isinstance(lineshape,renderable):
            self.renderable_components.append(lineshape)
            lineshape.set_twospect(self)
        else:
            raise TypeError("Argument must be renderable")
        
    
    def render_data(self):
        data = None
        
        for rc in self.renderable_components:
            if data is not None:
                data = data + rc.render_data()
            else:
                data = rc.render_data()
                
        return data

    def plot_data(self,data,part="real",vrange=(0.0,1.0),
                  normalize=True,ncontours=20,title=""):
        """Plots 2D spectrum 
        
        
        
        
        
        """ 
        
        # number of contours will be used to display the data
        Nstep = ncontours
        if normalize:
            # maximum value of the data can be used normalize data
            dmax = numpy.max(numpy.real(data))
            data = data/dmax
        if vrange is None:
            if normalize:
                urange = 1.0
            else:
                urange = numpy.max(numpy.real(data))
            lrange = numpy.min(numpy.real(data))
        else:
            lrange = vrange[0]
            urange = vrange[1]
                
                
        step = (urange-lrange)/Nstep
        levels = numpy.arange(lrange, urange, step)
        
        self.fig = plt.figure()


        if part == "imag":
            plt.contour(self.xaxis.data,self.yaxis.data,numpy.imag(data),
                        levels)
        else:
#            nx = len(self.xaxis.data)
#            ny = len(self.yaxis.data)
#            plt.imshow(numpy.real(data),
#                       interpolation='bilinear',
#                       origin='lower',
#                       cmap=cm.jet, extent=(self.xaxis.data[0], 
#                                              self.xaxis.data[nx-1],
#                                              self.yaxis.data[0],
#                                              self.yaxis.data[ny-1]))            
            plt.contourf(self.xaxis.data,self.yaxis.data,numpy.real(data),
                        levels,cmap=cm.jet)
            plt.colorbar()
            plt.contour(self.xaxis.data,self.yaxis.data,numpy.real(data),
                        levels,colors="k")
            plt.title(title)
        
        
        plt.show()        
    
    def savefig(self,filename):
        self.fig.savefig(filename)
        
        
    
class labsetup:
    """Laboratory set-up for non-linear spectroscopy
    
    
    
    """
    
    number_of_pulses = Integer("number_of_pulses")
    
    def __init__(self, nopulses = 3):
        
        self.number_of_pulses = nopulses
    
        self.M4 = numpy.array([[4.0, -1.0, -1.0],
                               [-1.0, 4.0, -1.0],
                               [-1.0,-1.0,  4.0]])/30.0
                        
    
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
    
                         
class renderable:

    def render_data(self):
        pass     
                   
                        
class lineshape_generator(renderable):

    def __init__(self):
        pass


    def render_data(self):
        return numpy.zeros((self.Nx,self.Ny))
        
    def set_twospect(self,tsp):
        self.twospect = tsp
        self.Nx = self.twospect.xaxis.noPoints
        self.Ny = self.twospect.yaxis.noPoints
        
    
    
class gaussian_lineshape_generator(lineshape_generator):

    xwidth = Float("xwidth")
    ywidth = Float("ywidth")
    xposition = Float("xposition")
    yposition = Float("yposition")
    prefactor = Float("prefactor")
    phase = Float("phase")
    correlation_time = Float("correlation_time")
    population_time = Float("population_time")
    corr = Float("corr")
    
    def __init__(self,xwidth=1.0,ywidth=1.0):
        """ 2D lineshape generator of Gaussian type
        
        The lineshape has a form of 
        
        f(x,y) = exp[-((x-xposition)/xwidth)**2 - ((y-yposition)/ywidth)**2
        + corr*(x-xposition)*(y-yposition)
        *exp(-population_time/correlation_time)]
        
    
        Parameters
        ----------
    
        xwidth
            Width of the Gaussian lineshape in the x-direction.
        
        ywidth 
            Width of the Gaussian lineshape in the y-direction.
        
        Properties
        ----------
        
        prefactor
            Prefactor which multiplies the lineshape which otherwise 
            has maximum values of 1 at its center.
            
        xposition
            Position of the lineshape maximum on x-axis.
            
        yposition
            Position of the lineshape maximum on y-axis.
            
        phase
            Phase of the lineshape. The real values standard lineshape
            gets multiplied by a factor exp(-i*phase).
            
        correlation_time
            This is the lifetime of x-y correlation of the lineshapes.
            
        population_time
            Population time in which the 2D spectrum is taken.
            
        corr
            This factor is set to 2.0 for a perfect correlation of x 
            and y axis. It is smaller than 2.0 in realistic cases and can be
            set to zero to forbid correlation altogether.
        
        
        """        
        self.xwidth = xwidth
        self.ywidth = ywidth
        self.corr = 0.0
        self.phase = 0.0
        self.correlation_time = 1.0e10       
        
        self.prefactor = 1.0
        self.xposition = 0.0
        self.yposition = 0.0
        self.population_time = 0.0
        self.liouville_pathways = None
        
        self.frequency_domain = True
        
        
    def set_liouville_pathways(self,lps):
        self.liouville_pathways = lps
        
    def append_liouville_pathways(self,lps):
        if self.liouville_pathways is None:
            self.liouville_pathways = lps
        else:
            for lp in lps:
                self.liouville_pathways.append(lp)
        
    def render_lineshape(self,ptype="R"):
        data = numpy.zeros((self.Nx,self.Ny))
        
        kx = 0

        y = self.twospect.yaxis.data 
        fy = numpy.exp(-((y-self._yposition)/self._ywidth)**2)
        
        cc = 1.0
                
        for x in self.twospect.xaxis.data:
                        
            data[:,kx] = self._prefactor*fy\
              *numpy.exp(-((x-self._xposition)/self._xwidth)**2)\
              *numpy.exp(numpy.exp(-self.population_time/self.correlation_time)
                        *cc*(x-self._xposition)*
                        (y-self._yposition)/(self._xwidth*self._ywidth))
              
            kx += 1

        return data*numpy.exp(-1j*self.t2_phase)*self.t2_amplitude
        
                
        
    def render_liouville_pathways(self,lps):
        
        pref_tol = 1.0e-4        
        
        # scan prefactors
        prefs = numpy.zeros(len(lps))
        k = 0
        for lp in lps:
            prefs[k] = abs(lp.pref)
            k += 1
            
        pmax = numpy.max(prefs)
        pref_tol = pmax*pref_tol
        
        data = numpy.zeros((self.Nx,self.Ny))
        k = 0
        l = 0
        for lp in lps:
            
            if abs(lp.pref) > pref_tol:
                
                self.prefactor = numpy.float(numpy.real(lp.pref))
            
                # Frequency of the first interval
                if lp.pathway_type == "R":
                    self.xposition = -lp.get_interval_frequency(0)
                elif lp.pathway_type == "NR":
                    self.xposition = lp.get_interval_frequency(0)
                else:
                    raise Exception("Unsupported pathway type")
                    
                # Frequency of the third interval
                self.yposition = lp.get_interval_frequency(2+lp.relax_order)
                
                #----------------------------------------------
                # Population interval is characterized only by 
                # phase evolution
                #----------------------------------------------
                
                # Frequency of the second (population) interval
                omega_2 = lp.get_interval_frequency(1)

                # Phase of the population contribution                
                self.t2_phase = omega_2*\
                             cm2int*self.population_time
                self.t2_amplitude = 1.0

                # render the corresponding 2D data                        
                data = data + self.render_lineshape(ptype=lp.pathway_type)
                
                l += 1
                
            k += 1
        
        print("Pathways eliminated: ", l, " out of ", k, " rendered")        
        
        return data
        
    
    def render_data(self):
        
        if self.liouville_pathways:
            return self.render_liouville_pathways(self.liouville_pathways)
        else:
            return self.render_lineshape()
            

class aggregate_lineshape_generator(gaussian_lineshape_generator):
    """ 
    
    
    """
    
    Uet = None
    Uex = None
    
    def __init__(self,aggregate,Rate=0.000001):
        self.aggregate = aggregate
        self.last_t2 = -1.0
        self.Rate = Rate
        super().__init__(xwidth=200.0,ywidth=200.0)
    
    def render_lineshape(self,ptype="R"):
        
        #return super().render_lineshape(ptype=ptype)
        
        data = numpy.zeros((self.Nx,self.Ny))

        # define time dependent grid depending on self.twospect.yaxis.data
        # and self.twospect.xaxis.data  

        #
        #  Calculate time-dependent signal as a function of t1, t2, t3
        #
        
        # Fourier transform it into a line shape
        
        # Coherence pathways only decay
        
        # Population pathways are calculated in a Markovian way
        
        
        cc = 1.0
        
        y = self.twospect.yaxis.data 
        fy = numpy.exp(-((y-self._yposition)/self._ywidth)**2)
        
        kx = 0        
        for x in self.twospect.xaxis.data:
                        
            data[:,kx] = self._prefactor*fy\
              *numpy.exp(-((x-self._xposition)/self._xwidth)**2)\
              *numpy.exp(numpy.exp(-self.population_time/self.correlation_time)
                        *cc*(x-self._xposition)*
                        (y-self._yposition)/(self._xwidth*self._ywidth))
              
            kx += 1


        return data*self.t2_ampphase     
        
        

    def render_liouville_pathways(self,lps):
        
        pref_tol = 1.0e-4        
        
        # scan prefactors
        prefs = numpy.zeros(len(lps))
        k = 0
        for lp in lps:
            prefs[k] = abs(lp.pref)
            k += 1
            
        pmax = numpy.max(prefs)
        pref_tol = pmax*pref_tol
        
        print(pmax,pref_tol)
        
        data = numpy.zeros((self.Nx,self.Ny))
        k = 0
        l = 0
        for lp in lps:
            
            if abs(lp.pref) > pref_tol:
                
                self.prefactor = numpy.float(numpy.real(lp.pref))
            
                self.set_Us(self.population_time,Rate=self.Rate)
                
                # Frequency of the first interval                
                if lp.pathway_type == "R":
                    self.xposition = -lp.get_interval_frequency(0)
                elif lp.pathway_type == "NR":
                    self.xposition = lp.get_interval_frequency(0)
                elif lp.pathway_type == "R_E":
                    self.xposition = -lp.get_interval_frequency(0)
                elif lp.pathway_type == "NR_E":
                    self.xposition = lp.get_interval_frequency(0)
                else:
                    raise Exception("Unsupported pathway type")
                    
                # Frequency of the third interval
                self.yposition = lp.get_interval_frequency(2+lp.relax_order)

                #----------------------------------------------
                # Population interval is characterized  
                # by a corresponding evolution super operator
                # element, which comprises both the phase evolution
                # and an amplitude evolution
                #----------------------------------------------
                
                if lp.pathway_type in ["R","NR"]:
                    
                    if lp.popt_band == 1:
                        # population or electronic coherence propagation
                        # is after the second event
                        ln = lp.states[1,0]
                        rn = lp.states[1,1]                             
                             
                        self.t2_ampphase = self.Uex[ln-self.aggregate.Nb[0],
                                                rn-self.aggregate.Nb[0]]
                        
                    else:                        
                        # Frequency of the second (population) interval
                        omega_2 = lp.get_interval_frequency(1)
                        self.t2_ampphase = numpy.exp(-1j*omega_2*
                             cm2int*self.population_time)

                elif lp.pathway_type in ["R_E", "NR_E"]:

                    sf = lp.relaxations[0][0]
                    si = lp.relaxations[0][1]
                    
                    ag = sf[0]
                    bg = sf[1]
                    ae = si[0]
                    be = si[1]
                                        
                    self.t2_ampphase = self.Uet[ag,bg,ae-self.aggregate.Nb[0],
                                                      be-self.aggregate.Nb[0]]                            
                            

                # render the corresponding 2D data                        
                data = data + self.render_lineshape(ptype=lp.pathway_type)
                
                l += 1
                
            k += 1
        
        print("Pathways eliminated: ", l, " out of ", k, " rendered")        
        
        return data

    def set_Us(self,t,Rate=0.0):
        
        if (self.Uet is None) or (self.last_t2 != t):
            #print("Setting Us")
            self.Uet, self.Uex = \
            self.aggregate.ETICS_evolution_operator(
                                  self.population_time,Rate,[1])
            self.last_t2 = t     
            
    class evolution_lineshape_generator():
        pass
        
        
       
