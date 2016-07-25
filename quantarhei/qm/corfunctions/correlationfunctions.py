# -*- coding: utf-8 -*-

import numpy
import scipy.interpolate

from ...core.dfunction import DFunction
from ...core.units import kB_intK

from ...core.managers import UnitsManaged

class CorrelationFunction(DFunction,UnitsManaged):
    """Provides typical Bath correlation functions
        
    Types of correlation function are provided:
    -----------------------------
    OverdampedBrownian-HighTemperature : 
        OverdampedBrownian oscillator in high temperature limit
            
    OverdampedBrownian :
        General overdampedBrownian oscillator
       
    """       

    
    allowed_types = ("OverdampedBrownian-HighTemperature",\
                     "OverdampedBrownian","Value-defined")    
    
    def __init__(self,timeAxis,params,values=None,domain='Time'):
        self.valueAxis = timeAxis
        self.timeAxis = self.valueAxis
        self.params = params
        if domain == 'Time':
            self.time_domain = True
        else:
            self.time_domain = False
        
        self.is_composed = False
        
        
        try:
            ftype = params["ftype"]
            if ftype in CorrelationFunction.allowed_types:
                self.ftype = ftype
            else:
                raise Exception("Unknown Correlation Function Type")
                 
            # we need to save the defining energy units
            self.energy_units = self.manager.get_current_units("energy")
                 
             
        except:
            raise Exception
            
        if ((self.ftype == "OverdampedBrownian-HighTemperature") 
             and self.time_domain):

            
            T    = params["T"]
            tc   = params["cortime"]
            # use the units in which params was defined
            lamb = self.manager.iu_energy(params["reorg"],
                                          units=self.energy_units)    
            
            
            kBT = kB_intK*T
            
            t = self.timeAxis.time
            
            cc = 2.0*lamb*kBT*numpy.exp(-t/tc) \
                  + 1.0j*(lamb/tc)*numpy.exp(-t/tc)
            
            self.data = cc
            self.lamb = lamb
            self.T = T
            
            self.cutoff_time = 5.0*tc
            
        elif ((self.ftype == "OverdampedBrownian") and self.time_domain):
            
            T    = params["T"]
            tc   = params["cortime"]
            lamb = self.manager.iu_energy(params["reorg"],
                                          units=self.energy_units)    
            
            try:
                N = params["matsubara"]
            except:
                N = 10
            
            kBT = kB_intK*T            
            
            t = self.timeAxis.time           
           
            cc = (lamb/(tc*numpy.tan(1.0/(2.0*kBT*tc))))*numpy.exp(-t/tc) \
                  + 1.0j*(lamb/tc)*numpy.exp(-t/tc)
                  
            cc += \
            (4.0*lamb*kBT/tc)*CorrelationFunction.__matsubara(kBT,tc,N,t)
                  
            self.data = cc
            self.lamb = lamb
            self.T    = T
                 
            self.cutoff_time = 5.0*tc
            
        # FIXME: here I should be able to set a frequency domain spectral dens
        elif ((self.ftype == "Value-defined") and self.time_domain):
            
            lamb = self.manager.iu_energy(params["reorg"],
                                          units=self.energy_units)    
            self.lamb = lamb
            try:
                tcut = params["cutoff-time"]
            except:
                tcut = self.timeAxis.tmax
            self.cutoff_time = tcut
            
            if len(values) == self.timeAxis.length:
                self.data = values
            else:
                raise Exception("Incompatible values")
            
        else:
            raise Exception("Unknown correlation function type of"+
                            "type domain combination.")
        
    def __matsubara(kBT,tc,N,t):
        ms = 0.0
        
        nu = 2.0*numpy.pi*kBT
        for ii in range(0,N):
            nn = ii+1
            ms += nu*nn*numpy.exp(-nu*nn*t)/((nu*nn)**2-(1.0/tc)**2)
        
        return ms
        
    def get_temperature(self):
        """Returns the temperature of the correlation function
        
        """
        return self.T
        
        
    def add(self,cf):
        """Adds another correlation function to the current one"""
        if ((self.timeAxis == cf.self.timeAxis) and (self.T == cf.T)):
            self.lamb += cf.lamb
            self.data += cf.data
            if self.cutoff_time < cf.cutoff_time:
                self.cutoff_time = cf.cutoff_time
            self.is_composed = True
            
        else:
            raise Exception()
        
    def copy(self):
        """Creates a copy of the current correlation function"""
        return CorrelationFunction(self.timeAxis,self.params)
        
    def convert_2_spectral_density(self):
        """Converts the data from time domain to frequency domain 
        
        
        """
        # if the correlation function is in the time domain, we can convert it
        if self.time_domain:
            
            tt = self.timeAxis.time
            ff = self.data
            
            # use the symetry of the correlation function
            for i in range(len(tt)//2):
                j = len(tt)-i-1
                ff[j] = numpy.conj(self.data[i+1])

            # frequencies            
            om = (2.0*numpy.pi)*numpy.fft.fftfreq(self.timeAxis.length,
                                   self.timeAxis.dt)
            # values
            fo = numpy.fft.fft(ff)

            # rotate the values and frequencies to have them in 
            # increasing order
            om = numpy.fft.fftshift(om)
            fo = numpy.fft.fftshift(fo)*self.timeAxis.dt

            # interpolate the spectral density
            sr = scipy.interpolate.UnivariateSpline(om,
                        numpy.real(fo),s=0)

            # save the data
            self.interp_data = sr
            self.data = fo
            self.time_domain = False
            

#            om = om/cm2int
#            fo = fo
#
#            f,(ax1,ax2) = plt.subplots(2)
#            ax1.plot(tt,numpy.real(ff),'-k')
#            ax1.plot(tt,numpy.imag(ff),'-g')
#            ax1.set_xlim([0.0,300.0])
#            ax2.plot(om,numpy.real(fo),'-k')                        
#            ax2.plot(om,numpy.imag(fo),'-g')
#            ax2.set_xlim([-300,300])
#                        
#            fr = numpy.zeros(400)
#            fr2 = numpy.zeros(400)
#            points = numpy.zeros(400)
#            for i in range(400):
#                j = i - 200
#                points[i] = j 
#                if (j <= 0):
#                    fr[i] = sr(j*cm2int)
#                    fr2[i] = sr(j*cm2int)
#                else:
#                    pref = numpy.exp(numpy.real(j*cm2int)/(kB_intK*self.T))
#                    print(i,j,pref,kB_intK*300.0/cm2int)
#                    #pref = 1.0
#                    fr[i] = pref*sr(-j*cm2int)
#                    fr2[i] = sr(-j*cm2int)
#            #fr2 = sr(points)
#            ax2.plot(points,fr,'-r')
#            ax2.plot(points,fr2,'-m')
#            
#            plt.show()
#            
        else:
            # nothing happens if the function is already in frequency domain
            pass
        
        
        
def c2g(timeaxis,coft):
    """ Converts correlation function to lineshape function
        
    Explicit numerical double integration of the correlation
    function to form a lineshape function.

    Parameters
    ----------

    timeaxis : cu.oqs.time.TimeAxis
        TimeAxis of the correlation function
            
    coft : complex numpy array
        Values of correlation function given at points specified
        in the TimeAxis object
            
        
    """
        
    ta = timeaxis
    rr = numpy.real(coft)
    ri = numpy.imag(coft)
    sr = scipy.interpolate.UnivariateSpline(ta.time,
                        rr,s=0).antiderivative()(ta.time)
    sr = scipy.interpolate.UnivariateSpline(ta.time,
                        sr,s=0).antiderivative()(ta.time)
    si = scipy.interpolate.UnivariateSpline(ta.time,
                        ri,s=0).antiderivative()(ta.time)
    si = scipy.interpolate.UnivariateSpline(ta.time,
                        si,s=0).antiderivative()(ta.time)
    gt = sr + 1j*si
    return gt
    
    
def c2h(timeaxis,coft):

    ta = timeaxis
    rr = numpy.real(coft)
    ri = numpy.imag(coft)
    sr = scipy.interpolate.UnivariateSpline(ta.time,
                        rr,s=0).antiderivative()(ta.time)
    si = scipy.interpolate.UnivariateSpline(ta.time,
                        ri,s=0).antiderivative()(ta.time)
    ht = sr + 1j*si
    
    return ht

def h2g(timeaxis,coft):
    return c2h(timeaxis,coft)

        