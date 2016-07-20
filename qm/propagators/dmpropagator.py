# -*- coding: utf-8 -*-

import numpy

from .propagations import DMPropagation

"""


    DENSITY MATRIX PROPAGATOR


"""    
class DMPropagator:
    
    def __init__(self,timeaxis,ham):
        self.timeaxis = timeaxis
        self.ham = ham
        
        
        
        self.Odt = self.timeaxis.time[1]-self.timeaxis.time[0]
        self.dt = self.Odt
        self.Nref = 1
        
        self.Nt = self.timeaxis.time.shape[0]
        
        N = self.ham.data.shape[0]
        self.N = N
        self.data = numpy.zeros((self.Nt,N,N),dtype=numpy.complex64)        


    def propagate(self,rhoi):
        
        return self.__propagate_short_exp(rhoi,L=4)
        
        
    def __propagate_short_exp(self,rhoi,L=4):
        """
              Short exp integration
        """
        
        pr = DMPropagation(self.timeaxis,rhoi)
        rho1 = rhoi.data
        rho2 = rhoi.data
        
        HH = self.ham.data        
        
        indx = 1
        for ii in self.timeaxis.time[1:self.Nt]:
            
                    
            for jj in range(0,self.Nref):
                
                for ll in range(1,L+1):
                    
                    pref = (self.dt/ll) 
                    
                    rho1 = -1j*pref*(numpy.dot(HH,rho1) \
                             - numpy.dot(rho1,HH) )
                             
                    rho2 = rho2 + rho1
                rho1 = rho2    
                
            pr.data[indx,:,:] = rho2                        
            indx += 1             
            
        return pr