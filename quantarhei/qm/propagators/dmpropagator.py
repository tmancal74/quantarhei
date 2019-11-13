# -*- coding: utf-8 -*-
"""


    DENSITY MATRIX PROPAGATOR


""" 

import numpy

from .dmevolution import DensityMatrixEvolution

   
class DMPropagator:
    
    def __init__(self, timeaxis, ham):
        self.timeaxis = timeaxis
        self.ham = ham
        
        self.Odt = self.timeaxis.data[1]-self.timeaxis.data[0]
        self.dt = self.Odt
        self.Nref = 1
        
        self.Nt = self.timeaxis.length
        
        N = self.ham.data.shape[0]
        self.N = N
        self.data = numpy.zeros((self.Nt,N,N),dtype=numpy.complex64)        


    def propagate(self,rhoi):
        
        return self._propagate_short_exp(rhoi,L=4)
        
        
    def _propagate_short_exp(self,rhoi,L=4):
        """
              Short exp integration
        """
        
        pr = DensityMatrixEvolution(self.timeaxis,rhoi)
        rho1 = rhoi.data
        rho2 = rhoi.data
        
        HH = self.ham.data        
        
        indx = 1
        for ii in self.timeaxis.data[1:self.Nt]:
            
                    
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