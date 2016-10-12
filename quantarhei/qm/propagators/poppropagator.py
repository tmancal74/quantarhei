# -*- coding: utf-8 -*-
import numpy


class PopulationPropagator:
    
    def __init__(self, timeaxis, RateMatrix=None):
        """ Propagator for a population vector 
        
        """
        self.timeAxis = timeaxis
        self.Nref = 1
        self.Nt = self.timeAxis.length
        self.dt = self.timeAxis.step
        
        if RateMatrix is not None:
            self.KK = RateMatrix
    
    def propagate(self, pini):
        if not isinstance(pini, numpy.ndarray):
            pini = numpy.array(pini)
            
        return self._propagate_short_exp(pini)    
        
        
    def _propagate_short_exp(self,pini,L=4):
        Nt = self.timeAxis.length
        pops = numpy.zeros((Nt,pini.shape[0]))

        rho1 = pini
        rho2 = pini
        
        indx = 1
        for ii in self.timeAxis.data[1:self.Nt]:
                                   
            for jj in range(0, self.Nref):
                
                for ll in range(1,L+1):
                    
                    pref = (self.dt/ll) 
                    rho1 = pref*numpy.dot(self.KK.data,rho1)
                             
                    rho2 = rho2 + rho1
                rho1 = rho2    
                
            pops[indx,:] = rho2                        
            indx += 1             
            
        return pops
        