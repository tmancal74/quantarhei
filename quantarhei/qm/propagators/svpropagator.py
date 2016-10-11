# -*- coding: utf-8 -*-
"""


    StateVector propagator


""" 

import numpy

from .statevectorevolution import StateVectorEvolution
from ..hilbertspace.evolutionoperator import EvolutionOperator

   
class StateVectorPropagator:
    
    def __init__(self, timeaxis, ham):
        self.timeaxis = timeaxis
        self.ham = ham
        
        self.Odt = self.timeaxis.data[1]-self.timeaxis.data[0]
        self.dt = self.Odt
        self.Nref = 1
        
        self.Nt = self.timeaxis.length
        
        N = self.ham.data.shape[0]
        self.N = N
        self.data = numpy.zeros((self.Nt,N),dtype=numpy.complex64)        


    def propagate(self, psii):
        
        return self._propagate_short_exp(psii,L=4)
        
    def get_evolution_operator(self):
        
        eop = 0.0
        return EvolutionOperator(self.timeaxis, data=eop)
        
        
    def _propagate_short_exp(self, psii, L=4):
        """
              Short exp integration
        """
        
        pr = StateVectorEvolution(self.timeaxis, psii)
        
        psi1 = psii.data
        psi2 = psii.data
        
        HH = self.ham.data        
        
        indx = 1
        for ii in self.timeaxis.data[1:self.Nt]:
            
                    
            for jj in range(0,self.Nref):
                
                for ll in range(1,L+1):
                    pref = (self.dt/ll) 
                    psi1 = -1j*pref*numpy.dot(HH,psi1)
                    psi2 = psi2 + psi1

                psi1 = psi2    
                
            pr.data[indx,:] = psi2                        
            indx += 1             
            
        return pr
