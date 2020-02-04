# -*- coding: utf-8 -*-
"""


    StateVector propagator


""" 

import numpy

from .statevectorevolution import StateVectorEvolution
from ..hilbertspace.evolutionoperator import EvolutionOperator
from ... import REAL

   
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
        

    def setDtRefinement(self, Nref):
        """
        The TimeAxis object specifies at what times the propagation
        should be stored. We can tell the propagator to use finer
        time step for the calculation by setting the refinement. The
        refinement is an integer by which the TimeAxis time step should
        be devided to get the finer time step. In the code below, we
        have dt = 10 in the TimeAxis, but we want to calculate with
        dt = 1
        
        >>> HH = numpy.array([[0.0, 0.0],[0.0,1.0]])
        >>> times = numpy.linspace(0,1000,10)
        >>> pr = StateVectorPropagator(HH,times)
        >>> pr.setDtRefinement(10)
        
        """
        self.Nref = Nref
        self.dt = self.Odt/self.Nref
        
        
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
        
        #
        # RWA is applied here
        #
        if self.ham.has_rwa:
            HH = self.ham.get_RWA_data()
        else:
            HH = self.ham.data      
        
        indx = 1
        for ii in range(1,self.Nt):
            
            for jj in range(0,self.Nref):
                
                for ll in range(1,L+1):
                    pref = (self.dt/ll)
                    psi1 = -1j*pref*numpy.dot(HH,psi1)
                    psi2 = psi2 + psi1

                psi1 = psi2    
                
            pr.data[indx,:] = psi2                        
            indx += 1       
            
        if self.ham.has_rwa:
            pr.is_in_rwa = True
            
        return pr
