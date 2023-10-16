# -*- coding: utf-8 -*-
"""


    StateVector propagator


""" 

import numpy

from .statevectorevolution import StateVectorEvolution
from ..hilbertspace.evolutionoperator import EvolutionOperator
#from ... import REAL


   
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
        
        >>> from quantarhei import Hamiltonian, TimeAxis
        >>> HH = Hamiltonian(data=numpy.array([[0.0, 0.0],[0.0,1.0]]))
        >>> times = TimeAxis(0,1000,10.0)
        >>> pr = StateVectorPropagator(times,HH)
        >>> pr.setDtRefinement(10)
        
        """
        self.Nref = Nref
        self.dt = self.Odt/self.Nref
        
        
    def propagate(self, psii, L=4, hfce=None, nonlinear=False):
        """Propagates the state vector
        
        
        Parameters
        ----------
        
        psii : 
            Initial state vector
        
        L : int
            Order of the exponential expansion used for the propagation
            
        hfce : function or None
            A function that returns the Hamiltonian at a given time
            
        nonlinear : bool
            If True, the function will accept the state vector as its argument   
        
        
        """
        if hfce is not None:
            
            # propagation with the Hamiltonian defined through a function
            
            if nonlinear:
                # non-linear version
                return self._propagate_short_exp_nonlin(psii, hfce, L=L)
                
            else:
                # just time dependence
                return self._propagate_short_exp_tdep(psii, hfce, L=L)

        # standard propagation with time independent Hamiltonian      
        return self._propagate_short_exp(psii, L=L)
        
        
    def get_evolution_operator(self):
        """Returns the evolution operator corresponding to the propagator
        
        
        """
        
        # FIXME: not yet implemented
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

    def _propagate_short_exp_nonlin(self, psii, hfce, L=4) :
        """
              Short exp integration with non-linear "Hamiltonian"
              
        """
        
        # here we will store the results
        pr = StateVectorEvolution(self.timeaxis, psii)
        
        Nel = 3
        Nvib = numpy.zeros(Nel, dtype=int)
        Nvib[0] = 2
        Nvib[1] = 2
        Nvib[2] = 2
        
        print(psii.data.shape[0], numpy.sum(Nvib))
        
        raise Exception()

        # we create larger vectors (size + no of electronic states),
        # store the normalization constants and normalized state vectors
        
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
                
                hh = hfce(HH, psi1)
                
                for ll in range(1,L+1):
                    pref = (self.dt/ll)
                    psi1 = -1j*pref*numpy.dot(hh,psi1)
                    psi2 = psi2 + psi1

                psi1 = psi2    
                
            pr.data[indx,:] = psi2                        
            indx += 1       
            
        if self.ham.has_rwa:
            pr.is_in_rwa = True
            
        return pr


    def _propagate_short_exp_tdep(self, psii, hfce, L=4):
        """
              Short exp integration with time-dependent Hamiltonian
              
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


