# -*- coding: utf-8 -*-
import time
import numpy
#from numba import njit

from .oqssvevolution import OQSStateVectorEvolution
from quantarhei.qm import TDRedfieldRateMatrix
import quantarhei as qr



class OQSStateVectorPropagator():
    
    
    def __init__(self, timeaxis=None, current_matrix=None, agg=None,
                 theory_type="Redfield"):
        
        
        self.TimeAxis = timeaxis
        self.Nt = self.TimeAxis.length
        self.Odt = self.TimeAxis.data[1]-self.TimeAxis.data[0]
        self.dt = self.Odt
        self.Nref = 1     
        self.KK = current_matrix
        
        self.theory_type = theory_type
        self.agg = agg
    

    def setDtRefinement(self, Nref):
        """
        The TimeAxis object specifies at what times the propagation
        should be stored. We can tell the propagator to use finer
        time step for the calculation by setting the refinement. The
        refinement is an integer by which the TimeAxis time step should
        be devided to get the finer time step. In the code below, we
        have dt = 10 in the TimeAxis, but we want to calculate with
        dt = 1
        
       # >>> from quantarhei import TimeAxis
       # >>> KK =numpy.array([[0.0, 0.0],[0.0,1.0]], dtype=float)
       # >>> times = TimeAxis(0,1000,10.0)
       # >>> pr = OQSStateVectorPropagator(times, KK)
       # >>> pr.setDtRefinement(10)
        
        """
        self.Nref = Nref
        self.dt = self.Odt/self.Nref
    
    
    def propagate(self, psii, L=4):
        
       """Short expansion of an exponention to integrate equations of motion
       
       
       Propagation with Hamiltonian only
       
       
       """

       pr = OQSStateVectorEvolution(self.TimeAxis, psii)
       
       _CALC(self.KK, psii.data, self.Nt, self.Nref, L, self.dt, pr.data)
       
       return pr
   
    
    def get_secular_dynamics(self, psii, L=4):
       """Here we obtain the secular dynamics for self-consistent methodology
       
       """
       
       pr = OQSStateVectorEvolution(self.TimeAxis, psii)
       
       if self.theory_type == "Redfield":
           
           ham = self.agg.get_Hamiltonian()
           sbi = self.agg.get_SystemBathInteraction()
           
           
           RR = TDRedfieldRateMatrix(ham, sbi)
       
           ppi = psii.data*psii.data 
           
           
           ppt = numpy.zeros((self.Nt, ham.dim), dtype=qr.REAL)
           ppt[0,:] = ppi
               
           
           #RR0 = numpy.zeros((ham.dim, ham.dim, self.Nt))
           #for it in range(self.Nt):
           #    RR0[:,:,it] = RR.data[it,:,:]
           
           #t1 = time.time()
           _CALC(RR.data, ppi, self.Nt, self.Nref, L, self.dt, ppt)
           #_CALC(RR0, ppi, self.Nt, self.Nref, L, self.dt, ppt)
           #t2 = time.time()
           #print("Done in ", t2-t1, "sec")
           pr.data = numpy.sqrt(ppt)
           
           return pr
           

    def get_c0(self, psii, L=4):
        """Returns the zero's order solution to the iterative problem 
        
        """
        return self.get_secular_dynamics(psii, L)

        
#@njit(cache=True)   
def _CALC(KK, psii, Nt, Nref, L, dt, out):
    """Propagation numerics 
    
    """
    
    psi1 = psii
    psi2 = psii
    
    indx = 1
    for ii in range(1, Nt):
        
        for jj in range(Nref):
            
            for ll in range(1,L+1):
               
                psi1 = (dt/ll)*numpy.dot(KK[ii, :,:],psi1)
                         
                psi2 += psi1
            psi1 = psi2    
            
        out[indx,:] = psi2                        
        indx += 1                       
       