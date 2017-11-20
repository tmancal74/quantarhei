# -*- coding: utf-8 -*-
import numpy

from ..propagators.rdmpropagator import ReducedDensityMatrixPropagator
from ..hilbertspace.operators import ReducedDensityMatrix
from ...core.time import TimeAxis

class SuperOperator:
    
    def __init__(self, data=None):
        
        self.data = data
        
    def apply(self, oper):
        
        oper.data = numpy.tensordot(self.data, oper.data)
        
        return oper
    
    

class EvolutionSuperOperator:
    """Class representing evolution super operator
    
    
    
    Parameters
    ----------
    
    time : TimeAxis
        TimeAxis obejct specifying the time points of the superoperator
    
    ham: Hamiltonian
        Hamiltonian of the system
        
    relt: relaxation tensor
        Relaxation tensor of the system
        
    
    """
    
    def __init__(self, time=None, ham=None, relt=None):
        
        
        self.time = time
        self.ham = ham
        self.relt = relt
        
        self.dense_time = None
        
    def set_dense_dt(self, Nt):
        """Set a denser time axis for calculations between two points of the superoperator
        
        Parameters
        ----------
        
        Nt : int
            Number of steps between two points of the superoperator to be used
            for numerical propagation
            
        """
        self.dense_time = TimeAxis(0.0, Nt, self.time.step/Nt)
        
        
    def update_dense_time(self, i):
        """Update the start time of the dense_time
        
        
        """
        start = self.time.data[i]
        
        self.dense_time.start = start
        self.dense_time = TimeAxis(self.dense_time.start, 
                                   self.dense_time.length,
                                   self.dense_time.step) 
        
        
        
    def calculate(self):

        dim = self.ham.dim
        Nt = self.time.length
        self.data = numpy.zeros((Nt, dim, dim, dim, dim),
                                dtype=numpy.complex128)

        # zero time value (unity superoperator)
        for i in range(dim):
            for j in range(dim):
                self.data[0,i,j,i,j] = 1.0

        # first propagation
        self.update_dense_time(0)
        prop = ReducedDensityMatrixPropagator(self.dense_time, self.ham,
                                              self.relt)
        ctime = self.dense_time  
        for n in range(dim):
            for m in range(dim):
                rhonm0 = ReducedDensityMatrix(dim=dim)
                rhonm0.data[n,m] = 1.0
                rhot = prop.propagate(rhonm0)
                
                self.data[1, :, :, n, m] = rhot.data[ctime.length-1, :, :]

        
        for ti in range(1, Nt):
            self.update_dense_time(ti)
            prop = ReducedDensityMatrixPropagator(self.dense_time, self.ham,
                                                  self.relt)
            ctime = self.dense_time  
                
            # the rest of calculations
            print("evolutionsuperoperator.py: Cycle", ti, "of", Nt)
            for n in range(dim):
                for m in range(dim):
                    rhonmt = ReducedDensityMatrix(dim=dim)
                    rhonmt.data = self.data[ti-1, :, :, n, m]
                    rhot = prop.propagate(rhonmt)
                    
                    self.data[ti, :, :, n, m] = rhot.data[ctime.length-1, :, :]            

                    
                    
    def at(self, time):

        ti, dt = self.time.locate(time)

        return SuperOperator(data=self.data[ti, :, :, :, :])

                  
                