# -*- coding: utf-8 -*-
import numpy

from ..propagators.rdmpropagator import ReducedDensityMatrixPropagator
from ..hilbertspace.operators import ReducedDensityMatrix
from ...core.time import TimeAxis

class SuperOperator:
    
    def __init__(self, dim=None, data=None, real=False):
        
        if dim is not None:
            self.dim = dim
            if real:
                self.data = numpy.zeros((dim, dim, dim, dim),
                                        dtype=numpy.float64)
            else:
                self.data = numpy.zeros((dim, dim, dim, dim), 
                                        dtype=numpy.complex128)
        elif data is not None:
            self.data = data
            self.dim = data.shape[0]
        
    def apply(self, oper):
        
        oper.data = numpy.tensordot(self.data, oper.data)
        
        return oper
    
    
class SOpUnity(SuperOperator):
    
    def __init__(self, dim=None, data=None):
        
        if dim is not None:
            super().__init__(dim, data=None)
            
        elif data is not None:
            dim = data.shape[0]
            super().__init__(dim, data=None)
            
        for i in range(self.dim):
            for j in range(self.dim):
                self.data[i,j,i,j] = 1.0
            
        
        
        
    

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
        
    
    def _estimate_propagation_time(self, dt, Nt, dim, ti, n, m):
        
        #tottime = dt*Nt*(dim**2)
        
        remtime = dt*(Nt-ti)*(dim**2) + (dim-n)*(dim-m)*dt
        return remtime
        
        
    def calculate(self):
        import time

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
        oldtime = time.time()
        for n in range(dim):
            for m in range(dim):
                curtime = time.time()
                dt = curtime - oldtime
                oldtime = curtime
                remtime = self._estimate_propagation_time(dt, Nt, dim, 0, n, m)
                print("First propagation: ", n, m, ": remaining time ", remtime)
                rhonm0 = ReducedDensityMatrix(dim=dim)
                rhonm0.data[n,m] = 1.0
                rhot = prop.propagate(rhonm0)
                
                self.data[1, :, :, n, m] = rhot.data[ctime.length-1, :, :]

        
        for ti in range(1, Nt):
            print("Time step no: ", ti)
            self.update_dense_time(ti)
            prop = ReducedDensityMatrixPropagator(self.dense_time, self.ham,
                                                  self.relt)
            ctime = self.dense_time  
                
            # the rest of calculations
            print("evolutionsuperoperator.py: Cycle", ti, "of", Nt)
            for n in range(dim):
                for m in range(dim):
                    print("Propagation: ", n, m)
                    rhonmt = ReducedDensityMatrix(dim=dim)
                    rhonmt.data = self.data[ti-1, :, :, n, m]
                    rhot = prop.propagate(rhonmt)
                    
                    self.data[ti, :, :, n, m] = rhot.data[ctime.length-1, :, :]            

                    
                    
    def at(self, time):

        ti, dt = self.time.locate(time)

        return SuperOperator(data=self.data[ti, :, :, :, :])

                  
                