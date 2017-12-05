# -*- coding: utf-8 -*-
"""
*******************************************************************************


    QUANTArhei: Open Quantum System Theory for Molecular Systems 
    ============================================================
    
    (c) 2016 Tomáš Mančal
    
    Charles University
    Faculty of Mathematics and Physics
    Ke Karlovu 5
    CZ-121 16 Prague 2
    Czech Repubic


    For support contact the author at : mancal@karlov.mff.cuni.cz
    
    
*******************************************************************************

    Evolution Superoperator module

"""

# standard library imports
import time

# dependencies imports
import numpy

# quantarhei imports
from ..propagators.rdmpropagator import ReducedDensityMatrixPropagator
from ..hilbertspace.operators import ReducedDensityMatrix
from ...core.time import TimeAxis
from ...core.saveable import Saveable

import quantarhei as qr

# FIXME: This class should be a base class for Relaxation tensors
class SuperOperator:
    """Class representing superoperators
    
    
    This class represents operators on the space of Hilbert space operators.
    Usually, we refer to such operators as superoperators. 
    
    Parameters
    ----------
    
    dim : int
        Dimension of the superoperator
        
    data : array
        Data of the superoperator
        
    real : bool
        Is this data real? False if they are complex
    
    """
    
    def __init__(self, dim=None, data=None, real=False):
        
        if dim is not None:
            self.dim = dim
            if real:
                self.data = numpy.zeros((dim, dim, dim, dim),
                                        dtype=qr.REAL)
            else:
                self.data = numpy.zeros((dim, dim, dim, dim), 
                                        dtype=qr.COMPLEX)
        elif data is not None:
            self.data = data
            self.dim = data.shape[0]
        
    def apply(self, oper):
        
        oper.data = numpy.tensordot(self.data, oper.data)
        
        return oper
    
    
class SOpUnity(SuperOperator):
    """Class representing a unity superoperator
    
    
    Parameters
    ----------
    
    dim : int
        Dimension of the unity superoperator
        
    data : array
        If data is specified, only their dimension is used to construct
        a unity superoperator
        
    """
    
    def __init__(self, dim=None, data=None):
        
        if dim is not None:
            super().__init__(dim, data=None)
            
        elif data is not None:
            dim = data.shape[0]
            super().__init__(dim, data=None)
            
        # initialize the data
        for i in range(self.dim):
            for j in range(self.dim):
                self.data[i,j,i,j] = 1.0
            
        
        
        
    

class EvolutionSuperOperator(Saveable):
    """Class representing evolution superoperator
    
    
    
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
        
    
    def _estimate_remaining_loops(self, Nt, dim, ti, n, m):
        """
        
        """
        return (dim-n)*dim 
        
    def _init_progress(self):

        self.ccount = 0
        self.strtime = time.time()
        self.oldtime = self.strtime
        self.remlast = 0
        
    def _progress(self, Nt, dim, ti, n, m):
        
        curtime = time.time()
        #dt = curtime - self.oldtime
        Dt = curtime - self.strtime
        self.oldtime = curtime
        self.ccount += 1
        
        dt_past = Dt/(self.ccount*dim)
        
        remloops = self._estimate_remaining_loops(Nt, dim, ti, n, m)
        
        remtime = int(remloops*(dt_past))
        
        txt1 = "Propagation cycle "+str(ti)+" of "+str(Nt)+ \
               " : ("+str(n)+" of "+str(dim)+")"
        if True:
        #if remtime <= self.remlast:
            txt2 = " : remaining time "+str(remtime)+" sec"
        #else:
        #    txt2 = " : remaining time ... calculating"
        tlen1 = len(txt1)
        tlen2 = len(txt2)
        
        rem = 65 - tlen1 - tlen2
        txt = "\r"+txt1+(" "*rem)+txt2
        print(txt, end="\r")
        
        self.remlast = remtime
        
        
    def calculate(self, show_progress=False):
        """Calculates the data of the evolution superoperator
        
        
        """

        dim = self.ham.dim
        Nt = self.time.length
        self.data = numpy.zeros((Nt, dim, dim, dim, dim),
                                dtype=qr.COMPLEX)

        # zero time value (unity superoperator)
        for i in range(dim):
            for j in range(dim):
                self.data[0,i,j,i,j] = 1.0
            
        if show_progress:
            print("Calculating evolution superoperator ")
        self._init_progress()
        
        # first propagation
        self.update_dense_time(0)
        prop = ReducedDensityMatrixPropagator(self.dense_time, self.ham,
                                              self.relt)
        ctime = self.dense_time  

        for n in range(dim):
            if show_progress:
                self._progress(Nt, dim, 0, n, 0)
            for m in range(dim):
                rhonm0 = ReducedDensityMatrix(dim=dim)
                rhonm0.data[n,m] = 1.0
                rhot = prop.propagate(rhonm0)
                
                self.data[1, :, :, n, m] = rhot.data[ctime.length-1, :, :]

        Udt = self.data[1,:,:,:,:]
        for ti in range(1, Nt):
            if show_progress:
                print("Self propagation: ", ti, "of", Nt)            
           
            self.data[ti,:,:,:,:] = numpy.tensordot(Udt,self.data[ti-1,:,:,:,:])
            
            
            
            if False:
                self.update_dense_time(ti)
                prop = ReducedDensityMatrixPropagator(self.dense_time, self.ham,
                                                      self.relt)
                ctime = self.dense_time  
                    
                # the rest of calculations
                for n in range(dim):
                    self._progress(Nt, dim, ti, n, 0)
                    for m in range(dim):
                        
                        rhonmt = ReducedDensityMatrix(dim=dim)
                        rhonmt.data = self.data[ti-1, :, :, n, m]
                        rhot = prop.propagate(rhonmt)
                        
                        self.data[ti, :, :, n, m] = rhot.data[ctime.length-1, :, :]            

        if show_progress:
            print("...done")
                    
                    
    def at(self, time):

        ti, dt = self.time.locate(time)

        return SuperOperator(data=self.data[ti, :, :, :, :])

                  
                