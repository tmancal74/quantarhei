# -*- coding: utf-8 -*-

import numpy
import matplotlib.pylab as plt

from ...core.matrixdata import MatrixData
from ...core.time import TimeAxis

from ...utils.types import BasisManagedComplexArray
from ...core.managers import BasisManaged

class StateVectorEvolution(MatrixData, BasisManaged): 
    
    data = BasisManagedComplexArray("data")
    
    def __init__(self, timeaxis, psii):
        
        if not isinstance(timeaxis, TimeAxis):
            raise Exception
                
        self.TimeAxis = timeaxis
        self.psi_i = psii
        self._data = numpy.zeros((timeaxis.length, psii.data.shape[0]),
                                 dtype=numpy.complex128)
        self.dim = psii.data.shape[0]
        self.data[0,:] = psii.data


    def plot(self, show=True, ptype="real"):
        
        if ptype == "real":
            for i in range(self.data.shape[1]):
                plt.plot(self.TimeAxis.data, numpy.real(self.data[:,i]))
        elif ptype == "imag":
            for i in range(self.data.shape[1]):
                plt.plot(self.TimeAxis.data, numpy.imag(self.data[:,i]))
        elif ptype == "abs":
            for i in range(self.data.shape[1]):
                plt.plot(self.TimeAxis.data, numpy.abs(self.data[:,i]))
        elif ptype == "square":
            for i in range(self.data.shape[1]):
                plt.plot(self.TimeAxis.data, numpy.abs(self.data[:,i])**2)
            
            
        if show:
            plt.show()
            

    def transform(self, SS, inv=None):
        """Transformation of the operator by a given matrix
        
        
        This function transforms the Operator into a different basis, using
        a given transformation matrix.
        
        Parameters
        ----------
        
        SS : matrix, numpy.ndarray
            transformation matrix
            
        inv : matrix, numpy.ndarray
            inverse of the transformation matrix
            
        """        
        if inv is None:
            S1 = numpy.linalg.inv(SS)
        else:
            S1 = inv

        for nt in range(self.TimeAxis.length):
            self._data[nt,:] = numpy.dot(S1,self._data[nt,:])
        
            
            
    def __str__(self):
        out  = "\nquantarhei.StateVectorEvolution object"
        out += "\n======================================"
        out += "\nTimeAxis parameters: "        
        out += "\n    start = "+str(self.TimeAxis.start)
        out += "\n    length = "+str(self.TimeAxis.length)
        out += "\n    step = "+str(self.TimeAxis.step)
        out += "\nInitial StateVector data: "
        out += "\n    "+str(self.psi_i.data)
        #out += str(self.psi_i)
        out += "\ndata = \n"
        out += str(self.data)
        return out
        