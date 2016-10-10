# -*- coding: utf-8 -*-

import numpy
import matplotlib.pylab as plt

from ...core.matrixdata import MatrixData
from ...core.time import TimeAxis


class StateVectorEvolution(MatrixData): 
    
    
    def __init__(self, timeaxis, psii):
        
        if not isinstance(timeaxis, TimeAxis):
            raise Exception
            
            
        self.TimeAxis = timeaxis
        self.psi_i = psii
        self.data = numpy.zeros((timeaxis.length, psii.data.shape[0]),
                                 dtype=numpy.complex128)
        self.data[0,:] = psii.data


    def plot(self):
        
        for i in range(self.data.shape[1]):
            plt.plot(self.TimeAxis.data, numpy.real(self.data[:,i]))
            
            
            
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
        