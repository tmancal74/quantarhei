# -*- coding: utf-8 -*-

import numpy
import matplotlib.pylab as plt

from ...core.matrixdata import MatrixData
from ...core.time import TimeAxis

from ...utils.types import BasisManagedComplexArray
from ...core.managers import BasisManaged

from .dmevolution import DensityMatrixEvolution
from ..hilbertspace.operators import DensityMatrix

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


    def convert_from_RWA(self, ham, sgn=1):
        """Converts density matrix evolution from RWA to standard repre
        
        
        Parameters
        ----------
        
        ham : qr.Hamiltonian
            Hamiltonian with respect to which we construct RWA
            
        sgn : {1, -1}
            Forward (1) or backward (-1) conversion. Default sgn=1 corresponds
            to the function name. Backward conversion sgn=-1 is called from
            the inverse routine.
        """
        
        if (self.is_in_rwa and sgn == 1) or sgn == -1:
            
            HOmega = ham.get_RWA_skeleton()
            
            for i, t in enumerate(self.TimeAxis.data):
                # evolution operator
                Ut = numpy.exp(-sgn*1j*HOmega*t)
                # revert RWA
                rhot = numpy.dot(Ut,self.data[i,:])
                self.data[i,:] = rhot
                
        if sgn == 1:
            self.is_in_rwa = False


    def convert_to_RWA(self, ham):
        """Converts density matrix evolution from standard repre to RWA


        Parameters
        ----------
        
        ham : qr.Hamiltonian
            Hamiltonian with respect to which we construct RWA
            
        """        
        if not self.is_in_rwa:
            self.convert_from_RWA(ham, sgn=-1)
            self.is_in_rwa = True


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


    def get_DensityMatrixEvolution(self):
        """Constructs DensityMatrix from the present StateVector
        
        """
        
        rhot = DensityMatrixEvolution(timeaxis=self.TimeAxis)
        
        rhoi = DensityMatrix(dim=self.dim)
        for ii in range(self.dim):
            for jj in range(self.dim):
                rhoi.data[ii,jj] = self.data[0,ii]* \
                                   numpy.conj(self.data[0,jj])
        
        rhot.set_initial_condition(rhoi)
        
        #for tt in range(1, self.TimeAxis.length):
        for ii in range(self.dim):
            for jj in range(self.dim):
                rhot.data[1:,ii,jj] = self.data[1:,ii]* \
                                      numpy.conj(self.data[1:,jj])
        
        return rhot
        
            
            
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
        