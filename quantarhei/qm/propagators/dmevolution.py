# -*- coding: utf-8 -*-

import numpy
import matplotlib.pylab as plt

from ...core.matrixdata import MatrixData
from ...core.time import TimeAxis
from ...utils.types import BasisManagedComplexArray
from ...core.managers import BasisManaged
from ...core.saveable import Saveable

class DensityMatrixEvolution(MatrixData, BasisManaged, Saveable):
    
    data = BasisManagedComplexArray("data") 
    
    def __init__(self, timeaxis=None, rhoi=None, name=None):
        
        if timeaxis is not None:
            
            if name is not None:
                self.name = name
            else:
                self.name = ""
                
            self.TimeAxis = timeaxis
            
            if rhoi is not None:
                self.dim = rhoi.data.shape[1]        
                
                self.set_initial_condition(rhoi)
            else:
                self.dim = 0
            
            
    def set_initial_condition(self, rhoi):
        """
        
        
        """
        self.dim = rhoi.data.shape[1] 
        self._data = numpy.zeros((self.TimeAxis.length, self.dim, \
                    self.dim),dtype=numpy.complex128)
        self.data[0,:,:] = rhoi.data        
        

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
        if (self.manager.warn_about_basis_change):
                print("\nQr >>> DensityMatrixEvolution '%s' changes basis" %self.name)
        
        if inv is None:
            S1 = numpy.linalg.inv(SS)
        else:
            S1 = inv

        #S1 = scipy.linalg.inv(SS)                 
        for ii in range(self.TimeAxis.length):
            self._data[ii,:,:] = numpy.dot(S1,
                    numpy.dot(self._data[ii,:,:],SS))    
                    
                    
    def plot(self,populations=True,popselection="All",\
                  coherences=True, cohselection="All",how='-',
                  axis=None, show=True):
        """
            Plots selected data.
            Return figure so that it can be manipulated
        """
        population_sum = False
        #populations=False

        if how == '-':
            howi = ['-k','-r','-b','-g','-m','-y','-c']
        if how == '--':
            howi = ['--k','--r','--b','--g','--m','--y','--c']
            
        N = self.data.shape[1]

        if populations:
            for ii in range(1,N):
                kk = ii
                while kk > 6:
                    kk = kk - 7
                plt.plot(self.TimeAxis.data,
                         numpy.real(self.data[:,ii,ii]),howi[kk])
                
        
        if coherences:
            kk = 0
            for ii in range(0,N):
                for jj in range(ii+1,N):
                    if (ii != 0):
                        kk += 1
                        if kk > 6:
                            kk = kk - 7
                        plt.plot(self.TimeAxis.data,numpy.real(
                                    self.data[:,ii,jj]),howi[kk]) #how)
                
        ii = 0
        ss = numpy.zeros((self.TimeAxis.length),dtype=numpy.complex64)
        
        if population_sum:
            for t in self.TimeAxis.time:
                ss[ii] = numpy.trace(self.data[ii,:,:])
                ii +=1
        
            plt.plot(self.TimeAxis.time,numpy.real(ss),'-k')
        
        if axis is not None:
            plt.axis(axis)

        if show:
            plt.show()
        
        
    def _exportDataToText(self, file):
        """Saves textual data to a file

        """
        Nt = self.data.shape[0]
        N = self.data.shape[1]
        # all elements [real] + (all elements - diagonal) [imaginary] + time
        Na = N + 1 + N*(N-1) 

        out = numpy.zeros((Nt, Na), dtype=numpy.float64)   
        
        for i in range(Nt):
            #print("%%%%")
            # time
            out[i,0] = self.TimeAxis.data[i]
            #print(0)
            # populations
            for j in range(N):
               out[i,j+1] = numpy.real(self.data[i,j,j])
               #print(j+1)
            # coherences
            l = 0
            for j in range(N):
                for k in range(j+1,N):
                    out[i,N+1+l] = numpy.real(self.data[i,j,k])
                    #print(N+1+l)
                    l += 1
                    out[i,N+1+l] = numpy.imag(self.data[i,j,k])
                    #print(N+1+l)
                    l += 1
                    
        numpy.savetxt(file, out)


    def _importDataFromText(self, filename):
        """Imports textual data to a file

        """        

        out = numpy.loadtxt(filename)
        
        N = int(numpy.sqrt(out.shape[1] - 1))

        Nt = self.TimeAxis.length
        if Nt != out.shape[0]:
            raise Exception("Incompatibel number of time steps")
        
        din = numpy.zeros((N,N), dtype=numpy.float64)
        for i in range(N):
            din[i,i] = out[0,i+1]
            
        self.set_initial_condition(din)
        
        for i in range(Nt):
            # populations
            for j in range(N):
               self.data[i,j,j] = out[i,j+1]
            # coherences
            l = 0
            for j in range(N):
                for k in range(j+1,N):
                    self.data[i,j,k] = out[i,N+1+l] + 1.0j*out[i,N+2+l]
                    self.data[i,k,j] = numpy.conj(self.data[i,j,k])
                    l += 2
              
        
        
class ReducedDensityMatrixEvolution(DensityMatrixEvolution):
    pass

