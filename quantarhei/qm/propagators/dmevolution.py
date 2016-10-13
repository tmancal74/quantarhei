# -*- coding: utf-8 -*-

import numpy
import matplotlib.pylab as plt

from ...core.matrixdata import MatrixData
from ...core.time import TimeAxis
from ...utils.types import BasisManagedComplexArray
from ...core.managers import BasisManaged

class DensityMatrixEvolution(MatrixData, BasisManaged):
    
    data = BasisManagedComplexArray("data") 
    
    def __init__(self, timeaxis, rhoi, name=""):
        
        if not isinstance(timeaxis, TimeAxis):
            raise Exception("First argument has to be a TimeAxis")
            
        #if not isinstance(rhoi,ReducedDensityMatrix):
        #    raise Exception
        self.name = name
            
        self.TimeAxis = timeaxis
        self.dim = rhoi.data.shape[1]        
        
        self._data = numpy.zeros((timeaxis.length, rhoi.data.shape[0], \
                    rhoi.data.shape[1]),dtype=numpy.complex128)
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
                  coherences=True, cohselection="All",how='-'):
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
        

        return 1
        
    def set_name(self, name):
        self.name = name
        
        
class ReducedDensityMatrixEvolution(DensityMatrixEvolution):
    pass

