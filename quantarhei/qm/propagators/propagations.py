# -*- coding: utf-8 -*-

import numpy
import matplotlib.pylab as plt

from ...data.matrixdata import MatrixData
from ...domains.time import TimeAxis


class DMPropagation(MatrixData): #,TimeDependent):
    
    
    def __init__(self,timeaxis,rhoi):
        
        if not isinstance(timeaxis,TimeAxis):
            raise Exception
            
        #if not isinstance(rhoi,ReducedDensityMatrix):
        #    raise Exception
            
        self.TimeAxis = timeaxis
        self.data = numpy.zeros((timeaxis.length, rhoi.data.shape[0], \
                    rhoi.data.shape[1]),dtype=numpy.complex128)
        self.data[0,:,:] = rhoi.data
        
        
    def transform(self,SS):
        S1 = SS.T
        for ii in range(len(self.TimeAxis.data)):
            self.data[ii,:,:] = numpy.dot(S1,
                    numpy.dot(self.data[ii,:,:],SS))
                    
    
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
                plt.plot(self.TimeAxis.time,
                         numpy.real(self.data[:,ii,ii]),howi[kk])
                
        
        if coherences:
            kk = 0
            for ii in range(0,N):
                for jj in range(ii+1,N):
                    if (ii != 0):
                        kk += 1
                        if kk > 6:
                            kk = kk - 7
                        plt.plot(self.TimeAxis.time,numpy.real(
                                    self.data[:,ii,jj]),howi[kk]) #how)
                
        ii = 0
        ss = numpy.zeros((self.TimeAxis.length),dtype=numpy.complex64)
        
        if population_sum:
            for t in self.TimeAxis.time:
                ss[ii] = numpy.trace(self.data[ii,:,:])
                ii +=1
        
            plt.plot(self.TimeAxis.time,numpy.real(ss),'-k')
        

        return 1
        
        
