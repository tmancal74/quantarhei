# -*- coding: utf-8 -*-

import numpy
import matplotlib.pylab as plt

from ...core.matrixdata import MatrixData
#from ...core.time import TimeAxis
from ...utils.types import BasisManagedComplexArray
from ...core.managers import BasisManaged
from ...core.saveable import Saveable
from ..hilbertspace.operators import DensityMatrix
from ..hilbertspace.operators import ReducedDensityMatrix

class DensityMatrixEvolution(MatrixData, BasisManaged, Saveable):
    
    data = BasisManagedComplexArray("data") 
    
    def __init__(self, timeaxis=None, rhoi=None, is_in_rwa=False, name=None):
        
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
            
        self.is_in_rwa = is_in_rwa
                                    
            
    def set_initial_condition(self, rhoi):
        """
        
        
        """
        self.dim = rhoi.data.shape[1] 
        self._data = numpy.zeros((self.TimeAxis.length, self.dim, \
                    self.dim),dtype=numpy.complex128)
        self.data[0,:,:] = rhoi.data        

        
    def at(self, time):
        """Returns density matrix at a given time
        
        
        Parameters
        ----------
        
        time : float
            Time (in fs) at which the tensor should be returned
            
            
        """

        ti, dt = self.TimeAxis.locate(time)

        return DensityMatrix(data=self.data[ti, :, :])


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
            
    
    def normalize(self):
        
        tr = numpy.trace(self.data[0,:,:])
        self.data = self.data/tr


    def get_TimeAxis(self):
        """Returns the TimeAxis of the present evolution
        
        """
        return self.TimeAxis              
                    
    
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
                Ut = numpy.diag(numpy.exp(-sgn*1j*HOmega*t))
                # revert RWA
                rhot = numpy.dot(Ut,numpy.dot(self.data[i,:,:],
                                              numpy.conj(Ut)))
                self.data[i,:,:] = rhot
                
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
          
    
    def plot(self, populations=True, popselection="All", trace=False,
                   coherences=True, cohselection="All", how='-',
                   axis=None, show=True, legend=False, start=1):
        """
            Plots selected data.
            Return figure so that it can be manipulated
        """
        population_sum = False
        #populations=False

        if how == '-':
            howi = ['-k','-r','-b','-g','-m','-y','-c']
        if how == '--':
            howi = ['--k','--r','--b','--g','--m','--y','--c',]
            
        N = self.data.shape[1]
        #print("Plotting", N, "states")
        leg = []

        if populations:
            for ii in range(start,N):
                kk = ii
                while kk > 6:
                    kk = kk - 7
                plt.plot(self.TimeAxis.data,
                         numpy.real(self.data[:,ii,ii]),howi[kk])
                leg.append(str(ii))
                
        
        if trace:
            trc = numpy.zeros(self.TimeAxis.length, dtype=numpy.float64)
            for ti in range(self.TimeAxis.length):
                trc[ti] = numpy.trace(numpy.real(self.data[ti, :, :]))
                
            plt.plot(self.TimeAxis.data, trc, "-k")
             
        
        if coherences:
            kk = 0
            for ii in range(N):
                for jj in range(ii+1,N):
                    if (ii >= start):
                        kk += 1
                        if kk > 6:
                            kk = kk - 7
                        plt.plot(self.TimeAxis.data,numpy.real(
                                    self.data[:,ii,jj]),howi[kk]) #how)
                        leg.append(str(ii)+"-"+str(jj))
                
        ii = 0
        ss = numpy.zeros((self.TimeAxis.length),dtype=numpy.complex64)
        
        if population_sum:
            for t in self.TimeAxis.time:
                ss[ii] = numpy.trace(self.data[ii,:,:])
                ii +=1
        
            plt.plot(self.TimeAxis.time,numpy.real(ss),'-k')
 

        if legend:
            plt.legend(leg)
            
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
    
    def at(self, time):
        """Returns density matrix at a given time
        
        
        Parameters
        ----------
        
        time : float
            Time (in fs) at which the tensor should be returned
            
            
        """

        ti, dt = self.TimeAxis.locate(time)

        return ReducedDensityMatrix(data=self.data[ti, :, :])

