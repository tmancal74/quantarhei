# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot as plt

from ...propagators.dmevolution import ReducedDensityMatrixEvolution
from .... import COMPLEX, REAL

class IntegrodiffPropagator:
    """Solver of integrodifferential equations
    
    
    Parameters
    ----------
    
    timeaxis : TimeAxis
        Time axis on which the solution should be defined
        
    ham : Hamiltonian
        Hamiltonian of the problem
        
    kernel : Time dependent tensor 
        Time dependent tensor of the integration kernel
    
    inhom : Time dependent operator
        Inhomogeneous terms - a time dependent operator


    """
    
    def __init__(self, timeaxis, ham, kernel=None,
                 inhom=None, fft=True, save_fft_kernel=False):
        
        self.timeaxis = timeaxis
        self.ham = ham
        self.kernel = kernel
        self.inhom = inhom
        self.fft = fft
        
        if self.fft:

            self.resolv = None
            
            # FFT on timeaxis twice as long as defines (we add negative times)
            tlen = 2*self.timeaxis.length
            tstart = self.timeaxis.data[0]
            tstop = self.timeaxis.step*(tlen-1)
            
            tt = numpy.linspace(tstart, tstop, tlen)
            
            gamma = 1.0/self.timeaxis.data[self.timeaxis.length-1]
            self.gamma = gamma
            
            # test of the FFT
            test = True
            if test:
                
                gam = 1.0/5.0
                
                ff = numpy.cos(2.0*numpy.pi*tt/3.0)
                plt.plot(tt, ff)
                plt.show()
                
                fcor = ff*numpy.exp(-gam*tt)
                fo = numpy.fft.fft(fcor)
                
                plt.plot(tt, fo)
                plt.show()
                
                f1cor = numpy.fft.ifft(fo)
                f1 = f1cor*numpy.exp(gam*tt)
        
                plt.plot(tt, ff, "-b")
                plt.plot(tt, f1, "--r")
                plt.show()        
    
            # FFT of the kernel
            N1 = self.ham.dim
             
            self.om = numpy.zeros(len(tt), REAL)
            om = self.om
            
            # calculate frequencies
            self.om[:] = 1.0
            
            unity = numpy.zeros(( N1, N1, N1, N1), dtype=REAL)
            
            # Superoperator unity                  
            for ii in range(N1):
                for jj in range(N1):
                    for kk in range(N1):
                        for ll in range(N1):
                            if (ii == kk) and (jj == ll):
                                unity[ii,jj,kk,ll] = 1.0
             
            # Kronecker delta                  
            delta = numpy.zeros((N1, N1), dtype=REAL)
            for ii in range(N1):
                delta[ii,ii] = 1.0
            
            # Liouvillian
            LL = numpy.zeros((N1, N1, N1, N1), dtype=REAL)
            for ii in range(N1):
                for jj in range(N1):
                    for kk in range(N1):
                        for ll in range(N1):
                            LL[ii,jj,kk,ll] = ham.data[ii,kk]*delta[jj,ll] \
                                             -ham.data[ll,jj]*delta[kk,ii]
            if self.kernel is not None:
                with_kernel = True
            else:
                with_kernel = False
                
            if with_kernel:
                # if kernel is present, we use it for calculation
                MM = numpy.zeros((tlen, N1, N1, N1, N1), dtype=COMPLEX)
                MM[:self.timeaxis.length,:,:,:,:] = self.kernel
                
                # we renormalize the kernel by a e^{-gam*t} decay
                for tm in tt:
                    MM[tm,:,:,:,:] = MM[tm,:,:,:,:]*numpy.exp(-gamma*tm)
                    
                MM = numpy.fft.fft(MM, axis=0)
                if save_fft_kernel:
                    self.fftKernel = MM
                    
                # this is now going over frequencies 
                self.resolv = numpy.zeros(MM.shape, COMPLEX)
                for io in range(len(tt)):
                    self.resolv[io,:,:,:,:] = (-1j*om[io] + gamma)*unity \
                                               +1j*LL + MM[io,:,:,:,:]
                    A = self.resolv[io,:,:,:,:].reshape(N1**2, N1**2)
                    A = numpy.linalg.inv(A)                          
                    self.resolv[io,:,:,:,:] = A.reshape(N1, N1, N1, N1)
            else:
                # if kernel is not present, we calculate with Hamiltonian only
                self.resolv = numpy.zeros((tlen, N1, N1, N1, N1), COMPLEX)
                for io in range(len(tt)):
                    self.resolv[io,:,:,:,:] = (-1j*om[io] + gamma)*unity \
                                               +1j*LL
                    A = self.resolv[io,:,:,:,:].reshape(N1**2, N1**2)
                    A = numpy.linalg.inv(A)
                    self.resolv[io,:,:,:,:] = A.reshape(N1, N1, N1, N1)                          
            
        else:
            
            # prepare propagation in time domain
            pass
        

    def propagate(self, rhoi):
        """Propagates initial condition of an integro-differential eq.


        """
        N1 = rhoi.data.shape[0]
        rhot = ReducedDensityMatrixEvolution(self.timeaxis, rhoi)

        if self.fft:
            
            Nt = len(self.om)
            
            rhOm = numpy.zeros((Nt, N1**2), dtype=COMPLEX)
            
            rho0 = rhoi.data.reshape(N1**2)
            
            for ii in range(Nt):
                G = self.resolv[ii,:,:,:,:]
                Gr = G.reshape(N1**2, N1**2)
                rhOm[ii,:] = numpy.dot(Gr, rho0)
                
            rhOm = numpy.fft.ifft(rhOm, axis=0)
    
            tlen = 2*self.timeaxis.length
            tstart = self.timeaxis.data[0]
            tstop = self.timeaxis.step*(tlen-1)
            
            tt = numpy.linspace(tstart, tstop, tlen)
    
            # lift the gamma damping
            for ii in range(self.timeaxis.length):
                rhot.data[ii,:,:] = rhOm[ii,:].reshape(N1, N1)* \
                                    numpy.exp(self.gamma*tt[ii])
            
            return rhot
        
        else:
            
            # propagation in time domain
            
            return rhot
        
        
        
        
