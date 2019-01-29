# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot as plt

from ...propagators.dmevolution import ReducedDensityMatrixEvolution
from ...hilbertspace.operators import ProjectionOperator
from ...liouvillespace.systembathinteraction import SystemBathInteraction
from ...liouvillespace.lindbladform import LindbladForm
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

    fft : bool
        If True, propagation is solved through Fourier transform method
        Default is True
    
    save_fft_kernel : bool
        If True, the FFTed integration kernel is saved
        
    timefac : int
        Factor by which the time of the FTT is lengthed with respect 
        to the submitted TimeAxis. The added part of the time axis
        is zero-padded. Default values is 3.
        
    decay_fraction : float
        Fraction of the TimeAxis in which the decay should be e^-1. When 
        the max time of the TimeAxis is Tmax, and decay_fraction=2, then
        the decay time is Tmax/decay_fraction. Default values if 2.0
        
    correct_short_time : bool
        If set True, short time propagation will be calculated in time
        domain. Frequency domain calculation usually shows some oscillations
        around time zero.
    
    correction_length : float
        How long the correction should be (from zero)

    """
    
    def __init__(self, timeaxis, ham, kernel=None,
                 inhom=None, fft=True, save_fft_kernel=False,
                 timefac=3, decay_fraction=2.0,
                 correct_short_time=False, correction_length=0.0):
        
        self.timeaxis = timeaxis
        self.ham = ham
        self.kernel = kernel
        self.inhom = inhom
        self.fft = fft
        
        if self.fft:

            self.resolv = None
            
            # FFT on timeaxis twice as long as defines (we add negative times)
            tlen = timefac*self.timeaxis.length
            tstart = self.timeaxis.data[0]
            tstop = tstart + self.timeaxis.step*(tlen-1)
            
            tt = numpy.linspace(tstart, tstop, tlen)
            
            gamma = decay_fraction/self.timeaxis.data[self.timeaxis.length-1]
            self.gamma = gamma
            
            # test of the FFT
            test = False
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
            # calculate frequencies
            self.om[:] = (2.0*numpy.pi)* \
                         numpy.fft.fftfreq(tlen, self.timeaxis.step)
            om = self.om           
            #self.om[:] = numpy.fft.fftfreq(tlen, self.timeaxis.step)
            
            unity = numpy.zeros((N1, N1, N1, N1), dtype=REAL)
            
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
                for tm in range(len(tt)):
                    MM[tm,:,:,:,:] = MM[tm,:,:,:,:]*numpy.exp(-gamma*tt[tm])
                    
                MM = numpy.fft.ifft(MM, axis=0)*self.timeaxis.step*tlen*2.0
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
            
            rhOm = numpy.zeros((Nt, N1, N1), dtype=COMPLEX)
            
            rho0 = rhoi.data #.reshape(N1**2)
            
            for ii in range(Nt):
                G = self.resolv[ii,:,:,:,:]
                #Gr = G.reshape(N1**2, N1**2)
                rhOm[ii,:,:] = numpy.tensordot(G, rho0)            
            
            rhOm = numpy.fft.fft(rhOm, axis=0) \
                    *((self.om[1]-self.om[0])/(2.0*numpy.pi))
    
            # lift the gamma damping
            Nt = self.timeaxis.length
            for ii in range(1,Nt):
                rhot.data[ii,:,:] = rhOm[ii,:,:]* \
                numpy.exp(self.gamma*self.timeaxis.data[ii])
            
            return rhot
        
        else:
            
            # propagation in time domain

            rho1 = rhoi.data
            rho2 = rhoi.data
        
            self.Nref = 1
            L = 4
            dt = self.timeaxis.step
            indx = 1
            for ii in range(1,self.timeaxis.length):
            
                for jj in range(0, self.Nref):
                    
                    for ll in range(1,L+1):
                        
                        rho1 = (dt/ll)* \
                               self._right_hand_side(ii-1, rho1, rhot)
                                 
                        rho2 = rho2 + rho1
                    rho1 = rho2    
                    
                rhot.data[indx,:,:] = rho2                        
                indx += 1             
            
            return rhot
        
        
        
    def _convolution_with_kernel(self, tn, rhot):
        """Convolution of the density matrix with integration kernel
        
        
        """
        dim = rhot.data.shape[1]
        rho = numpy.zeros((dim, dim), dtype=COMPLEX)
        for nn in range(tn):
            nr = tn - nn
            rho += self.timeaxis.step* \
                   numpy.tensordot(self.kernel[nn,:,:,:,:],rhot.data[nr,:,:])
            
        return rho


    def _right_hand_side(self, tn, rho, rhot):
        
        ham = self.ham.data
        drho = -1j*(numpy.dot(ham, rho) - numpy.dot(rho, ham))
        #drho = numpy.tensordot(self.kernel, rho)
        drho += -self._convolution_with_kernel(tn, rhot)
        
        return drho

        
    def _test_kernel(self):
        """Returns a simple kernel for tests
        
        """
        
        dim = self.ham.dim
        Nt = self.timeaxis.length
        MM = numpy.zeros((Nt, dim, dim, dim, dim), dtype=COMPLEX)
        gamma = 1.0/10.0
        
        if dim == 2:
            K01 = ProjectionOperator(0,1,dim)
            K10 = ProjectionOperator(1,0,dim)
            
            sys_ops = [K01, K10]
            rates = [1.0/30.0, 1.0/20.0]
            sbi = SystemBathInteraction(sys_operators=sys_ops, rates=rates)
            
            lbf = LindbladForm(self.ham, sbi, as_operators=False)
            #return lbf.data
            
            for ti in range(Nt):
                tm = self.timeaxis.data[ti]
                MM[ti,:,:,:,:] = -lbf.data*numpy.exp(-gamma*tm)
                
            return MM
            