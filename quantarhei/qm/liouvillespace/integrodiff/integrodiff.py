# -*- coding: utf-8 -*-
import types
import numpy


from ...propagators.dmevolution import ReducedDensityMatrixEvolution
from ...liouvillespace.liouvillian import Liouvillian
from ...liouvillespace.supopunity import SOpUnity
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
    
    def __init__(self, timeaxis, ham, kernel=None, cutoff_time=-1,
                 inhom=None, fft=True, save_fft_kernel=False,
                 timefac=3, decay_fraction=2.0,
                 correct_short_time=False, correction_length=0.0):
        
        self.timeaxis = timeaxis
        self.ham = ham
        self.kernel = kernel
        
        self.kernel_is_tdependent = False
        if isinstance(self.kernel, types.FunctionType):
            # kernel is updatable
            self.kernel_is_tdependent = True
        
        self._kernel = self.kernel  # we make the _kernel same as kernel if 
                                    # the kernel does not depend on time
                                    # independently
                                    
        if self.kernel_is_tdependent:
            # get the kernel at t=0.0
            self._kernel = self.kernel(0.0)
            
        
        if cutoff_time > 0 and cutoff_time <= \
                               self._kernel.shape[0]*timeaxis.step:
            self.kernel_cutoff = min(int(cutoff_time/timeaxis.step),
                                     self._kernel.shape[0])
        elif self._kernel is not None:
            self.kernel_cutoff = self._kernel.shape[0]
        else:
            self.kernel_cutoff = -1
            
        self.inhom = inhom
        
        if self.kernel_is_tdependent:
            self.fft = False
        else:
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
                   
    
            # FFT of the kernel
            N1 = self.ham.dim
             
            self.om = numpy.zeros(len(tt), REAL)
            # calculate frequencies
            self.om[:] = (2.0*numpy.pi)* \
                         numpy.fft.fftfreq(tlen, self.timeaxis.step)
            om = self.om           
            
            # Superoperator unity  
            unity = SOpUnity(dim=self.ham.dim).data                
            
            # Liouvillian
            LL = Liouvillian(self.ham).data
            
            if self._kernel is not None:
                with_kernel = True
            else:
                with_kernel = False
                
            if with_kernel:
                
                # if kernel is present, we use it for calculation
                MM = numpy.zeros((tlen, N1, N1, N1, N1), dtype=COMPLEX)
                MM[:self.timeaxis.length,:,:,:,:] = self._kernel
                
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
            #
            # Propagation by FFT method
            #
            
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
            
            #
            # propagation in time domain by short exponential approximation
            #
            
            rho1 = rhoi.data
            rho2 = rhoi.data
        
            self.Nref = 1  # we need to calculate integral from the saved
                           # values, and correspondingly, dt has to be small
            
            L = 4
            self.last_tn = -1
            dt = self.timeaxis.step
            indx = 1
            
            if self._kernel is None:
                
                #   
                # Propagation without kernel
                #
                for ii in range(1,self.timeaxis.length):
                
                    for jj in range(0, self.Nref):
                        
                        for ll in range(1, L+1):
                            
                            rho1 = (dt/ll)* \
                                   self._no_kernel_rhs(ii, rho1, rhot)
                                     
                            rho2 = rho2 + rho1
                        rho1 = rho2    
                        
                    rhot.data[indx,:,:] = rho2                        
                    indx += 1                             
                
            else:
            
                #
                # Propagation with the kernel
                #
                for ii in range(1,self.timeaxis.length):
                
                    for jj in range(0, self.Nref):
                        
                        for ll in range(1, L+1):
                            
                            rho1 = (dt/ll)* \
                                   self._right_hand_side(ii, rho1, rhot)
                                     
                            rho2 = rho2 + rho1
                        rho1 = rho2    
                        
                    rhot.data[indx,:,:] = rho2   
                   
                    indx += 1  
                    
                    #print(indx, ii)
            
            return rhot
        
        
        
    def _convolution_with_kernel(self, tn, rho_in, rhot):
        """Convolution of the density matrix with integration kernel
        
        
        """
        dim = rhot.data.shape[1]
        
        
        #self.kernel_cutoff = self._kernel.shape[0]
        
        if tn == self.last_tn:
            
            rho = self.last_int
            
        else:
            
            if self.kernel_is_tdependent:
                # updating the kernel if it depends on t independently
                self._kernel = self.kernel(tn)
                
            rho = numpy.zeros((dim, dim), dtype=COMPLEX)
            self.last_int = numpy.zeros((dim, dim), dtype=COMPLEX)
            
            for nn in range(1, min(tn+1, self.kernel_cutoff)):
                nr = tn - nn
                rho += self.timeaxis.step* \
                   numpy.tensordot(self._kernel[nn,:,:,:,:],rhot.data[nr,:,:])
            
            self.last_int[:,:] = rho
            
        rho += \
            self.timeaxis.step*numpy.tensordot(self._kernel[0,:,:,:,:],rho_in)
        
        self.last_tn = tn
        
        return rho


    def _right_hand_side(self, tn, rho, rhot):
        """ Right-hand side of the master equation 
        
        """
        ham = self.ham.data
        drho = -1j*(numpy.dot(ham, rho) - numpy.dot(rho, ham))
        drho += -self._convolution_with_kernel(tn, rho, rhot)
        
        return drho


    def _no_kernel_rhs(self, tn, rho, rhot):
        """ Right-hand side of the master equation without the kernel 
        
        """
        #M0 = self.kernel
        
        ham = self.ham.data
        drho = -1j*(numpy.dot(ham, rho) - numpy.dot(rho, ham))
        #drho += -self.timeaxis.stop*M0*numpy.sum(rhot.data[:,:,:],axis=0)
        
        return drho
        
    
    def check_solution(self, rhot):
        """Masures the difference between the right- and left-hand-side 
           of the equation
        
        """
        N1 = self.ham.dim
        diff = numpy.zeros((self.timeaxis.length, N1, N1),
                           dtype=COMPLEX)
        for tn in range(self.timeaxis.length):
            
            # time derivative
            if tn == 0:
                lhs = (rhot.data[tn+1]-rhot.data[tn])/self.timeaxis.step
            elif tn == self.timeaxis.length-1:
                lhs = (rhot.data[tn]-rhot.data[tn-1])/self.timeaxis.step
            else:
                lhs = (rhot.data[tn+1]-rhot.data[tn-1])/ \
                      (2.0*self.timeaxis.step)
                     
            # right-hand-side
            rho = rhot.data[tn, :,:]
            rhs = self._right_hand_side(tn, rho, rhot)
            
            diff[tn, :,:] = lhs - rhs
        
        return diff
    

        
             