# -*- coding: utf-8 -*-
import numpy
import scipy

from .redfieldtensor import RedfieldRelaxationTensor
from ...core.time import TimeDependent

from ...core.managers import Manager
from ... import REAL

class TDRedfieldRelaxationTensor(RedfieldRelaxationTensor, TimeDependent):

    # FIXME: mimick the time-independent case
    Lm = None  # we isolate the operators defined by inheritance as time-independent
    Ld = None
    Km = None
    
    def _implementation(self, ham, sbi):
        """ Reference implementation, completely in Python
        
        Implementation of Redfield relaxation tensor according to 
        
        V. May and O. Kuehn, Charge and Energy Transfer Dynamics in Molecular
        System, Wiley-VCH, Berlin, 2000, 1st edition, chapter 3.8.2.
        In particular we refer to Eq. (3.8.13) on page 132
        
        We assume the system-bath interaction operator in form of Eq. (3.5.30)
        with the bath part specified through two-point correlation functions.
        We construct operators K_{m} introduced in Eq. (3.5.30) and 
        operators \Lambda_{m} of Eq. (3.8.11). 
        
        We do not delete the imaginary part of the tensor (as is done later 
        in Section 3.8.3 to get the so-called "Multi-level Redfield
        Equations"). Such a deletion can be done later manually. 
        
        
        """
        #
        # dimension of the Hamiltonian (includes excitons
        # with all multiplicities specified at its creation)
        #
        Na = ham.dim #data.shape[0]
        
        # time axis
        ta = sbi.TimeAxis
        
        #
        # is this beyond single excitation band?
        #
        multi_ex = False
        
        # figure out if the aggregate specifies more than one exciton band
        if sbi.aggregate is not None:
            agg = sbi.aggregate
            if agg.mult > 1:
                multi_ex = True 
        
        #
        # shorten the interval of integration if a cut-off time is set
        #
        if self._has_cutoff_time:
            # index of the cut-off time on the time axis
            tcut = ta.nearest(self.cutoff_time)
            # select the section of the time axis up to the cut-off time
            tm = ta.data[0:tcut]
            # length of the section corresponds to the index of cut-off time
            length = tcut
        else:
            # if cut-off time is not set then we take the whole time axis
            tm = ta.data
            # and the length corresponds to the length of the time axis
            length = ta.length

        # time step
        dt = ta.step

        #
        # Get eigenenergies and transformation matrix of the Hamiltonian
        #
        if True:
            hD, SS = numpy.linalg.eigh(ham.data)   
               
        #
        #  Find all transition frequencies
        # 
        Om = numpy.zeros((Na, Na))
        for a in range(Na):
            for b in range(Na):
                Om[a,b] = hD[a] - hD[b]
                
        # number of baths - one per monomer            
        Nb = sbi.N

        self.Nt = length
        Nt = self.Nt
        
        #
        # Site K_m operators 
        #

        if not multi_ex:
            Km = numpy.zeros((Nb, Na, Na), dtype=numpy.float64) 
        else:
            Km = numpy.zeros((2*Nb, Na, Na), dtype=numpy.float64)
        # Transform site operators       
        S1 = scipy.linalg.inv(SS)
        #FIXME: SBI should also be basis controlled
        for ns in range(Nb): 
            Km[ns,:,:] = numpy.dot(S1, numpy.dot(sbi.KK[ns,:,:],SS))

        if multi_ex:

            try:
                ii = sbi.system.twoex_state[0,0]
                #print(ii)
                #if ii >= 0:
                #    print("Something is wrong")
                #    raise Exception("twoex_state property ill-defined.")
            except:
                raise Exception("System requires the twoex_state property.")

            # There are also two-exciton bath which can also be converted into one per monomer
            Nb2 = sbi.N

            # projectors to two-exciton states
            K2 = numpy.zeros((Nb2,Na,Na), dtype=REAL)
            i_start = ham.rwa_indices[2]  # here the two-exciton bands starts
            ii = 0
            for ms in range(i_start, ham.dim):
                K2[ii,ms,ms] = 1.0
                ii += 1
            
            # here we do two-exciton projectors
            for ns in range(Nb): 
                for ms in range(Nb):
                    if ms != ns:
                        n2ex = sbi.system.twoex_state[ms,ns] - Nb - 1  # we number the K2 projectors from zero
                        Km[Nb + ns,:,:] += numpy.dot(S1, numpy.dot(K2[n2ex,:,:],SS))

        #
        # \Lambda_m operator
        #

        # Integrals of correlation functions from the set 
        # --- Lm = numpy.zeros((Nt, Nb, Na, Na), dtype=numpy.complex128)
        if not multi_ex:
            Lm = numpy.zeros((Na, Na, Nb, Nt), dtype=numpy.complex128)
        else:
            Lm = numpy.zeros((Na, Na, 2*Nb, Nt), dtype=numpy.complex128)


        #print("FAST INTEGRATION")
        eexp = numpy.exp(-1.0j*Om[:,:, numpy.newaxis]*tm[numpy.newaxis,numpy.newaxis,:])

        for ms in range(Nb):
            #for ns in range(Nb):
                ns = ms
                
                # correlation function of site ns (if ns == ms)
                # or a cross-correlation function of sites ns and ms
                # ms counts baths and it starts from 0, but the excited states are counted from 1
                rc1 = sbi.get_coft(ms+1, ns+1)

                print("MAX:", ms, numpy.max(numpy.abs(rc1)) )
                rc = numpy.einsum("t,ijt->ijt", rc1[0:length], eexp)

                cc_mnab = _integrate_last_axis(rc, dt)
                Lm[:,:,ms,:] = numpy.einsum("ijt,ij->ijt",cc_mnab,Km[ns,:,:])

        if multi_ex:
            for ms in range(Nb):
                ns = ms
                #print("two-ex", ns, ":")
                rc1 = sbi.get_coft(ms+1, ns+1)
                rc = numpy.einsum("t,ijt->ijt", rc1[0:length], eexp)
                cc_mnab = _integrate_last_axis(rc, dt)
                Lm[:,:,Nb + ms,:] = numpy.einsum("ijt,ij->ijt",cc_mnab,Km[Nb + ms,:,:])

        # create the Hermite conjuged version of \Lamnda_m
        if not multi_ex:
            Ld = numpy.zeros((Na, Na, Nb, Nt), dtype=numpy.complex128)
            Nthrough = Nb
        else:
            Ld = numpy.zeros((Na, Na, 2*Nb, Nt), dtype=numpy.complex128)
            Nthrough = 2*Nb

        for tt in range(Nt):
            for ms in range(Nthrough):
                Ld[:, :,ms,tt] += numpy.conj(numpy.transpose(Lm[:,:,ms,tt]))        
            
        if self.as_operators:
            
            # save the operators - propagation methods must know about them
            self.Km = Km
            self.Lm = Lm
            self.Ld = Ld
            
        else:
            
            # save the relaxation tensor
            RR = self._convert_operators_2_tensor(Km, Lm, Ld)
            

            if True:
                self.data = RR
                self._data_initialized = True

        self._is_initialized = True
        self.is_time_dependent = True

    
    def _convert_operators_2_tensor(self, Km, Lm, Ld):
        """Converts operator representation to the tensor one
        
        Convertes operator representation of the Redfield tensor
        into a truely tensor representation
        
        Parameters
        ----------
        
        Km : 3D array
            K_m operators
            
        Lm : 3D array
            \Lambda_m operators
            
        Ld : 3D array
            Hermite conjuget \Lambda_m operators
            
        """    
        
        Na = self.Hamiltonian.data.shape[0]
        Nb = self.SystemBathInteraction.N
        Nt = self.Nt
        
        if self.SystemBathInteraction.system.mult > 1:
            Nb = 2*Nb

        RR = numpy.zeros((Nt, Na, Na, Na, Na), dtype=numpy.complex128)
        
        for m in range(Nb):
            #print("m =", m ,"of", Nb)
            for tt in range(Nt):
                KmLm = numpy.dot(Km[m,:,:],Lm[:,:,m,tt])
                LdKm = numpy.dot(Ld[:,:,m,tt],Km[m,:,:])
                for a in range(Na):
                    for b in range(Na):
                        for c in range(Na):
                            RR[tt,a,b,c,b] -= KmLm[a,c]
                            RR[tt,a,b,a,c] -= LdKm[c,b]
                            for d in range(Na):
                            
                                RR[tt,a,b,c,d] += (Km[m,a,c]*Ld[d,b,m,tt] 
                                                + Lm[a,c,m,tt]*Km[m,d,b])
                                #if b == d:
                                #    RR[tt,a,b,c,d] -= KmLm[a,c] 
                                #if a == c:
                                #    RR[tt,a,b,c,d] -= LdKm[d,b]

        """
        for mi in range(Nb):
            #print("m =", m ,"of", Nb)
            m = mi + Nb
            for tt in range(Nt):
                KmLm = numpy.dot(Km[m,:,:],Lm[:,:,m,tt])
                LdKm = numpy.dot(Ld[:,:,m,tt],Km[m,:,:])
                for a in range(Na):
                    for b in range(Na):
                        for c in range(Na):
                            for d in range(Na):
                            
                                RR[tt,a,b,c,d] += (Km[m,a,c]*Ld[d,b,m,tt] 
                                                + Lm[a,c,m,tt]*Km[m,d,b])
                                if b == d:
                                    RR[tt,a,b,c,d] -= KmLm[a,c] 
                                if a == c:
                                    RR[tt,a,b,c,d] -= LdKm[d,b]
        """

        return RR
        
        
    def transform(self, SS, inv=None):
        """Transformation of the tensor by a given matrix
        
        
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
        dim = SS.shape[0]

        if not self._data_initialized:
            for tt in range(self.Nt):
                for m in range(self.Km.shape[0]):
                    self.Lm[:, :,m,tt] = \
                    numpy.dot(S1,numpy.dot(self.Lm[:, :, m, tt],SS))
                    self.Ld[:, :,m,tt] = \
                    numpy.dot(S1,numpy.dot(self.Ld[:, :, m, tt],SS))            
            for m in range(self.Km.shape[0]):
                self.Km[m, :, :] = numpy.dot(S1,numpy.dot(self.Km[m, :, :],SS))
                
            return
        
        
        if (self.manager.warn_about_basis_change):
                print("\nQr >>> Relaxation tensor '%s' changes basis" %self.name)
           
        
        for tt in range(self.Nt):
            for c in range(dim):
                for d in range(dim):
                    self._data[tt,:,:,c,d] = \
                        numpy.dot(S1,numpy.dot(self._data[tt,:,:,c,d],SS))
                
            for a in range(dim):
                for b in range(dim):
                    self._data[tt,a,b,:,:] = \
                        numpy.dot(S1,numpy.dot(self._data[tt,a,b,:,:],SS))

            
    def secularize(self):
        """Secularizes the relaxation tensor


        """
        if self.as_operators:
            raise Exception("Cannot be secularized in an opeator form")
            
        else:
            N = self.data.shape[1]
            for ii in range(N):
                for jj in range(N):
                    for kk in range(N):
                        for ll in range(N):
                            if not (((ii == jj) and (kk == ll)) 
                                or ((ii == kk) and (jj == ll))) :
                                    self.data[:,ii,jj,kk,ll] = 0
                                    


def _integrate(f, dt):
    """Cummulative simpson rule for integration (even number of points also handled)
    
    """
    N = len(f)
    if N < 2:
        raise ValueError("Need at least two points for integration")
    
    F = numpy.zeros_like(f)
    
    # Integrate in blocks of 2 intervals (3 points) using Simpson's rule
    last = N - 1 if N % 2 == 1 else N - 2  # leave out last interval if N is even
    
    for i in range(2, last + 1, 2):
        F[i] = F[i - 2] + (dt / 3) * (f[i - 2] + 4 * f[i - 1] + f[i])
        F[i - 1] = (F[i - 2] + F[i]) / 2  # midpoint estimate

    # Handle last interval with trapezoidal rule if needed
    if N % 2 == 0:
        F[-1] = F[-2] + 0.5 * dt * (f[-2] + f[-1])
    
    return F


def _integrate_last_axis(f, dt):
    """
    Cumulative Simpson integration over the last axis of a 2D or higher NumPy array.

    Parameters:
        f : ndarray
            Input array with shape (..., N), where integration is over the last axis (time).
        dt : float
            Uniform time step between samples.

    Returns:
        F : ndarray
            Cumulative integral of f along the last axis, with the same shape as f.
    """
    if f.ndim < 2:
        raise ValueError("Input array must be at least 2D (e.g. shape (i, t))")

    N = f.shape[-1]
    if N < 2:
        raise ValueError("Need at least two points along the last axis")

    F = numpy.zeros_like(f)

    last = N - 1 if N % 2 == 1 else N - 2  # leave last interval for trapezoid if needed

    idxs = numpy.arange(2, last + 1, 2)
    for i3 in idxs:
        F[..., i3] = F[..., i3 - 2] + (dt / 3) * (
            f[..., i3 - 2] + 4 * f[..., i3 - 1] + f[..., i3]
        )
        F[..., i3 - 1] = (F[..., i3 - 2] + F[..., i3]) / 2  # midpoint

    if N % 2 == 0:
        # Final trapezoidal step
        F[..., -1] = F[..., -2] + 0.5 * dt * (f[..., -2] + f[..., -1])

    return F