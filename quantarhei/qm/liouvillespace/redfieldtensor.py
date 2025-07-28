# -*- coding: utf-8 -*-
"""
    Redfield Tensor
    
    
    Class Details
    -------------
    
"""
#import time
import numpy
import scipy

from .systembathinteraction import SystemBathInteraction
from ..hilbertspace.hamiltonian import Hamiltonian

from .relaxationtensor import RelaxationTensor
from ...core.managers import  energy_units
from ...core.parallel import block_distributed_range
from ...core.parallel import start_parallel_region, close_parallel_region
from ...core.parallel import distributed_configuration
#from ...core.managers import BasisManaged
from ...utils.types import BasisManagedComplexArray

from ... import REAL, COMPLEX

import quantarhei as qr

class RedfieldRelaxationTensor(RelaxationTensor):
    """Redfield Relaxation Tensor
        
        
    Parameters
    ----------
    
    ham : Hamiltonian
        Hamiltonian of the system.
        
    sbi : SystemBathInteraction
        Object specifying the system-bath interaction
        
    initialize : bool
        If True, the tensor will be imediately calculated
        
    cutoff_time : float
        Time in femtoseconds after which the integral kernel in the 
        definition of the relaxation tensor is assummed to be zero.
            
    as_operators : bool
        If True the tensor will not be constructed. Instead a set of
        operators whose application is equal to the application of the
        tensor will be defined and stored
            
    Methods
    -------
        
    secularize()
        Deletes non-secular terms in the relaxation tensor. If the tensor
        is represented by operators (i.e. the as_operators atribute is
        True), the methods raises an Exception. Use `convert_2_tensor`
        first.
        
    initialize()
        Initializes the Redfield tensor if created with parameter
        `initialize=False`
            
    convert_2_tensor()
        When created with `as_operators=True` the tensor is internally
        represented by a set of operators. This method converts this
        representation to a tensor
            
            
    """

    #data = BasisManagedComplexArray("data")
    Km = BasisManagedComplexArray("Km")
    Lm = BasisManagedComplexArray("Lm")
    Ld = BasisManagedComplexArray("Ld")    

    def __init__(self, ham, sbi, initialize=True,
                 cutoff_time=None, as_operators=False,
                 name=""):
                     
        self._initialize_basis()
        
        #
        # Check the types of the arguments
        #
    
        # Hamiltonian
        if not isinstance(ham, Hamiltonian):
            raise Exception("First argument must be a Hamiltonian")
            
        # SystemBathInteraction
        if not isinstance(sbi, SystemBathInteraction):
            if sbi is not None:
                raise Exception("Second argument must be of" +
                            " the SystemBathInteraction type")

        self.Hamiltonian = ham
        self.SystemBathInteraction = sbi
        
        self.dim = self.Hamiltonian.dim
        self.name = name
        
        self._data_initialized = False
        self._is_initialized = False
        self._has_cutoff_time = False
        self.as_operators = as_operators
        
        if not self.as_operators:
            self.data = numpy.zeros((self.dim,self.dim,self.dim,self.dim),
                                    dtype=numpy.complex128)
        # set cut-off time
        if cutoff_time is not None:
            self.cutoff_time = cutoff_time
            self._has_cutoff_time = True            
            
        # initialize the tensor right at creation
        if initialize:  
            
            try: 
                from ...implementations.qm.liouvillespace.redfieldtensor \
                        import redfieldtensor
            except:
                pass
            
            
            with energy_units("int"):
                self._implementation(ham, sbi)   
                

        self.Iterm = None
        self.has_Iterm = False


    def apply(self, oper, copy=True):
        """Applies the relaxation tensor on a superoperator
        
        """
        
        if self.as_operators:
            
            print("Applying Relaxation tensor")
            if copy:
                import copy
                oper_ven = copy.copy(oper)
            else:
                oper_ven = oper
            
            rho1 = oper.data
            # apply tensor to data
            
            Lm = self.Lm            
            Km = self.Km
            Ld = self.Ld
            
            Kd = numpy.zeros(Km.shape, dtype=numpy.float64)
            Nm = Km.shape[0]
            ven = numpy.zeros(oper.data.shape, dtype=numpy.complex128)
            for mm in range(Nm):
                Kd[mm, :, :] = numpy.transpose(Km[mm, :, :])
            
                ven += (
                numpy.dot(Km[mm,:,:],numpy.dot(rho1, Ld[:,:,mm]))
                +numpy.dot(Lm[:,:,mm],numpy.dot(rho1, Kd[mm,:,:]))
                -numpy.dot(numpy.dot(Kd[mm,:,:],Lm[:,:,mm]), rho1)
                -numpy.dot(rho1, numpy.dot(Ld[:,:,mm],Km[mm,:,:])))
                
            oper_ven.data = ven
                
            return oper_ven
            
        else:
            
            return super().apply(oper, copy=copy)


    def get_population_rate(self, N, M):
        """Returns the relaxation rate between states N -> M
        
        
        """

        if self.as_operators:

            rt = 0.0
            for ii in range(self.Lm.shape[2]):
                rt += 2.0*numpy.real(self.Lm[N,M,ii]*self.Km[ii,N,M])
            
            return rt
        
        else:
            return super().get_population_rate(N, M)
        
        
    def get_dephasing_rate(self, N, M):
        """Returns the dephasing rate of a coherence between states N and M
        
        
        """
        if self.as_operators:

            rt = 0.0
            for ii in range(self.Lm.shape[2]):

                A_a = numpy.dot(self.Km[ii,N,:],self.Lm[:,N,ii])
                A_b = numpy.dot(self.Ld[M,:,ii],self.Km[ii,:,M])
                KL = self.Km[ii,N,N]*self.Ld[M,M,ii]
                LK = self.Lm[N,N,ii]*self.Km[ii,M,M]
                rt += A_a + A_b - KL - LK
            
            return rt 
                   
        else:  
            return super().get_dephasing_rate(N, M)


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
        
        if self.as_operators:
            
            if (self.manager.warn_about_basis_change):
                print("\nQr >>> Operators of relaxation"+
                      " in Redfield tensor '%s' changes basis" %self.name)
        
            if inv is None:
                S1 = numpy.linalg.inv(SS)
            else:
                S1 = inv

            _Lm = numpy.zeros_like(self._Lm, dtype=COMPLEX)
            _Ld = numpy.zeros_like(self._Ld, dtype=COMPLEX)
            _Km = numpy.zeros_like(self._Km, dtype=COMPLEX)
            for m in range(self._Km.shape[0]):
                _Lm[:,:,m] = numpy.dot(S1,numpy.dot(self._Lm[:,:,m], SS))
                _Ld[:,:,m] = numpy.dot(S1,numpy.dot(self._Ld[:,:,m], SS))
                _Km[m,:,:] = numpy.dot(S1,numpy.dot(self._Km[m,:,:], SS))
     
            self._Lm = _Lm
            self._Ld = _Ld
            self._Km = _Km
            
        else:
    
            super().transform(SS)
            

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
        
        qr.log_detail("Reference time-independent Redfield tensor calculation")
        #print("Reference Redfield implementation ...")
        #
        # dimension of the Hamiltonian (includes excitons
        # with all multiplicities specified at its creation)
        #
        Na = ham.dim #data.shape[0]
        
        # time axis
        ta = sbi.TimeAxis
        dt = ta.step
        
        #
        # is this beyond single excitation band?
        #
        multi_ex = False
        
        # figure out if the aggregate or molecule specifies more than
        # one exciton band
        if sbi.aggregate is not None:
            agg = sbi.aggregate
            if agg.mult > 1:
                multi_ex = True 
        elif sbi.molecule is not None:
            mol = sbi.molecule
            if mol.mult > 1:
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

        #
        # Get eigenenergies and transformation matrix of the Hamiltonian
        #
        if True:
            # FIXME: here we need to access ham._data (we want to protect basis)
            #
            # THIS ASSUMES WE ARE IN SITE BASIS
            # FIXME: devise a mechanism to ensure this!!!!
            #
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

        #
        # Site K_m operators 
        #
        Nex = Nb
        if multi_ex:
            Nex = 2*Nb 

        Km = numpy.zeros((Nex, Na, Na), dtype=numpy.float64) 
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
        if not multi_ex:     
            Lm = numpy.zeros((Na, Na, Nb), dtype=numpy.complex128)
        else:
            Lm = numpy.zeros((Na, Na, 2*Nb), dtype=numpy.complex128)
            Nex = 2*Nb

        #######################################################################
        # PARALLELIZATION
        #######################################################################

        start_parallel_region()
        for ms in block_distributed_range(0, Nb): #range(Nb):
            qr.log_quick("Calculating bath component", ms, "of", Nb, end="\r")
            #print(ms, "of", Nb)
            #for ns in range(Nb):

            ns = ms
            
            # correlation function of site ns (if ns == ms)
            # or a cross-correlation function of sites ns and ms
            
            #FIXME: reaching correct correlation function is a nightmare!!!
            rc1 = sbi.CC.get_coft(ms, ns)  
                    
            self._guts_Cmplx_Splines(ms, Lm, Km, Na, Om, length, rc1, tm)
             
        # perform reduction of Lm
        qr.log_quick()
        distributed_configuration().allreduce(Lm, operation="sum")
        close_parallel_region()

        #######################################################################
        #  END PARALLELIZATION
        #######################################################################
        eexp = numpy.exp(-1.0j*Om[:,:, numpy.newaxis]*tm[numpy.newaxis,numpy.newaxis,:])

        if multi_ex:
            for ms in range(Nb):
                ns = ms
                #print("two-ex", ns, ":")
                rc1 = sbi.get_coft(ms+1, ns+1)
                rc = numpy.einsum("t,ijt->ijt", rc1[0:length], eexp)
                cc_mnab = _integrate_last_axis(rc, dt)
                Lm[:,:,Nb + ms] = numpy.einsum("ij,ij->ij",cc_mnab[:,:,-1],Km[Nb + ms,:,:])

        # create the Hermite conjuged version of \Lamnda_m
        if not multi_ex:
            Ld = numpy.zeros((Na, Na, Nb), dtype=numpy.complex128)
            Nex = Nb
        else:
            Ld = numpy.zeros((Na, Na, 2*Nb), dtype=numpy.complex128)
            Nex = 2*Nb

        for ms in range(Nex):
            Ld[:, :, ms] += numpy.conj(numpy.transpose(Lm[:,:,ms]))        
            
        self._post_implementation(Km, Lm, Ld)
        
        qr.log_detail("... Redfield done")
        
    def _post_implementation(self, Km, Lm, Ld):
        """When components of the tensor are calculated, should they be
        saved or converted into full tensor?
        
        """
            
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


    def _guts_Cmplx_Splines(self, ms, Lm, Km, Na, Om, length, rc1, tm):
        
        for a in range(Na):
            for b in range(Na):
                
                # argument of the integration
                eexp = numpy.exp(-1.0j*Om[a,b]*tm) 
                rc = rc1[0:length]*eexp
                
                # spline integration instead of FFT
                rr = numpy.real(rc)
                ri = numpy.imag(rc)
                sr = scipy.interpolate.UnivariateSpline(tm,
                            rr, s=0).antiderivative()(tm)
                si = scipy.interpolate.UnivariateSpline(tm,
                            ri, s=0).antiderivative()(tm)
                        
                # we take the last value (integral to infinity)
                cc_mnab = (sr[length-1] + 1.0j*si[length-1]) 

                # \Lambda_m operators
                Lm[a,b,ms] += cc_mnab*Km[ms,a,b] 
      
        
    def _guts_Cmplx_Sum(self, ms, Lm, Km, Na, Om, length, rc1, tm):
        
        import scipy

        dt = tm[1]-tm[0]
        for a in range(Na):
            for b in range(Na):
                
                # argument of the integration
                eexp = numpy.exp(-1.0j*Om[a,b]*tm) 
                rc = rc1[0:length]*eexp
            
                # we take integral to infinity
                #cc_mnab = numpy.sum(rc)*dt
                cc_mnab = scipy.integrate.trapz(rc, dx=dt)

                # \Lambda_m operators
                Lm[a,b,ms] += cc_mnab*Km[ms,a,b] 
                
                
    def _guts_Cmplx_FFT(self, ms, Lm, Km, Na, Om, length, rc1, tm):

        for a in range(Na):
            for b in range(Na):
                
                # argument of the integration
                eexp = numpy.exp(-1.0j*Om[a,b]*tm) 
                rc = rc1[0:length]*eexp
                
                # spline integration instead of FFT
                rr = numpy.real(rc)
                ri = numpy.imag(rc)
                sr = scipy.interpolate.UnivariateSpline(tm,
                            rr, s=0).antiderivative()(tm)
                si = scipy.interpolate.UnivariateSpline(tm,
                            ri, s=0).antiderivative()(tm)
                        
                # we take the last value (integral to infinity)
                cc_mnab = (sr[length-1] + 1.0j*si[length-1]) 

                # \Lambda_m operators
                Lm[a,b,ms] += cc_mnab*Km[ms,a,b]   
                
            
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
        Nb = Lm.shape[2] #self.SystemBathInteraction.N
        
        RR = numpy.zeros((Na, Na, Na, Na), dtype=numpy.complex128)
        
        #######################################################################
        # PARALLELIZATION
        #######################################################################
        
        #from ...implementations.cython.loopit import loopit
        

        start_parallel_region()
        #tt1 = time.time()

        for m in block_distributed_range(0,Nb): #range(Nb):
            
            Kd = numpy.transpose(Km[m,:,:])
#            KdLm = numpy.dot(Kd,Lm[m,:,:])
#            LdKm = numpy.dot(Ld[m,:,:],Km[m,:,:])
#            for a in range(Na):
#                print(a)
#                for b in range(Na):
#                    for c in range(Na):
#                        for d in range(Na):
#                            
#                            RR[a,b,c,d] += (Km[m,a,c]*Ld[m,d,b] 
#                                            + Lm[m,a,c]*Kd[d,b])
#                            if b == d:
#                                RR[a,b,c,d] -= KdLm[a,c] 
#                            if a == c:
#                                RR[a,b,c,d] -= LdKm[d,b]
                                
            _loopit(Km, Kd, Lm, Ld, Na, RR, m)
        
        # perform reduction of the RR
        distributed_configuration().allreduce(RR, operation="sum")
        #tt2 = time.time()
        
        close_parallel_region()
        #######################################################################
        # END PARALLELIZATION
        #######################################################################
        
        #print(tt2-tt1)
        return RR


    def initialize(self):
        """Initializes the Redfield tensor with values 
        
        """
        self._implementation(self.Hamiltonian,
                             self.SystemBathInteraction)


    def convert_2_tensor(self):
        """Converts internally the operator representation to a tensor one
        
        Converst the representation of the relaxation tensor through a set
        of operators into a tensor representation.
        
        """
        
        if self.as_operators:
            
            RR = self._convert_operators_2_tensor(self.Km, self.Lm, self.Ld)
            if True:
                self.data = RR
                self._data_initialized = True
                                                         
            self.as_operators = False
         
 

def _loopit(Km, Kd, Lm, Ld, Na, RR, m):

    
    KdLm = numpy.dot(Kd,Lm[:,:,m])
    LdKm = numpy.dot(Ld[:,:,m],Km[m,:,:])
    for a in range(Na):
        for b in range(Na):
            for c in range(Na):
                for d in range(Na):
                    
                    RR[a,b,c,d] += (Km[m,a,c]*Ld[d,b,m] 
                                    + Lm[a,c,m]*Kd[d,b])
                    if b == d:
                        RR[a,b,c,d] -= KdLm[a,c] 
                    if a == c:
                        RR[a,b,c,d] -= LdKm[d,b]
            

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