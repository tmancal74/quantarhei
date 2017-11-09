# -*- coding: utf-8 -*-
"""
    Quantarhei package (http://www.github.com/quantarhei)

    redfield tensor module
    
*******************************************************************************

    REDFIELD RELAXATION TENSOR

*******************************************************************************
"""
import numpy
import scipy

from .systembathinteraction import SystemBathInteraction
from ..hilbertspace.hamiltonian import Hamiltonian

from .relaxationtensor import RelaxationTensor
from ...core.managers import  energy_units
from ...core.parallel import block_distributed_range
from ...core.parallel import start_parallel_region, close_parallel_region
from ...core.parallel import distributed_configuration

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
    
    
    def __init__(self, ham, sbi, initialize=True,
                 cutoff_time=None, as_operators=False,
                 name=""):
                     
        super().__init__()
        
        #
        # Check the types of the arguments
        #
    
        # Hamiltonian
        if not isinstance(ham, Hamiltonian):
            raise Exception("First argument must be a Hamiltonian")
            
        # SystemBathInteraction
        if not isinstance(sbi, SystemBathInteraction):
            raise Exception("Second argument must be of" +
                            " the SystemBathInteraction type")

        self.Hamiltonian = ham
        self.SystemBathInteraction = sbi
        
        self.dim = self.Hamiltonian.dim
        self.name = name
        
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
        #print("Reference Redfield implementation ...")
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

        Km = numpy.zeros((Nb, Na, Na), dtype=numpy.float64) 
        # Transform site operators       
        S1 = scipy.linalg.inv(SS)
        #FIXME: SBI should also be basis controlled
        for ns in range(Nb): 
            Km[ns,:,:] = numpy.dot(S1, numpy.dot(sbi.KK[ns,:,:],SS))
        
        #
        # \Lambda_m operator
        #
        
        # Integrals of correlation functions from the set      
        Lm = numpy.zeros((Nb, Na, Na), dtype=numpy.complex128)
        
        #######################################################################
        # PARALLELIZATION
        #######################################################################

        start_parallel_region()
        for ms in block_distributed_range(0,Nb): #range(Nb):
            #for ns in range(Nb):
            if not multi_ex:
                ns = ms
                
                # correlation function of site ns (if ns == ms)
                # or a cross-correlation function of sites ns and ms
                
                #FIXME: reaching correct correlation function is a nightmare!!!
                rc1 = sbi.CC.get_coft(ms, ns) 
                
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
                        Lm[ms,a,b] += cc_mnab*Km[ns,a,b] 
             
        # perform reduction of Lm
        distributed_configuration().allreduce(Lm, operation="sum")
        close_parallel_region()

        #######################################################################
        #  END PARALLELIZATION
        #######################################################################
                        
        # create the Hermite conjuged version of \Lamnda_m
        Ld = numpy.zeros((Nb, Na, Na), dtype=numpy.complex128)
        for ms in range(Nb):
            Ld[ms, :, :] += numpy.conj(numpy.transpose(Lm[ms,:,:]))        
            
        self._post_implementation(Km, Lm, Ld)
        
        
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
        
        RR = numpy.zeros((Na, Na, Na, Na), dtype=numpy.complex128)
        
        #######################################################################
        # PARALLELIZATION
        #######################################################################

        start_parallel_region()
        for m in block_distributed_range(0,Nb): #range(Nb):
            Kd = numpy.transpose(Km[m,:,:])
            KdLm = numpy.dot(Kd,Lm[m,:,:])
            LdKm = numpy.dot(Ld[m,:,:],Km[m,:,:])
            for a in range(Na):
                for b in range(Na):
                    for c in range(Na):
                        for d in range(Na):
                            
                            RR[a,b,c,d] += (Km[m,a,c]*Ld[m,d,b] 
                                            + Lm[m,a,c]*Kd[d,b])
                            if b == d:
                                RR[a,b,c,d] -= KdLm[a,c] 
                            if a == c:
                                RR[a,b,c,d] -= LdKm[d,b]
                                
        # perform reduction of the RR
        distributed_configuration().allreduce(RR, operation="sum")

        close_parallel_region()
        #######################################################################
        # END PARALLELIZATION
        #######################################################################
        
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
            
            self.data = self._convert_operators_2_tensor(self.Km,
                                                         self.Lm, self.Ld)
                                                         
            self.as_operators = False
         
                      