# -*- coding: utf-8 -*-
import numpy
import scipy

from .redfieldtensor import RedfieldRelaxationTensor
from ...core.time import TimeDependent

class TDRedfieldRelaxationTensor(RedfieldRelaxationTensor, TimeDependent):
    
    
    
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
        Lm = numpy.zeros((Nt, Nb, Na, Na), dtype=numpy.complex128)
        for ms in range(Nb):
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
                        # #### cc_mnab = (sr[length-1] + 1.0j*si[length-1]) 
                        cc_mnab = (sr + 1.0j*si)
                        
                        # \Lambda_m operators
                        Lm[:,ms,a,b] += cc_mnab*Km[ns,a,b] 
             
        
        # create the Hermite conjuged version of \Lamnda_m
        Ld = numpy.zeros((Nt, Nb, Na, Na), dtype=numpy.complex128)
        for tt in range(Nt):
            for ms in range(Nb):
                Ld[tt, ms, :, :] += numpy.conj(numpy.transpose(Lm[tt,ms,:,:]))        
            
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
        Nt = self.Nt
        
        RR = numpy.zeros((Nt, Na, Na, Na, Na), dtype=numpy.complex128)
        
        for m in range(Nb):
            #print("m =", m ,"of", Nb)
            for tt in range(Nt):
                KmLm = numpy.dot(Km[m,:,:],Lm[tt,m,:,:])
                LdKm = numpy.dot(Ld[tt,m,:,:],Km[m,:,:])
                for a in range(Na):
                    for b in range(Na):
                        for c in range(Na):
                            for d in range(Na):
                            
                                RR[tt,a,b,c,d] += (Km[m,a,c]*Ld[tt,m,d,b] 
                                                + Lm[tt,m,a,c]*Km[m,d,b])
                                if b == d:
                                    RR[tt,a,b,c,d] -= KmLm[a,c] 
                                if a == c:
                                    RR[tt,a,b,c,d] -= LdKm[d,b]
        
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
                    self.Lm[tt, m, :, :] = \
                    numpy.dot(S1,numpy.dot(self.Lm[tt, m, :, :],SS))
                    self.Ld[tt, m, :, :] = \
                    numpy.dot(S1,numpy.dot(self.Ld[tt, m, :, :],SS))            
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