# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot as plt
import scipy.interpolate as interp

from ..corfunctions.correlationfunctions import c2g
from ...core.managers import energy_units
from .tdfoerstertensor import TDFoersterRelaxationTensor
from .tdfoerstertensor import _td_reference_implementation
from ... import COMPLEX, REAL


class NEFoersterRelaxationTensor(TDFoersterRelaxationTensor):
    """Weak resonance coupling relaxation tensor by Foerster theory
    
    Non-equilibrium version according to 
    
    J. Seibt and T. Mančal, J. Chem. Phys. 146 (2017) 174109
    
    
    """
    def __init__(self, ham, sbi, as_kernel=False, 
                 initialize=True, cutoff_time=None):
        """Initiation is the same as for the TDFoerster tensor
        
        """
        self.as_kernel = as_kernel
        super().__init__(ham, sbi, initialize, cutoff_time)
        
  
    def initialize(self):
        
        tt = self.SystemBathInteraction.TimeAxis.data
        Nt = len(tt)
        
        self.nsc = True
        
        
        #
        # Tensor data
        #
        Na = self.dim
        self.data = numpy.zeros((Nt,Na,Na,Na,Na),dtype=COMPLEX)
        
        with energy_units("int"):

            # Hamiltonian matrix
            HH = self.Hamiltonian.data

            sbi = self.SystemBathInteraction
            Na = self.dim
            
            # line shape functions
            gt = numpy.zeros((Na, sbi.TimeAxis.length),
                             dtype=numpy.complex64)
    
            # SBI is defined with "sites"
            for ii in range(1, Na):
                gt[ii,:] = c2g(sbi.TimeAxis, sbi.CC.get_coft(ii-1,ii-1))
                

            
            # reorganization energies
            ll = numpy.zeros(Na)
            for ii in range(1, Na):
                ll[ii] = sbi.CC.get_reorganization_energy(ii-1,ii-1)


            if self.nsc:
                self.data =  \
                    self.td_reference_implementation(Na, Nt, HH, tt, gt, ll)

            else:
                KK = self.td_reference_implementation(Na, Nt, HH, tt, gt, ll)
                for a in range(Na):
                    for b in range(Na):
                        if a != b:
                            self.data[:,a,a,b,b] = KK[:,a,b]
                            
                self.updateStructure()
                self.add_dephasing()
    


    def td_reference_implementation(self, Na, Nt, HH, tt, gt, ll):
        """ Overloaded implementation method replacing integration kernel
        
        """
        
        if self.nsc:
            self._initial_term_pre_nsc(Na, Nt, HH, tt, gt, ll)
            #
            # Here we calculate the time-dependent non-secular tensor
            #
            if not self.as_kernel:
                return _nsc_reference_implementation(Na, Nt, HH, tt,
                                                gt, ll, _nsc_fintegral)
            #
            # Here we prepare the kernel so that it can be returned
            # on demand for any value of propagation time t
            #
            else:
                return _nsc_reference_implementation(Na, Nt, HH, tt,
                                                gt, ll, _nsc_fkernel)
        else:
            self._initial_term_pre(Na, Nt, HH, tt, gt, ll)
            return _td_reference_implementation(Na, Nt, HH, tt,
                                            gt, ll, _ne_fintegral)
    

    def get_kernel_tensor(self, ti):
        """Return the relaxation tensor kernel at a time given by the index ti
        
        """
        pass



    def _initial_term_pre(self, Na, Nt, HH, tt, gt, ll):
        """ Quantities to easily calculate initial term with known initial
            condition
            
        """
        II = numpy.zeros((Nt, Na, Na), dtype=COMPLEX)
        
        JJ = HH.data
        
        for aa in range(Na):
            for mm in range(Na):
                if mm != aa:
                    II[:,aa,mm] = 2.0*JJ[aa,mm]*numpy.imag( \
                                  numpy.exp(-1j*(JJ[mm,mm]-JJ[aa,aa])*tt) \
                                  *numpy.exp(-numpy.conj(gt[aa,:])-gt[mm,:])) 
        
        self.II = II
        self.has_Iterm = True



    def initial_term(self, rhoi):
        """ Inhomogeneous (initial condition) term 
            of the non-equilibrium Foerster theory
            
        """
        
        if self.nsc:
            self.initial_term_nsc(rhoi)
            return 
        
        II = numpy.zeros(self.II.shape, dtype=COMPLEX)
        Na = self.II.shape[1]
        
        for aa in range(Na):
            for mm in range(Na):
                II[:,aa,aa] += self.II[:,aa,mm]*rhoi.data[mm,aa]
        
        self.Iterm = II
        self.has_Iterm = True
  
        
    def _initial_term_pre_nsc(self, Na, Nt, HH, tt, gt, ll):
        """ Quantities to easily calculate initial term in non-secular theory
            with known initial condition
            
        """
        II = numpy.zeros((Nt, Na, Na), dtype=COMPLEX)
        
        JJ = HH.data
        
        for aa in range(Na):
            for bb in range(Na):
                bb != aa
                II[:,aa,bb] += \
                    JJ[aa,bb]*numpy.exp(-1j*(JJ[bb,bb]-JJ[aa,aa])*tt) \
                             *numpy.exp(-numpy.conj(gt[aa,:])-gt[bb,:]) 
        
        self.II = II
        self.has_Iterm = True 

        
        
    def initial_term_nsc(self, rhoi): 
        """ Inhomogeneous (initial condition) term 
            of the effective non-secular non-equilibrium Foerster theory
            
        """        
        II = numpy.zeros(self.II.shape, dtype=COMPLEX)
        Na = self.II.shape[1]
        
        for aa in range(Na):
            for bb in range(Na):
                for cc in range(Na):
                    II[:,aa,bb] += -1j*self.II[:,aa,cc]*rhoi.data[cc,bb] \
                                   +1j*rhoi.data[aa,cc]*self.II[:,cc,bb]
        
        self.Iterm = II
        self.has_Iterm = True
        
        

def _kernel_at_t(ti, tt, gtd, gta, ed, ea, ld):
    """ Two-time kernel to be integrated 
    
    
    """ 
    Nt = tt.shape[0]
    gtd_i = gtd[0:ti+1]
    gtd_m = numpy.zeros(Nt, dtype=COMPLEX)
    gtd_m[0:ti+1] = numpy.flip(gtd_i) 
    
    prod = numpy.exp(-gtd-gta +1j*(ed-ea)*tt) \
        *numpy.exp(-2.0*1j*numpy.imag(gtd_m)) \
        *numpy.exp(2.0*1j*numpy.imag(gtd[ti]))
        
    return prod


def _nsc_kernel_at_t(ti, tt, aa, bb, cc, dd, HH, gt):
    """ Two-time kernel to be integrated 
    
    
    """ 
    exp = numpy.exp
    conj = numpy.conj
    Nt = tt.shape[0]

    # expressions for t-tau
    gtb_i = gt[bb, 0:ti+1]
    gta_i = gt[aa, 0:ti+1]
    gtb_m = numpy.zeros(Nt, dtype=COMPLEX)
    gta_m = numpy.zeros(Nt, dtype=COMPLEX)
    gtb_m[0:ti+1] = numpy.flip(gtb_i) 
    gta_m[0:ti+1] = numpy.flip(gta_i)    

  
    ea = HH[aa,aa]
    eb = HH[bb,bb]
    ec = HH[cc,cc]
    ed = HH[dd,dd]
    
    if False: #(bb==cc) and (aa==dd):
        # manually simplified expression for population rates
    
        prod = exp(2.0*1j*numpy.imag(gt[bb,ti])  \
                  -2.0*1j*numpy.imag(gtb_m)  \
                  - gt[bb,:] - gt[aa,:] + 1j*(ea-eb)*tt)
        
    else:
        # general expression for all indices
        
        dl = numpy.eye(HH.shape[0], dtype=REAL)
        
        tt_i = tt[ti] 
        prod = exp( - conj(gt[aa, ti] + gt[cc, :])                  \
               - gt[bb, ti] - gt[dd, :]                             \
               + dl[aa,bb]*(+ conj(gt[aa, ti]) + gt[aa, ti])        \
               + dl[aa,cc]*(- conj(gt[aa, :]) + gta_m - gt[aa, ti]) \
               + dl[aa,dd]*(+ conj(gt[aa, :]) + gt[aa, ti] - gta_m) \
               + dl[bb,cc]*(+ conj(gt[bb, :]) + gt[bb, ti] - gtb_m) \
               + dl[bb,dd]*(- conj(gt[bb, :]) - gt[bb, ti] + gtb_m) \
               + dl[cc,dd]*(conj(gt[cc, :]) + gt[cc, :])            \
               + 1j*((ea-eb)*tt_i)+1j*((ec-ed))*tt)        
        
    return prod


def _integrate_kernel(tt, fce):
    """ Spline integration of a complex function
    
    """
    preal = numpy.real(fce)
    pimag = numpy.imag(fce)
    splr = interp.UnivariateSpline(tt,
                           preal, s=0).antiderivative()(tt)
    spli = interp.UnivariateSpline(tt,
                           pimag, s=0).antiderivative()(tt)
    inte = splr + 1j*spli
    return inte        


def _integrate_kernel_to_t(ti, tt, fce, margin=10):
    """ Spline partial integration of a complex function

    A function of variables tau is integrated from zero to t.

                 t
                /
                |
         f(t) = | d tau fce(tau)
                |
                /
                0

    """
    import scipy.interpolate as interp
    import numpy

    ti_min = margin
    ti_eff = max(ti_min, ti) + 1
    fce_ti = fce[0:ti_eff]
    tt_ti = tt[0:ti_eff]

    preal = numpy.real(fce_ti)
    pimag = numpy.imag(fce_ti)
    splr = interp.UnivariateSpline(tt_ti,
                           preal, s=0).antiderivative()(tt_ti)
    spli = interp.UnivariateSpline(tt_ti,
                           pimag, s=0).antiderivative()(tt_ti)
    inte = splr + 1j*spli

    return inte


def _integrate_three_to_t(ti, tt, fce_t, fce_tau, fce_ttau):
    """ Integrate the product of three functions

    The product of functions of variables t, tau and t-tau, respectively,
    is integrated over the variable tau from zero to t.

                 t
                /
                |
         f(t) = | d tau  fce_t(t)*fce_tau(tau)*fce_ttau(t-tau)
                |
                /
                0

    """
    import numpy

    # convert the function of t-tau into a function of tau
    fce_flip = numpy.zeros(fce_ttau.shape, dtype=COMPLEX)
    fce_flip[0:ti+1] = numpy.flip(fce_ttau[0:ti+1])

    # construct the function to be integrated
    fce = fce_t[ti]*fce_tau*fce_flip

    # integrate
    return _integrate_kernel_to_t(ti, tt, fce)[ti]



def _ne_fintegral(tt, gtd, gta, ed, ea, ld):
    """Time dependent non-equilibrium Foerster integral
    
    
    Parameters
    ----------
    tt : numpy array
        Time 
        
    gtd : numpy array
        lineshape function of the donor transition

    gta : numpy array
        lineshape function of the acceptor transition 
        
    ed : float
        Energy of the donor transition
        
    ea : float
        Energy of the acceptor transition

    ld : float
        Reorganization energy of the donor             

    Returns
    -------
    ret : float
        The value of the Foerster integral            
    
    """

    Nt = tt.shape[0]
    hoft = numpy.zeros(Nt, dtype=COMPLEX)
    
   
    for ti in range(Nt):
        
        #
        # Here we calculate two-time integration kernel 
        #
        prod = _kernel_at_t(ti, tt, gtd, gta, ed, ea, ld)
        

        #
        # the kernel is integrated by splines
        #
        #inte = _integrate_kernel(tt, prod)
        inte = _integrate_kernel_to_t(ti, tt, prod)

        hoft[ti] = inte[ti]

    ret = 2.0*numpy.real(hoft)
    
    
    return ret


def _nsc_reference_implementation(Na, Nt, HH, tt, gt, ll, fce):
                                     
    #
    # Rates between states a and b
    # 
    KK = numpy.zeros((Nt,Na,Na,Na,Na), dtype=COMPLEX)
    fKK = numpy.zeros((Nt,Na,Na,Na,Na), dtype=COMPLEX)
    #fR = numpy.zeros((Nt,Na,Na,Na), dtype=COMPLEX)
    RR = numpy.zeros((Nt,Na,Na), dtype=COMPLEX)
    
    JJ = numpy.zeros(HH.shape, dtype=REAL)
    JJ[:,:] = HH[:,:]
    for ii in range(Na):
        JJ[ii,ii] = 0.0
    
    #
    # Integrals
    #
    for a in range(Na):
        for b in range(Na):
            for c in range(Na):
                for d in range(Na):
                    fKK[:,a,b,c,d] = fce(tt, d, b, a, c, HH, gt)

    #
    # Operator part of the non-eq. Foerster
    #
    for a in range(Na):
        for b in range(Na):
            for c in range(Na):
                RR[:,a,b] -= JJ[a,c]*JJ[c,b]*fKK[:,c,c,b,a] 
               
    #
    # Tensor elements
    #
    for a in range(Na):
        for b in range(Na):
            #b = a
            for c in range(Na):
                for d in range(Na):
                    #d = c
                    KK[:,a,b,c,d] += JJ[a,c]*JJ[d,b]*fKK[:, a, b, c, d] \
                        + numpy.conj(JJ[b,d]*JJ[c,a]*fKK[:, b, a, d, c])
                    if b == d:
                        KK[:,a,b,c,d] += RR[:,a,c]
                    if a == c:
                        KK[:,a,b,c,d] += numpy.conj(RR[:,b,d])

                        
    return KK  


def _nsc_fintegral(tt, a, b, c, d, HH, gt):
    """Time dependent non-secular effective non-equilibrium Foerster integral
    
    
    Parameters
    ----------
    tt : numpy array
        Time 
        
    gtd : numpy array
        lineshape function of the donor transition

    gta : numpy array
        lineshape function of the acceptor transition 
        
    ed : float
        Energy of the donor transition
        
    ea : float
        Energy of the acceptor transition

    ld : float
        Reorganization energy of the donor             

    Returns
    -------
    ret : float
        The value of the Foerster integral            
    
    """

    Nt = tt.shape[0]
    hoft = numpy.zeros(Nt, dtype=COMPLEX)
    
    
    for ti in range(Nt):
        
        #
        # Here we calculate two-time integration kernel 
        #
        prod = _nsc_kernel_at_t(ti, tt, a, b, c, d, HH, gt)

        #
        # the kernel is integrated by splines
        #
        #inte = _integrate_kernel(tt, prod)
        inte = _integrate_kernel_to_t(ti, tt, prod)

        hoft[ti] = inte[ti]

    return hoft
    
 