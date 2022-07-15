# -*- coding: utf-8 -*-
import numpy
#import matplotlib.pyplot as plt
import scipy.interpolate as interp

from .tdfoerstertensor import TDFoersterRelaxationTensor
from .tdfoerstertensor import _td_reference_implementation
from ... import COMPLEX

class NEFoersterRelaxationTensor(TDFoersterRelaxationTensor):
    """Weak resonance coupling relaxation tensor by Foerster theory
    
    Non-equilibrium version according to 
    
    J. Seibt and T. Mančal, J. Chem. Phys. 146 (2017) 174109
    
    
    """
    def __init__(self, ham, sbi, initialize=True, cutoff_time=None):
        """Initiation is the same as for the TDFoerster tensor
        
        """
        super().__init__(ham, sbi, initialize, cutoff_time)
        


    def td_reference_implementation(self, Na, Nt, HH, tt, gt, ll):
        """ Overloaded implementation method replacing integration kernel
        
        """
        
        self.Iterm = _initial_term(Na, Nt, HH, tt, gt, ll)
        self.has_Iterm = True
        
        return _td_reference_implementation(Na, Nt, HH, tt,
                                            gt, ll, _ne_fintegral)
    


def _initial_term(Na, Nt, HH, tt, gt, ll):
    """ Inhomogeneous (initial condition) term 
        of the non-equilibrium Foerster theory
        
    """
    return numpy.zeros((Nt, Na, Na), dtype=COMPLEX)


def _kernel_at_t(ti, tt, gtd, gta, ed, ea, ld):
    """ Two-time kernel to be integrated 
    
    
    """ 
    Nt = tt.shape[0]
    gtd_i = gtd[0:ti]
    gtd_m = numpy.zeros(Nt, dtype=COMPLEX)
    gtd_m[0:ti] = numpy.flip(gtd_i) 
    
    prod = numpy.exp(-gtd-gta +1j*(ed-ea)*tt) \
        *numpy.exp(-2.0*1j*numpy.imag(gtd_m)) \
        *numpy.exp(2.0*1j*numpy.imag(gtd[ti]))
        
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
    

def _integrate_kernel_to_t(ti, tt, fce):
    """ Spline partial integration of a complex function
    
    """
    ti_min = 10
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



 