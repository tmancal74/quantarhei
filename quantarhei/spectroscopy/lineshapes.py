# -*- coding: utf-8 -*-
#
# This module defines several normalized lineshapes for 1D and 2D spectroscopy
#
import numpy
from scipy import special

from .. import REAL
from .. import COMPLEX

###############################################################################
#
#    1D absorptive lineshapes
#
###############################################################################

def gaussian(omega, cent, delta):
    """Normalized Gaussian line shape
    
    """
    return numpy.sqrt(numpy.log(2.0)/numpy.pi)\
                     *numpy.exp(-numpy.log(2.0)*((omega-cent)/delta)**2) \
                     /delta


def lorentzian(omega, cent, gamma):
    """Normalized Lorenzian line shape
    
    """
    return (gamma/numpy.pi)/((omega-cent)**2 + gamma**2)


def lorentzian_im(omega, cent, gamma):
    """Imaginary part of a normalized Lorenzian line shape
    
    """
    return 1j*((omega-cent)/numpy.pi)/((omega-cent)**2 + gamma**2)


def voigt(omega, cent, delta, gamma=0.0):
    """Normalized Voigt line shape for absorption
    
    """    
    z = (omega - cent + 1j*gamma)*numpy.sqrt(numpy.log(2.0))/delta
    
    return numpy.sqrt(numpy.log(2.0))*\
                      numpy.real(special.wofz(z)) \
                      /(numpy.sqrt(numpy.pi)*delta)


def cvoigt(omega, cent, delta, gamma=0.0):
    """Complex normalized Voigt line shape
    
    """
    a = (delta**2)/(4.0*numpy.log(2))
    z = (gamma - 1j*(omega - cent))/(2.0*numpy.sqrt(a))
    
    
    return numpy.real(special.erfcx(z))*numpy.sqrt(numpy.pi/a)/2.0


###############################################################################
#
#    2D lineshapes
#
###############################################################################


def gaussian2D(omega1, cent1, delta1, omega2, cent2, delta2, corr=0.0):
    """Two-dimensional complex Gaussian lineshape   
    
    """
    gamma1 = 0.0
    gamma2 = 0.0
    return voigt2D(omega1, cent1, delta1, gamma1,
                   omega2, cent2, delta2, gamma2, corr=corr)


def voigt2D(omega1, cent1, delta1, gamma1, 
            omega2, cent2, delta2, gamma2, corr=0.0):
    """Two-dimensional complex Voigt lineshape
    
    """
    if corr == 0.0:
        
        N1 = omega1.shape[0]
        N2 = omega2.shape[0]
        
        dat1 = cvoigt(omega1, cent1, delta1, gamma1)
        dat2 = cvoigt(omega2, cent2, delta2, gamma2)
        
        data = numpy.zeros((N1, N2), dtype=COMPLEX)  

        for k in range(N1):
            data[:, k] = dat1[k]*dat2[:]
        
    else:
        
        raise Exception("Not implemented yet")
        
    return data        


def lorentzian2D(omega1, cent1, gamma1, omega2, cent2, gamma2, corr=0.0):
    """Two-dimensional complex Lorentzian lineshape
    
    """
    
    if corr == 0.0:
        
        N1 = omega1.shape[0]
        N2 = omega2.shape[0]
        
        dat1 = lorentzian(omega1, cent1, gamma1) + \
               lorentzian_im(omega1, cent1, gamma1)
        dat2 = lorentzian(omega2, cent2, gamma2) + \
               lorentzian_im(omega2, cent2, gamma2)
        
        data = numpy.zeros((N1, N2), dtype=COMPLEX)  

        for k in range(N1):
            data[k, :] = dat1[k]*dat2[:]
        
    else:
        
        raise Exception("Not implemented yet")
        
    return data
