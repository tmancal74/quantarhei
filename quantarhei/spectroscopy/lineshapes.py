# -*- coding: utf-8 -*-
#
# This module defines several normalized lineshape functions
#
import numpy
from scipy import special

from .. import REAL

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
    

def voigt(omega, cent, delta, gamma=0.0):
    """Normalized Voigt line shape
    
    """    
    z = (omega - cent + 1j*gamma)*numpy.sqrt(numpy.log(2.0))/delta
    
    return numpy.sqrt(numpy.log(2.0))*\
                      numpy.real(special.wofz(z)) \
                      /(numpy.sqrt(numpy.pi)*delta)


###############################################################################
#
#    2D lineshapes
#
###############################################################################


def gaussian2D(omega1, cent1, delta1, omega2, cent2, delta2, corr=0.0):
    """Two-dimensional Gaussian lineshape
    
    
    
    """
    
    if corr == 0.0:
        
        N1 = omega1.shape[0]
        N2 = omega2.shape[0]
        
        dat1 = gaussian(omega1, cent1, delta1)
        dat2 = gaussian(omega2, cent2, delta2)
        
        data = numpy.zeros((N1, N2), dtype=REAL)
        
        for k in range(N1):
            data[k, :] = dat1[k]*dat2[:]
        
    else:
        
        raise Exception("Not implemented yet")
        
    return data