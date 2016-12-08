# -*- coding: utf-8 -*-

import numpy
import scipy.constants as const
from ..core.units import eps0_int

def dipole_dipole_interaction(r1, r2, d1, d2, epsr):
    """ Calculates interaction between two dipoles
    
    
    Parameters
    ----------
    
    r1 : array like
        position of molecule 1
        
    r2 : array like 
        position of molecule 2
        
    d1 : array like
        dipole moment of the first molecule
        
    d2 : array like
        dipole moment of the second molecule
        
    epsr : float
        relative permitivity of the environment
        
    
    """
    
    R = r1 - r2
    RR = numpy.sqrt(numpy.dot(R,R))
    
    prf = 1.0/(4.0*const.pi*eps0_int)
    
    cc = (numpy.dot(d1,d2)/(RR**3)
        - 3.0*numpy.dot(d1,R)*numpy.dot(d2,R)/(RR**5))
    
    return prf*cc/epsr    
    
    

