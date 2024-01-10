# -*- coding: utf-8 -*-
"""
    
    State vector of the open quantum system


"""
import numpy

from ... import REAL, COMPLEX

class OQSStateVector():
    """Represents a quantum mechanical state of an open system 
    
    In this representation we keep state vector coefficients separate
    from the state of the bath representation.


    Parameters
    ----------


    Examples
    --------

    >>> psi = OQSStateVector(2)
    >>> print(psi.dim)
    2


    >>> vec = numpy.zeros((1,3), dtype=REAL)
    >>> psi = OQSStateVector(data=vec)
    Traceback (most recent call last):
    ...
    Exception: Data has to be a vector


    """    
    
    
    def __init__(self, dim=None, data=None):
        
        self._initialited = False
        
        if data is not None:
            
            ddat = numpy.array(data)
            
            if len(ddat.shape) > 1:
                raise Exception("Data has to be a vector")
                
            if dim != ddat.shape[0]:
                print("Dimension specification differes from data: ignored.")
                
            self.data = ddat 
            self.dim = ddat.shape[0]
            self._initialited = True
            
        elif dim is not None:
            
            self.dim = dim
            self.data = numpy.zeros(self.dim, dtype=REAL)
            self._initialited = True
            
            
    def norm(self):
        """Norm of the state vector 
        
        """
        
        return numpy.dot(self.data,self.data)


    def puredot(self, psi):
        """Dot product concerning only the system part
        
        """
        
        return numpy.dot(self.data, psi.data)
        
        