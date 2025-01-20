# -*- coding: utf-8 -*-
import numpy


"""






"""

class LiouvillePathway:
    """
    
    
    
    Parameters
    ----------
    
    
    ptype : str
        Type of the Liouville pathway
    
    
    """
    
    
    def __init__(self, ptype):
        
        self.ptype = ptype
        
        #
        # Determine the response type 
        # Types are:
        #   R ... rephasing
        #   NR .. non-rephasing
        #
        if ptype == "R2g":
            self.rtype = "R"
        else:
            # unknown response type
            self.rtpye = "U"
            
            
    def set_dipoles(self, d1, d2=None, d3=None, d4=None):
        """Sets the transition dipole moments of the response
        
        Parameters:
        -----------
        
        d1 : vector
            First transition dipole moment of the response. 
            
        d2 : vector or None
            The second transition dipole moment of the response. If None,
            it is set to the value of the first transition dipole moment.
        
        d3 : vector or None
            The third transition dipole moment of the response. If None,
            it is set to the value of the first transition dipole moment.

        d4 : vector or None
            The fourth transition dipole moment of the response. If None,
            it is set to the value of the first transition dipole moment.  
            
        """
        pass
        
    
    def set_frequencies(self, omega1, omega2):
        """Sets the frequencies of the response 
        
        
        """
        pass
    
    

            
            
            
            
class ResponseFunction(LiouvillePathway):
    """Non-linear response function
    
    
    
    
    """
    
    def calculate_matrix(self, lab, sys, it2, t1s, t3s, rwa):
        """Calculates the matrix of response values in t1 and t3 times
        
        
        Parameters
        ----------
        
        lab :
            
        sys : ... None
            System for which the response is calculated. If None, we assume
            that this response is defined as a stand alone. Its dipole factors
            and frequencies are defined by its constructor.
            
        it2 : integer
            Integer representing the t2 time
            
        t1s : 
            
        t2s :
            
        rwa: real
            Rotating wave frequency
            
        rmin : 
            
            
        
        """
        
        # create the dipole factor
        dip = 1.0
        
        # construct the phase factors in t1 and t3
        et13 = numpy.outer(numpy.exp(1j*rwa*t1s),numpy.exp(1j*rwa*t3s))
        
        t2 = t1s[it2]   # t2 time is taken from the t1 list
        
        return dip*self.func(t2,t1s,t3s,et13)
        
    
    def set_evaluation_function (self, func):
        """Sets the function to be evaluated to get the response matrix
        
        """
        self.func = func
        
  
    def set_auxliary_arguments(self, args):
        """Sets additional arguments and their values for evaluation calls
        
        """
        
        pass    
        
        