# -*- coding: utf-8 -*-
import numpy
from ..core.managers import Manager
from .response_implementations import get_implementation

"""
    This packege contains two methodologies for calculating non-linear
    response. The one based on the LiouvillPathway class is now deprecated.
    
    The new methodology is based on NonLinearResponse class





"""
class NonLinearResponse:
    """Non-linear response function 
    
    
    
    
    """
    
    def __init__(self, lab, system, diagram, t1s, t2s, t3s):
        
        # info about pulse polarizations
        self.lab = lab
        
        # info about energies, dipolemoments and rwa 
        self.sys = system
        
        # which response to calculate; the function to calculate the respose
        self.diag = diagram 
        self.func = get_implementation(self.diag)
        
        # what times to calculate for
        self.t1s = t1s
        self.t2s = t2s
        self.t3s = t3s
        
        # devise a way to pass it to the response calculation
    
    
    def calculate_matrix(self, t2): 
        """Calculates the matrix of response values in t1 and t3 times
        
        
        Parameters
        ----------
        
        t2 : float
            Waiting time for which the response is calculated
            
            
        
        """
        return self.func(t2, self.t1, self.t3, self.lab, self.system)





###############################################################################
#
#   DEPRECATED CODE BELOW
#
###############################################################################


class LiouvillePathway:
    """
    
    
    
    Parameters
    ----------
    
    
    ptype : str
        Type of the Liouville pathway
    
    
    """
    
    
    def __init__(self, ptype):
        
        self.ptype = ptype
        
        self.sign = 1.0
        #
        # Determine the response type 
        # Types are:
        #   R ... rephasing
        #   NR .. non-rephasing
        #
        if ptype == "R2g":
            self.rtype = "R"
        elif ptype == "R1g":
            self.rtype = "NR"
        elif ptype == "R3g":
            self.rtype = "R"
        elif ptype == "R4g":
            self.rtype = "NR"
        elif ptype == "R1f":
            self.rtype = "R"
            self.sign = -1.0
        elif ptype == "R2f":
            self.rtype = "NR"
            self.sign = -1.0
        else:
            # unknown response type
            self.rtype = "U"
            
        self._frequencies_set = False
        self._rwa_set = False
        
        self.F4n = numpy.zeros(3)
        
            
            
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
        d = numpy.zeros((4,3), dtype=float)
        
        d[0,:] = d1
        if d2 is None:
            d[1,:] = d1
        else:
            d[1,:] = d2
        if d3 is None:
            d[2,:] = d1
        else:
            d[2,:] = d3
        if d4 is None:
            d[3,:] = d1
        else:
            d[3,:] = d4
        
        self.F4n[0] = numpy.dot(d[3,:],d[2,:])*numpy.dot(d[1,:],d[0,:])
        self.F4n[1] = numpy.dot(d[3,:],d[1,:])*numpy.dot(d[2,:],d[0,:])
        self.F4n[2] = numpy.dot(d[3,:],d[0,:])*numpy.dot(d[2,:],d[1,:]) 
        self.dd = d
        
    
    def set_frequencies(self, omega1, omega3):
        """Sets the frequencies of the response 
        
        
        """
        
        if self._frequencies_set:
            raise Exception("Frequencies of are already set.")
            
        self._omega1 = Manager().convert_energy_2_internal_u(omega1) 
        self._omega3 = Manager().convert_energy_2_internal_u(omega3)
        
        self._frequencies_set = True
    
    
    def get_frequencies(self):
        """Returns the two main frequencies of the response
        
        """
        
        if self._frequencies_set:
            fr = (Manager().convert_energy_2_current_u(self._omega1),
                  Manager().convert_energy_2_current_u(self._omega3))
        
            return fr
        else:
            raise Exception("Frequencies not set.")
    
    
    def set_rwa(self, rwa):
        """Sets the RWA frequency
        
        """
        
        if not self._frequencies_set:
            raise Exception("Frequencies must be set before setting RWA.")
            
        if not self._rwa_set:
            self._rwa = Manager().convert_energy_2_internal_u(rwa)
            self._omega1 = self._omega1 - self._rwa
            self._omega3 = self._omega3 - self._rwa
            
        else:
            raise Exception("RWA cannot be set twice. Reset first.")
            
        self._rwa_set = True


    def reset_rwa(self):
        """Resets the RWA setting
        

        """
        if self._frequencies_set and self._rwa_set:
            self._omega1 = self._omega1 + self._rwa
            self._omega3 = self._omega3 + self._rwa   
            
            self._rwa = 0.0
            self._rwa_set = False
        
        
    def get_rwa(self):
        """Returns the RWA frequency
        
        
        """
        return Manager().convert_energy_2_internal_u(self._rwa)

            
            
            
            
class ResponseFunction(LiouvillePathway):
    """Non-linear response function
    
    
    
    
    """
    
    def calculate_matrix(self, lab, sys, t2, t1s, t3s, rwa):
        """Calculates the matrix of response values in t1 and t3 times
        
        
        Parameters
        ----------
        
        lab :
            
        sys : ... None
            System for which the response is calculated. If None, we assume
            that this response is defined as a stand alone. Its dipole factors
            and frequencies are defined by its constructor.
            
        t2 : real
            The t2 time
            
        t1s : 
            
        t2s :
            
        rwa: real
            Rotating wave frequency
            
        rmin : 
            
            
        
        """
        
        om3 = self._omega3
        om1 = self._omega1
        
        
        # create the dipole factor
        dip = self.sign*numpy.dot(lab.F4eM4,self.F4n)
        
        # construct the phase factors in t1 and t3
        if self.rtype == "R":
            et13 = numpy.outer(numpy.exp(-1j*om3*t3s),
                               numpy.exp(1j*om1*t1s))
        elif self.rtype == "NR":
            et13 = numpy.outer(numpy.exp(-1j*om3*t3s),
                               numpy.exp(-1j*om1*t1s))            
        
        
        val = self.func(t2, t3s, t1s, *self.args)
        
        return dip*val*et13
        
    
    def set_evaluation_function (self, func):
        """Sets the function to be evaluated to get the response matrix
        
        """
        self.func = func
        
  
    def set_auxliary_arguments(self, args):
        """Sets additional arguments and their values for evaluation calls
        
        
        Parameters:
        -----------
        
        args : tuple
            A tuple of arguments that will be passed to the function after
            the expected standard arguments
            
        """
        
        self.args = args    
        
        