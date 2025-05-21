# -*- coding: utf-8 -*-
import numpy
from ..core.managers import Manager
from ..qm.propagators.poppropagator import PopulationPropagator
from .response_implementations import get_implementation
from .. import REAL

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
        if self.diag in ["R1g","R4g", "R1f", "R1g_scM0g", "R1f_scM0g", "R1f_scM0e"]:
            self.rtype = "NR"
        elif self.diag in ["R2g", "R3g", "R2f", "R2g_scM0g", "R2f_scM0g", "R2f_scM0e"]:
            self.rtype = "R"
            
        self.func = get_implementation(self.diag)
        
        # what times to calculate for
        self.t1s = t1s
        self.t2s = t2s
        self.t3s = t3s

        # setting zero rates for the case they are not set from outside
        KK = numpy.zeros((system.Nb[1]+system.Nb[0],system.Nb[1]+system.Nb[0]),dtype=REAL)
        self.set_rate_matrix(KK)
    
    
    def calculate_matrix(self, t2): 
        """Calculates the matrix of response values in t1 and t3 times
        
        
        Parameters
        ----------
        
        t2 : float
            Waiting time for which the response is calculated
            
            
        
        """
        
        # identify the index of the present t2 time
        out = self.t2s.locate(t2)
        t2i = out[0]
        
        # population decay factors at t2
        U0t2 = self.U0_t2[:,t2i]
        
        return self.func(t2, self.t1s.data,
                             self.t3s.data, 
                             self.lab, self.sys, 
                             (self.U0_t1,self.U0_t3,U0t2, self.Uee[:,:,t2i]),
                             self.KK)


    def set_rwa(self, rwa):
        """Sets rotating wave approximation frequency
        
        """
        pass  # rwa is set through the system class, at least for now


    def set_rate_matrix(self, KK):
        """Sets the rate matrix and the corresponding evolution coefficients
        
        """
        if KK.shape[0] == KK.shape[1]:
            self.KK = KK
            
            self.U0_t2 = numpy.zeros((KK.shape[0], self.t2s.length),
                                  dtype=REAL)
            self.U0_t1 = numpy.zeros((KK.shape[0], self.t1s.length),
                                  dtype=REAL)
            self.U0_t3 = numpy.zeros((KK.shape[0], self.t3s.length),
                                  dtype=REAL)  
            
            if KK.shape[0] != self.sys.Ntot: 
                if self.sys.mult == 2:
                    self.U0fe_t3 = numpy.zeros((self.sys.Nb[2], 
                                                KK.shape[0], self.t3s.length),
                                                dtype=REAL)
                else:
                    raise Exception("Relaxation matrix has a wrong size: "
                                    +str(KK.shape[0]))                    
                    

        else:
            raise Exception("Square matrix must be submitted")
         
        # time independent rate matrix
        if len(KK.shape) == 2:
            
            #
            # Relaxation caused dephasing for single excitons
            #
            for aa in range(KK.shape[0]):
                if KK[aa,aa] <= 0.0:
                    self.U0_t2[aa,:] = numpy.exp(0.5*KK[aa,aa]*self.t2s.data)
                    self.U0_t1[aa,:] = 1.0 #numpy.exp(0.5*KK[aa,aa]*self.t1s.data)
                    self.U0_t3[aa,:] = 1.0 #numpy.exp(0.5*KK[aa,aa]*self.t3s.data)
                else:
                    raise Exception("Depopulation rate must be negative.")              
               
            #
            # Relaxation caused dephasing for double-excitons
            #
                    
               
            #
            # Finding population evolution matrix
            #
            
            # FIXME: Make sure it works with all t2s 
            prop = PopulationPropagator(self.t2s, self.KK)


            self.Uee, cor = prop.get_PropagationMatrix(self.t2s,
                                                  corrections=0)
              
            self.U1_t2 = cor[0]
              
            
        if len(KK.shape) == 3:
            # time dependent rate matrix
            pass


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
        
        