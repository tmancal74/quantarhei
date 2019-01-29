# -*- coding: utf-8 -*-
import numpy

from ...core.saveable import Saveable
from ...core.dfunction import DFunction


class fcstorage(Saveable):
    """FC factor look-up class
    
    Once Frank-Condon factors for some value of the shift are calculated
    they can be stored here, and retrieved when needed again
    """
    
    def __init__(self):
        """ Constructor """
        self._shifts = []
        self._fcs = []
        
    
    def lookup(self,shift):
        """ Returns true if the FC factors for a given shift are available """
        if self._shifts.count(shift) > 0:
            return True
        return False
    
    def index(self,shift):
        """ Returns an index of the FC factors with a given shift """
        return self._shifts.index(shift)
    
    def add(self,shift,fcmatrix):
        """ Adds the matrix of FC factors to the storage """
        self._shifts.append(shift)
        self._fcs.append(fcmatrix)
        
    def get(self,ii):
        """ Returns a stored FC matrix """
        return self._fcs[ii]
        
    


class operator_factory(Saveable):
    """Class providing useful operators
    
    
    
    Creation and anihilation operators    
    """    
    def __init__(self, N=100):
        # we choose a number of state 
        # to represent all operators
        self.N = N  
        
        
    def anihilation_operator(self):
        N = self.N
        aa = numpy.zeros((N,N),dtype=numpy.float) # matrix N x N full of zeros

        for ng in range(N):
            for mg in range(N):
                if ng == mg - 1:
                    aa[ng,mg] = numpy.sqrt(numpy.real(mg))
            
        return aa
    
    def creation_operator(self):
        N = self.N
        ad = numpy.zeros((N,N),dtype=numpy.float)
            
        for ng in range(N):
            for mg in range(N):
                if ng == mg + 1:
                    ad[ng,mg] = numpy.sqrt(numpy.real(mg+1))

        return ad
        
    def shift_operator(self,dd_):
        """Calculates the Shift Operator based on the size N_ of the basis
        of states and the shift dd_."""
        
        N_ = self.N
        aa = self.anihilation_operator()
        ad = self.creation_operator()
        
        # construct the Shift Operator
        Dd_large = numpy.zeros((N_,N_),dtype=numpy.float)
        Dd_large = dd_*(ad-aa)/numpy.sqrt(2.0)

        # Diagonalize and obtain transformation matrix
        A,S = numpy.linalg.eig(Dd_large)
        S1 = numpy.linalg.inv(S)
    
        # Exponentiate
        Dd_large = numpy.diag(numpy.exp(A))
    
        # Transform back and reduce to the lower number of states
        return numpy.real(numpy.dot(S,numpy.dot(Dd_large,S1)))
        
        
    def unity_operator(self):
        
        ones = numpy.ones(self.N, dtype=numpy.float)
        ret = numpy.diag(ones)
        return ret
        

    
class qrepresentation:
    """Coordinate representation of the HO wavefunctions
    
    
    """
    
    def __init__(self, qaxis):
        
        self.qaxis = qaxis
        self.ho_eigenfce_generated = False
        
    
    def generate_ho_eigenfunctions(self):
        """Generated q-representation of HO eigenfunctions
        
        """
        self.ho_eigenfce_generated = True


    def get_ho_eigenfunction(self, N):
        """Returns q-representation of HO eigenfunction
        
        """
        if self.ho_eigenfce_generated:
            pass
        else:
            raise Exception("HO Eigenfunctions must be generated first")


    def get_ho_ground_state(self):
        """Returns the ground state wavefunction of the Harmonic oscillator
        
        
        """
        
        data = numpy.exp(-(self.qaxis.data**2)/2)\
            /numpy.sqrt(numpy.sqrt(numpy.pi))
        psi0 = DFunction(x=self.qaxis, y=data)
        
        return psi0
    
    
    def get_coherent_state(self, alpha):
        """Returns q-representation of a coherent state with a given alpha
        
        """
        
        ar = numpy.real(alpha)
        ai = numpy.imag(alpha)
        sq2 = numpy.sqrt(2.0)
        q = self.qaxis.data
        data = numpy.exp(-((q+sq2*ar)**2)/2)*\
                            numpy.exp(-1j*ar*ai + 1j*sq2*ai*q)\
                            /numpy.sqrt(numpy.sqrt(numpy.pi))
                            
        psi_alpha = DFunction(x=self.qaxis, y=data)        
        
        return psi_alpha


    def get_probability_distribution(self, wfce):
        """Returns probability distribution for a given wavefunction
        
        """
        
        return DFunction(x=wfce.axis, y=numpy.abs(wfce.data)**2)
        
        
