# -*- coding: utf-8 -*-

"""
*******************************************************************************

      MODIFIED REDFIELD RATE MATRIX

*******************************************************************************
"""  

import numpy

from ....core.implementations import implementation
from ....core.units import cm2int
from ....core.units import kB_intK

from ...hilbertspace.hamiltonian import Hamiltonian
from ...liouvillespace.systembathinteraction import SystemBathInteraction

class ModifiedRedfieldRateMatrix:
    """Modifield Redfield relaxation rate matrix
    
    Modified Redfield population relaxation rate matrix is calculated from the 
    Hamiltonian and system-system bath interation. The bath
    correlation functions are Fourier transformed by FFT and the negative
    frequency part is calculated from the expected thermodynamic
    symmetry. 
    
    Parameters
    ----------
    
    ham : Hamiltonian
        Hamiltonian object 
        
    sbi : SystemBathInteraction
        SystemBathInteraction object
        
    initialize : bool (default True)
        If true, the rates will be calculated when the object is created
        
    cutoff_time : float
        If cutoff time is specified, the tensor is integrated only up to the
        cutoff time
    
    
    """
    
    def __init__(self, ham, sbi, initialize=True, cutoff_time=None):
        
        if not isinstance(ham,Hamiltonian):
            raise Exception("First argument must be a Hamiltonian")
            
        if not isinstance(sbi,SystemBathInteraction):
            raise Exception
            
        self._is_initialized = False            
        self._has_cutoff_time = False
        
        if cutoff_time is not None:
            self.cutoff_time = cutoff_time
            self._has_cutoff_time = True            
            
        self.ham = ham
        self.sbi = sbi
        
        if initialize: 
            self._set_rates()          
            self._is_initialized = True
            
    def _set_rates(self):
        pass


def ssModifiedRedfieldRateMatrix(Na, Nc, Nt, Ee, SS, prt, gg, hh, cc, tt, ls,
                                 rtol, werror, RR):
    """Standard redfield rates
    
    
    Parameters
    ----------
    
    Na : integer
        Rank of the rate matrix, number of excitons
        
    Nc : integer
        Number of components of the interaction Hamiltonian (number of sites
        or number of distinct correlation functions)
    
    Ee : float array
        Eigen energies of the Hamiltonian
        
    SS : float array
        Transformation matrix 
        (components of the interaction operator in the exciton basis)
        
    prt : integer array
        Pointer between site index and correlation function index
    
    gg : float array
        Line shape functions
        
    hh : float array
        derivative of the line shape function
    
    cc : float array
        second derivatives of the line shape functions (correlation functions)
        
    tt : float array
        values of time
        
    ls : float array
        reorganization energies of the sites
        
    RR : real array
        Relaxation rate matrix (to be calculated and returned)
    
    rtol : float array
        tolerances
        
    werror : integer array
        warnings and errors
        
    """
    

    pass
                    
                    
