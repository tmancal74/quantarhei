# -*- coding: utf-8 -*-

import numpy


from ...core.parallel import block_distributed_range
from ...core.parallel import start_parallel_region
from ...core.parallel import close_parallel_region
from ...core.parallel import distributed_configuration



def ssRedfieldRateMatrix(Na, Nk, KI, cc, rtol, werror, RR, corrM = None):
    """Standard redfield rates
    
    
    Parameters
    ----------
    
    Na : integer
        Rank of the rate matrix, number of excitons
        
    Nk : integer
        Number of components of the interaction Hamiltonian
    
    KI : float array
        System parts of the interaction Hamiltonian components
        
    cc : float array
        Half of the Fourier transform of the correlation functions 
        at all transition frequencies
        
    RR : real array
        Relaxation rate matrix (to be calculated and returned)
    
    """
    
    if corrM is None:
        corrM = numpy.identity(Nk)
    
    # loop over components

    #dc = distributed_configuration() # Manager().get_DistributedConfiguration()
    
    #
    #  SERIAL VERSION
    #
    for i in range(Na):
        for j in range(Na):
            
            # calculate rates, i.e. off diagonal elements
            if i != j:                                
                        
                RR[i,j] += numpy.dot(cc[:,i,j]*KI[:,i,j],numpy.dot(corrM,KI[:,j,i]))
                    
    
    # FIXME: parallelization ignores werror
    #dc.allreduce(RR, operation="sum")        
    #close_parallel_region()
    
    #
    #  END PARALLELIZED
    #                   
     
    # calculate the diagonal elements (the depopulation rates)            
    for i in range(Na):
        for j in range(Na):
            
            if i != j:
                if RR[i,j] < 0.0:
                    werror[0] = -1
                    if numpy.abs(RR[i,j]) < rtol:
                        RR[i,j] = 0.0
                    else:
                        werror[1] = -1
                    
            if i != j:
                RR[j,j] -= RR[i,j]
                    
