# -*- coding: utf-8 -*-

import numpy


from ...core.parallel import block_distributed_range
from ...core.parallel import start_parallel_region
from ...core.parallel import close_parallel_region
from ...core.parallel import distributed_configuration



def ssRedfieldRateMatrix(Na, Nk, KI, cc, rtol, werror, RR):
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
    
    # loop over components

    dc = distributed_configuration() # Manager().get_DistributedConfiguration()
    
    #
    #  PARALLELIZED
    #
    start_parallel_region()
    for k in block_distributed_range(0, Nk):
    #for k in range(Nk):
        
        # interaction operator
        KK = KI[k,:,:]

        for i in range(Na):
            for j in range(Na):
                
                # calculate rates, i.e. off diagonal elements
                if i != j:                                
                            
                    RR[i,j] += (cc[k,i,j]*KK[i,j]*KK[j,i])
                    
    
    # FIXME: parallelization ignores werror
    dc.allreduce(RR, operation="sum")        
    close_parallel_region()
    
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
                    
