# -*- coding: utf-8 -*-
from .. import COMPLEX


def R1g(t2, t1, t3, lab, system):
    """ Returns a matrix of the respose function values for given t1 and t3 
    
    Parameters:
    -----------
    
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)
        
    t2 : float
        Value of the t2 (waiting) time of the response
        
    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)
        
    system : aggregate or molecule class
        An object storing all information about the system including 
        the values of the line shape functions.
    
    
    """
    import numpy as np

    gg = system.get_lineshape_functions()
    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()
    g = 1  # ground state index
    gg.create_data(reset={'t2':t2})

    Ne = En.shape[0]

    # dipole arrangemenent type: abba
    F4 = system.get_F4d('abba')
    dfac = np.einsum('i,abi->ab',lab.F4eM4,F4)

    ret = np.zeros((len(t1),len(t3)), dtype=COMPLEX)
    for aa in range(1,Ne):
        a = aa - 1
        for bb in range(1,Ne):
            b = bb -1

            ret += dfac[a,b]* \
            np.exp(
                -(np.einsum("i,ij", MM[b,a,:], gg[:,"t1"]))[:,None]
                + (np.einsum("i,ij", MM[b,a,:], gg[:,"t1+t2"]))[:,None]
                - (np.einsum("i,ijk", MM[a,a,:], gg[:,"t1+t2+t3"]))[:,:]
                - np.conj((np.einsum("i,i", MM[b,b,:], gg[:,"t2"]))[None,None])
                + np.conj((np.einsum("i,ij", MM[a,b,:], gg[:,"t2+t3"]))[None,:])
                - np.conj((np.einsum("i,ij", MM[a,b,:], gg[:,"t3"]))[None,:])
                -1j*(En[a]-En[g])*t1[:,None]
                -1j*(En[a]-En[b])*t2
                -1j*(En[a]-En[g])*t3[None,:]
            )

    return ret
           
                

def R2g(t2, t1, t3, system):
    pass


dc = dict()

dc["R1g"] = R1g
dc["R2g"] = R2g


def get_implementation(name):
    """Returns a dictionary of functions 
    
    
    """
    return dc[name]
