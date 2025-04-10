# -*- coding: utf-8 -*-
from .. import COMPLEX


def R1g(t2, t1, t3, lab, system, evol, KK):
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
    #g = 0  # ground state index
    gg.create_data(reset={'t2':t2}) 

    band1 = system.get_band(1)
    band0 = system.get_band(0)
    
    Ut2 = evol[2]
    
    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: abba
    F4 = system.get_F4d('abba')
    dfac = np.einsum('i,abi->ab',lab.F4eM4,F4)

    ret = np.zeros((len(t1),len(t3)), dtype=COMPLEX)
    for g in band0:
        for aa in band1:
            a = aa - 1
            for bb in band1:
                b = bb - 1
                    
                ret += dfac[a,b]*Ut2[a]*Ut2[b]* \
                    Ut1[a,:][:,None]*Ut3[a,:][None,:]*\
                np.exp(
                    -(np.einsum("i,ij", MM[b,a,:], gg[:,"t1"]))[:,None]
                    + (np.einsum("i,ij", MM[b,a,:], gg[:,"t1+t2"]))[:,None]
                    - (np.einsum("i,ijk", MM[a,a,:], gg[:,"t1+t2+t3"]))[:,:]
                    - np.conj((np.einsum("i,i", MM[b,b,:], gg[:,"t2"]))[None,None])
                    + np.conj((np.einsum("i,ij", MM[a,b,:], gg[:,"t2+t3"]))[None,:])
                    - np.conj((np.einsum("i,ij", MM[a,b,:], gg[:,"t3"]))[None,:])
                    -1j*(En[aa]-En[g]-rwa)*t1[:,None]
                    -1j*(En[aa]-En[bb])*t2
                    -1j*(En[aa]-En[g]-rwa)*t3[None,:]
                )

    return np.transpose(ret)
           
                

def R2g(t2, t1, t3, lab, system, evol, KK):
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
    #g = 0  # ground state index
    gg.create_data(reset={'t2':t2})

    band1 = system.get_band(1)
    band0 = system.get_band(0)
    
    Ut2 = evol[2]
    
    Ut1 = evol[0]
    Ut3 = evol[1]
    
    # dipole arrangemenent type: baba
    F4 = system.get_F4d('baba')
    dfac = np.einsum('i,abi->ab',lab.F4eM4,F4)

    ret = np.zeros((len(t1),len(t3)), dtype=COMPLEX)
    for g in band0:
        for aa in band1:
            a = aa - 1
            for bb in band1:  
                b = bb - 1
                
                ret += dfac[a,b]*Ut2[a]*Ut2[b]* \
                    Ut1[a,:][:,None]*Ut3[b,:][None,:]*\
                np.exp(
                    (np.einsum("i,i", MM[a,b,:], gg[:,"t2"]))[None,None]
                    - (np.einsum("i,ij", MM[b,b,:], gg[:,"t2+t3"]))[None,:]
                    - np.conj((np.einsum("i,ij", MM[a,a,:], gg[:,"t1+t2"]))[:,None])
                    - np.conj((np.einsum("i,ij", MM[a,b,:], gg[:,"t1"]))[:,None])
                    - np.conj((np.einsum("i,ij", MM[b,a,:], gg[:,"t3"]))[None,:])
                    + np.conj((np.einsum("i,ijk", MM[a,b,:], gg[:,"t1+t2+t3"]))[:,:])
                    -1j*(En[g]-En[aa]+rwa)*t1[:,None]
                    -1j*(En[bb]-En[aa])*t2
                    -1j*(En[bb]-En[g]-rwa)*t3[None,:]
                )
    
    return np.transpose(ret)                


def R3g(t2, t1, t3, lab, system, evol, KK):
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
    #g = 0  # ground state index
    gg.create_data(reset={'t2':t2})

    band1 = system.get_band(1)
    band0 = system.get_band(0)
    
    #Ut2 = evol[1]
    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: bbaa
    F4 = system.get_F4d('bbaa')
    dfac = np.einsum('i,abi->ab',lab.F4eM4,F4)

    ret = np.zeros((len(t1),len(t3)), dtype=COMPLEX)
    for g in band0:
        for aa in band1:
            a = aa - 1
            for bb in band1:
                b = bb - 1
  
                ret += dfac[a,b]* \
                    Ut1[a,:][:,None]*Ut3[b,:][None,:]*\
                np.exp(
                    -(np.einsum("i,ij", MM[b,b,:], gg[:,"t3"]))[None,:]
                    + np.conj((np.einsum("i,i", MM[b,a,:], gg[:,"t2"]))[None,None])
                    - np.conj((np.einsum("i,ij", MM[a,a,:], gg[:,"t1"]))[:,None])
                    - np.conj((np.einsum("i,ij", MM[a,b,:], gg[:,"t1+t2"]))[:,None])
                    - np.conj((np.einsum("i,ij", MM[b,a,:], gg[:,"t2+t3"]))[None,:])
                    + np.conj((np.einsum("i,ijk", MM[a,b,:], gg[:,"t1+t2+t3"]))[:,:])
                    -1j*(En[g]-En[aa]+rwa)*t1[:,None]
                    -1j*(En[g]-En[g])*t2
                    -1j*(En[bb]-En[g]-rwa)*t3[None,:]
                )

    return np.transpose(ret) 
                

def R4g(t2, t1, t3, lab, system, evol, KK):
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
    #g = 0  # ground state index
    gg.create_data(reset={'t2':t2})

    band1 = system.get_band(1)
    band0 = system.get_band(0)
    
    #Ut2 = evol[2]
    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: bbaa
    F4 = system.get_F4d('bbaa')
    dfac = np.einsum('i,abi->ab',lab.F4eM4,F4)

    ret = np.zeros((len(t1),len(t3)), dtype=COMPLEX)
    for g in band0:
        for aa in band1:
            a = aa - 1
            for bb in band1:
                b = bb - 1

                ret += dfac[a,b]* \
                    Ut1[a,:][:,None]*Ut3[b,:][None,:]*\
                np.exp(
                    -(np.einsum("i,i", MM[b,a,:], gg[:,"t2"]))[None,None]
                    - (np.einsum("i,ij", MM[a,a,:], gg[:,"t1"]))[:,None]
                    + (np.einsum("i,ij", MM[b,a,:], gg[:,"t1+t2"]))[:,None]
                    + (np.einsum("i,ij", MM[b,a,:], gg[:,"t2+t3"]))[None,:]
                    - (np.einsum("i,ij", MM[b,b,:], gg[:,"t3"]))[None,:]
                    - (np.einsum("i,ijk", MM[b,a,:], gg[:,"t1+t2+t3"]))[:,:]
                    -1j*(En[aa]-En[g]-rwa)*t1[:,None]
                    -1j*(En[g]-En[g])*t2[None,None]
                    -1j*(En[bb]-En[g]-rwa)*t3[None,:]
                )

    return np.transpose(ret)             


def R1f(t2, t1, t3, lab, system, evol, KK):
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
    gg.create_data(reset={'t2':t2})
    
    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()    

    band0 = system.get_band(0)
    band1 = system.get_band(1)
    band2 = system.get_band(2)

    N0 = system.Nb[0]
    N1 = system.Nb[1]
    
    Ut2 = evol[2]
    
    Ut1 = evol[0]
    Ut3 = evol[1]
    
    # dipole arrangemenent type: fbfaba
    F4 = system.get_F4d('fbfaba')
    dfac = np.einsum('i,fabi->fab',lab.F4eM4,F4)
    
    ret = np.zeros((len(t1),len(t3)), dtype=COMPLEX)
    for g in band0:
        for ff in band2:
            f = ff - N1 - N0
            p = ff - N0
            for aa in band1:
                a = aa - N0
                for bb in band1:
                    b = bb - N0

                    ret += -1.0*dfac[f,a,b]*Ut2[a]*Ut2[b]* \
                                            Ut1[a,:][:,None]*Ut3[b,:][None,:]*\
                    np.exp(
                        -(np.einsum("i,ij", MM[a,a,:], gg[:,"t1+t2"]))[:,None]
                        - (np.einsum("i,ij", MM[b,a,:], gg[:,"t1"]))[:,None]
                        - (np.einsum("i,ij", MM[b,a,:], gg[:,"t3"]))[None,:]
                        + (np.einsum("i,ij", MM[b,p,:], gg[:,"t3"]))[None,:]
                        + (np.einsum("i,ij", MM[p,a,:], gg[:,"t1+t2"]))[:,None]
                        + (np.einsum("i,ij", MM[p,a,:], gg[:,"t3"]))[None,:]
                        - (np.einsum("i,ij", MM[p,p,:], gg[:,"t3"]))[None,:]
                        + (np.einsum("i,ijk", MM[b,a,:], gg[:,"t1+t2+t3"]))[:,:]
                        - (np.einsum("i,ijk", MM[p,a,:], gg[:,"t1+t2+t3"]))[:,:]
                        + np.conj((np.einsum("i,i", MM[a,b,:], gg[:,"t2"])))
                        - np.conj((np.einsum("i,i", MM[p,b,:], gg[:,"t2"])))
                        - np.conj((np.einsum("i,ij", MM[b,b,:], gg[:,"t2+t3"]))[None,:])
                        + np.conj((np.einsum("i,ij", MM[p,b,:], gg[:,"t2+t3"]))[None,:])
                        -1j*(En[aa]-En[g]-rwa)*t1[:,None]
                        -1j*(En[aa]-En[bb])*t2[None,None]
                        -1j*(En[ff]-En[bb]-rwa)*t3[None,:]
                    )
    
    return np.transpose(ret)


def R2f(t2, t1, t3, lab, system, evol, KK):
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

    gg.create_data(reset={'t2':t2})

    band0 = system.get_band(0)
    band1 = system.get_band(1)
    band2 = system.get_band(2)

    N0 = system.Nb[0]
    N1 = system.Nb[1]
    
    Ut2 = evol[2]
    
    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: fafbba
    F4 = system.get_F4d('fafbba')
    dfac = np.einsum('i,fabi->fab',lab.F4eM4,F4)

    ret = np.zeros((len(t1),len(t3)), dtype=COMPLEX)
    
    for g in band0:
        for ff in band2:
            f = ff - N1 - N0
            p = ff - N0
            for aa in band1:
                a = aa - N0
                for bb in band1:
                    b = bb - N0
    
                    ret += -1.0*dfac[f,a,b]*Ut2[a]*Ut2[b]* \
                                            Ut1[a,:][:,None]*Ut3[a,:][None,:]*\
                        np.exp(
                        -(np.einsum("i,i", MM[b,b,:], gg[:,"t2"]))[None,None]
                        + (np.einsum("i,i", MM[p,b,:], gg[:,"t2"]))[None,None]
                        + (np.einsum("i,ij", MM[a,b,:], gg[:,"t2+t3"]))[None,:]
                        - (np.einsum("i,ij", MM[a,b,:], gg[:,"t3"]))[None,:]
                        + (np.einsum("i,ij", MM[a,p,:], gg[:,"t3"]))[None,:]
                        - (np.einsum("i,ij", MM[p,b,:], gg[:,"t2+t3"]))[None,:]
                        + (np.einsum("i,ij", MM[p,b,:], gg[:,"t3"]))[None,:]
                        - (np.einsum("i,ij", MM[p,p,:], gg[:,"t3"]))[None,:]
                        - np.conj((np.einsum("i,ij", MM[a,b,:], gg[:,"t1"]))[:,None])
                        + np.conj((np.einsum("i,ij", MM[a,b,:], gg[:,"t1+t2"]))[:,None])
                        - np.conj((np.einsum("i,ij", MM[a,p,:], gg[:,"t1+t2"]))[:,None])
                        - np.conj((np.einsum("i,ijk", MM[a,a,:], gg[:,"t1+t2+t3"]))[:,:])
                        + np.conj((np.einsum("i,ijk", MM[a,p,:], gg[:,"t1+t2+t3"]))[:,:])
                        -1j*(En[g]-En[aa]+rwa)*t1[:,None]
                        -1j*(En[bb]-En[aa])*t2[None,None]
                        -1j*(En[ff]-En[aa]-rwa)*t3[None,:]
                )

    return np.transpose(ret)


def R1g_RL0(t2, t1, t3, lab, system, evol, KK):
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
    #g = 0  # ground state index
    gg.create_data(reset={'t2':t2}) 

    band1 = system.get_band(1)
    band0 = system.get_band(0)
    
    Ut2 = evol[3]
    
    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: abba
    F4 = system.get_F4d('bbaa')
    dfac = np.einsum('i,abi->ab',lab.F4eM4,F4)

    ret = np.zeros((len(t1),len(t3)), dtype=COMPLEX)
    for g in band0:
        for aa in band1:
            a = aa - 1
            for bb in band1:
                b = bb - 1
                    
                if a != b:
                    ret += dfac[a,b]*Ut2[b,a]* \
                        Ut1[a,:][:,None]*Ut3[b,:][None,:]*\
                    np.exp(
                        -(np.einsum("i,ij", MM[a,a,:], gg[:,"t1"]))[:,None]
                        - np.conj((np.einsum("i,ij", MM[b,b,:], gg[:,"t3"]))[None,:])
                        -1j*(En[aa]-En[g]-rwa)*t1[:,None]
                        -1j*(En[bb]-En[g]-rwa)*t3[None,:]
                    )

    return np.transpose(ret)
           
                

def R2g_RL0(t2, t1, t3, lab, system, evol, KK):
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
    #g = 0  # ground state index
    gg.create_data(reset={'t2':t2})

    band1 = system.get_band(1)
    band0 = system.get_band(0)
    
    Ut2 = evol[3]
    
    Ut1 = evol[0]
    Ut3 = evol[1]
    
    # dipole arrangemenent type: baba
    F4 = system.get_F4d('bbaa')
    dfac = np.einsum('i,abi->ab',lab.F4eM4,F4)

    ret = np.zeros((len(t1),len(t3)), dtype=COMPLEX)
    for g in band0:
        for aa in band1:
            a = aa - 1
            for bb in band1:  
                b = bb - 1
                
                if a != b:
                    ret += dfac[a,b]*Ut2[b,a]* \
                        Ut1[a,:][:,None]*Ut3[b,:][None,:]*\
                    np.exp(
                        - np.conj((np.einsum("i,ij", MM[a,a,:], gg[:,"t1"]))[:,None])
                        - np.conj((np.einsum("i,ij", MM[b,b,:], gg[:,"t3"]))[None,:])
                        -1j*(En[g]-En[aa]+rwa)*t1[:,None]
                        -1j*(En[bb]-En[g]-rwa)*t3[None,:]
                    )
    
    return np.transpose(ret)                



def R1f_RL0(t2, t1, t3, lab, system, evol, KK):
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
    gg.create_data(reset={'t2':t2})
    
    # Mx = system.get_participation()
    MM = system.get_weighted_participation()
    En = system.get_eigenstate_energies()
    rwa = system.get_RWA_suggestion()    

    band0 = system.get_band(0)
    band1 = system.get_band(1)
    band2 = system.get_band(2)

    N0 = system.Nb[0]
    N1 = system.Nb[1]
    
    Ut2 = evol[3]
    
    Ut1 = evol[0]
    Ut3 = evol[1]
    
    # dipole arrangemenent type: fbfaba
    F4 = system.get_F4d('fbfbaa')
    dfac = np.einsum('i,fabi->fab',lab.F4eM4,F4)
    
    ret = np.zeros((len(t1),len(t3)), dtype=COMPLEX)
    for g in band0:
        for ff in band2:
            f = ff - N1 - N0
            p = ff - N0
            for aa in band1:
                a = aa - N0
                for bb in band1:
                    b = bb - N0
                    
                    if a != b:

                        ret += -1.0*dfac[f,a,b]*Ut2[a,b]* \
                                                Ut1[a,:][:,None]*Ut3[b,:][None,:]*\
                        np.exp(
                            - (np.einsum("i,ij", MM[a,a,:], gg[:,"t1"]))[:,None]
                            - (np.einsum("i,ij", MM[p,p,:], gg[:,"t3"]))[None,:]
                            -1j*(En[aa]-En[g]-rwa)*t1[:,None]
                            -1j*(En[ff]-En[bb]-rwa)*t3[None,:]
                        )
    
    return np.transpose(ret)


def R2f_RL0(t2, t1, t3, lab, system, evol, KK):
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

    gg.create_data(reset={'t2':t2})

    band0 = system.get_band(0)
    band1 = system.get_band(1)
    band2 = system.get_band(2)

    N0 = system.Nb[0]
    N1 = system.Nb[1]
    
    Ut2 = evol[3]
    
    Ut1 = evol[0]
    Ut3 = evol[1]

    # dipole arrangemenent type: fafbba
    F4 = system.get_F4d('fbfbaa')
    dfac = np.einsum('i,fabi->fab',lab.F4eM4,F4)

    ret = np.zeros((len(t1),len(t3)), dtype=COMPLEX)
    
    for g in band0:
        for ff in band2:
            f = ff - N1 - N0
            p = ff - N0
            for aa in band1:
                a = aa - N0
                for bb in band1:
                    b = bb - N0
                    
                    if b != a:
    
                        ret += -1.0*dfac[f,a,b]*Ut2[b,a]* \
                                                Ut1[a,:][:,None]*Ut3[a,:][None,:]*\
                            np.exp(
                            - (np.einsum("i,ij", MM[p,p,:], gg[:,"t3"]))[None,:]
                            - np.conj((np.einsum("i,ij", MM[a,a,:], gg[:,"t1"]))[:,None])
                            -1j*(En[g]-En[aa]+rwa)*t1[:,None]
                            -1j*(En[ff]-En[bb]-rwa)*t3[None,:]
                        )

    return np.transpose(ret)




#
# REGISTRATION COE
#

dc = dict()

dc["R1g"] = R1g
dc["R2g"] = R2g
dc["R3g"] = R3g
dc["R4g"] = R4g
dc["R1f"] = R1f
dc["R2f"] = R2f
dc["R1g_RL0"] = R1g_RL0
dc["R2g_RL0"] = R2g_RL0
dc["R1f_RL0"] = R1f_RL0
dc["R2f_RL0"] = R2f_RL0


def get_implementation(name):
    """Returns a dictionary of functions 
    
    
    """
    return dc[name]
