# -*- coding: utf-8 -*-
"""
    Time-dependent third order non-linear response of a three band
    multi-level system


"""
#import numpy

#import aceto.nr3td_fic as nr3td_fic
    

#
#
#  Non-transfering pathways
#
#
#

def nr3_r1g(lab, sys, it2, t1s, t3s, rwa, rmin, resp):
    """ Calculates R2g response function
    
    
    """
   
    
    nr3td_fic.nr3_r1g_fic(lab.orient_aver, sys.Ns, sys.om01, sys.nn01,
                        sys.dd01, sys.Kd01, sys.Kd11, sys.gofts, sys.fptn, 
                        sys.SS1, it2+1, t1s, t3s, rwa, rmin, resp)   
  

    
def nr3_r2g(lab, sys, it2, t1s, t3s, rwa, rmin, resp):
    """ Calculates R2g response function
    
    

    Paramaters
    ----------

    lab : lab_settings
        Laboratory settings (polarizations of laser beams etc.) expressed 
        through the lab_setting class
        
    sys : band_system
        System to be calculated on expressed through the band_system class
        
    it2 : integer
        Index of the so-called waiting or population time of non-linear 
        spectroscopic techniques
        
    t1s : float array
        Values of t1 time for which response should be calculated
        
    t3s : float array
        Values of t3 time for which response shuld be calculated
        
    rwa: float
        Rotating wave frequency
        
    rmin: float
        Minimal value of the dipole prefactor, relative to its max value,
        which is taken into account
        
    resp : complex 2d array
        Non-linear response 
        
    """
#
#    For debugging, check if all arrays are fortran continuous
#
#    print(lab.orient_aver.flags['F_CONTIGUOUS'])
#    print(sys.Ns.flags['F_CONTIGUOUS'])
#    print(sys.om01.flags['F_CONTIGUOUS'])
#    print(sys.nn01.flags['F_CONTIGUOUS'])
#    print(sys.dd01.flags['F_CONTIGUOUS'])
#    print(sys.Kd01.flags['F_CONTIGUOUS'])
#    print(sys.Kd11.flags['F_CONTIGUOUS'])
#    print(t1s.flags['F_CONTIGUOUS'])
#    print(t3s.flags['F_CONTIGUOUS'])
#    print(resp.flags['F_CONTIGUOUS'])
    
    
    nr3td_fic.nr3_r2g_fic(lab.orient_aver, sys.Ns, sys.om01, sys.nn01,
                        sys.dd01, sys.Kd01, sys.Kd11, sys.gofts, sys.fptn, 
                        sys.SS1, it2+1, t1s, t3s, rwa, rmin, resp)


def nr3_r2gt10(lab, sys, it2, t3s, rwa, rmin, resp):
    """ Calculates R2g response function
    
    

    Paramaters
    ----------

    lab : lab_settings
        Laboratory settings (polarizations of laser beams etc.) expressed 
        through the lab_setting class
        
    sys : band_system
        System to be calculated on expressed through the band_system class
        
    it2 : integer
        Index of the so-called waiting or population time of non-linear 
        spectroscopic techniques
        
    t3s : float array
        Values of t3 time for which response shuld be calculated
        
    rwa: float
        Rotating wave frequency
        
    rmin: float
        Minimal value of the dipole prefactor, relative to its max value,
        which is taken into account
        
    resp : complex 2d array
        Non-linear response 
        
    """
 
    
    nr3td_fic.nr3_r2gt10_fic(lab.orient_aver, sys.Ns, sys.om01, sys.nn01,
                        sys.dd01, sys.Kd01, sys.Kd11, sys.gofts, sys.fptn, 
                        sys.SS1, it2+1, t3s, rwa, rmin, resp)



def nr3_r3g(lab, sys, it2, t1s, t3s, rwa, rmin, resp):
    """ Calculates R3g response function
    
    
    """
    
#
#    For debugging, check if all arrays are fortran continuous
#
#    print(lab.orient_aver.flags['F_CONTIGUOUS'])
#    print(sys.Ns.flags['F_CONTIGUOUS'])
#    print(sys.om01.flags['F_CONTIGUOUS'])
#    print(sys.nn01.flags['F_CONTIGUOUS'])
#    print(sys.dd01.flags['F_CONTIGUOUS'])
#    print(sys.Kd01.flags['F_CONTIGUOUS'])
#    print(sys.Kd11.flags['F_CONTIGUOUS'])
#    print(sys.gofts.flags['F_CONTIGUOUS'])
#    print(sys.fptn.flags['F_CONTIGUOUS']) 
#    print(sys.SS1.flags['F_CONTIGUOUS'])
#    print(t1s.flags['F_CONTIGUOUS'])
#    print(t3s.flags['F_CONTIGUOUS'])
#    print(resp.flags['F_CONTIGUOUS'])
#    
#    it2_in = it2+1
#    
#    print(type(sys.Ns), sys.Ns.dtype)
    
    nr3td_fic.nr3_r3g_fic(lab.orient_aver, sys.Ns, sys.om01, sys.nn01,
                        sys.dd01, sys.Kd01, sys.Kd11, sys.gofts, sys.fptn, 
                        sys.SS1, it2+1, t1s, t3s, rwa, rmin, resp)



def nr3_r4g(lab, sys, it2, t1s, t3s, rwa, rmin, resp):
    """ Calculates R4g response function
    
    
    """
    
    
    nr3td_fic.nr3_r4g_fic(lab.orient_aver, sys.Ns, sys.om01, sys.nn01,
                        sys.dd01, sys.Kd01, sys.Kd11, sys.gofts, sys.fptn, 
                        sys.SS1, it2+1, t1s, t3s, rwa, rmin, resp)
    

def nr3_r1fs(lab, sys, it2, t1s, t3s, rwa, rmin, resp, plist=False):
    """ Calculates R1f* response function
    
    """

    if plist:
        nr3td_fic.nr3_r1fs_list_fic(lab.orient_aver, sys.Ns, sys.om01, 
                                    sys.om12, sys.nn01, sys.dd01, sys.nn12,
                                    sys.dd12, sys.Kd01, sys.Kd11, sys.Kd12,
                                    sys.gofts, sys.fptn, sys.SS1, sys.SS2,
                                    it2+1, t1s, t3s, rwa, rmin, resp)
        
    else:
        nr3td_fic.nr3_r1fs_fic(lab.orient_aver, sys.Ns, sys.om01, sys.om12, 
                           sys.nn01, sys.dd01, sys.nn12, sys.dd12, sys.Kd01,
                           sys.Kd11, sys.Kd12, sys.gofts, sys.fptn, 
                           sys.SS1, sys.SS2, it2+1, t1s, t3s, rwa, rmin, resp)
    

def nr3_r2fs(lab, sys, it2, t1s, t3s, rwa, rmin, resp, plist=False):
    """ Calculates R2f* response function
    
    """
    if plist:
        nr3td_fic.nr3_r2fs_list_fic(lab.orient_aver, sys.Ns, sys.om01,
                                    sys.om12, sys.nn01, sys.dd01, sys.nn12,
                                    sys.dd12, sys.Kd01, sys.Kd11, sys.Kd12,
                                    sys.gofts, sys.fptn, sys.SS1, sys.SS2,
                                    it2+1, t1s, t3s, rwa, rmin, resp)
    
    else:
        nr3td_fic.nr3_r2fs_fic(lab.orient_aver, sys.Ns, sys.om01, sys.om12, 
                           sys.nn01, sys.dd01, sys.nn12, sys.dd12, sys.Kd01,
                           sys.Kd11, sys.Kd12, sys.gofts, sys.fptn, 
                           sys.SS1, sys.SS2, it2+1, t1s, t3s, rwa, rmin, resp)


#
#
#  Energy transfer pathways
#
#
def nr3_r2g_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp):
    """ Calculates R2g response function
    
    

    Paramaters
    ----------

    lab : lab_settings
        Laboratory settings (polarizations of laser beams etc.) expressed 
        through the lab_setting class
        
    sys : band_system
        System to be calculated on expressed through the band_system class
        
    it2 : integer
        Index of the so-called waiting or population time of non-linear 
        spectroscopic techniques
        
    t1s : float array
        Values of t1 time for which response should be calculated
        
    t3s : float array
        Values of t3 time for which response shuld be calculated
        
    rwa: float
        Rotating wave frequency
        
    rmin: float
        Minimal value of the dipole prefactor, relative to its max value,
        which is taken into account
        
    resp : complex 2d array
        Non-linear response 
        
    """
    
    
    nr3td_fic.nr3_r2g_trans_fic(lab.orient_aver, sys.Ns, sys.om01, sys.nn01,
                        sys.dd01, sys.Kd01, sys.Kd11, sys.gofts, sys.fptn, 
                        sys.SS1, sys.Ueet2, it2+1, t1s, t3s, rwa, rmin, resp)

def nr3_r2g_trN(lab, sys, No, it2, t1s, t3s, rwa, rmin, resp):
    """ Calculates R2g response function with energy transfer 
    
    Calculates R2g response function with energy transfer to the No-th
    order of energy transfer with a complete treatment of the memory
    
    

    Paramaters
    ----------

    lab : lab_settings
        Laboratory settings (polarizations of laser beams etc.) expressed 
        through the lab_setting class
        
    sys : band_system
        System to be calculated on expressed through the band_system class
        
    No  : integer
        Maximum order of energy transfer to which we should calculate
        
    it2 : integer
        Index of the so-called waiting or population time of non-linear 
        spectroscopic techniques
        
    t1s : float array
        Values of t1 time for which response should be calculated
        
    t3s : float array
        Values of t3 time for which response shuld be calculated
        
    rwa: float
        Rotating wave frequency
        
    rmin: float
        Minimal value of the dipole prefactor, relative to its max value,
        which is taken into account
        
    resp : complex 2d array
        Non-linear response 
        
    """
    
    
    nr3td_fic.nr3_r2g_trn_fic(lab.orient_aver, sys.Ns, sys.om01, sys.nn01,
                        sys.dd01, sys.Kd01, sys.Kd11, sys.gofts, sys.fptn, 
                        sys.SS1, No, it2+1, t1s, t3s, rwa, rmin, resp)


def nr3_r1g_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp):
    """ Calculates R2g response function
    
    

    Paramaters
    ----------

    lab : lab_settings
        Laboratory settings (polarizations of laser beams etc.) expressed 
        through the lab_setting class
        
    sys : band_system
        System to be calculated on expressed through the band_system class
        
    it2 : integer
        Index of the so-called waiting or population time of non-linear 
        spectroscopic techniques
        
    t1s : float array
        Values of t1 time for which response should be calculated
        
    t3s : float array
        Values of t3 time for which response shuld be calculated
        
    rwa: float
        Rotating wave frequency
        
    rmin: float
        Minimal value of the dipole prefactor, relative to its max value,
        which is taken into account
        
    resp : complex 2d array
        Non-linear response 
        
    """
    
    
    nr3td_fic.nr3_r1g_trans_fic(lab.orient_aver, sys.Ns, sys.om01, sys.nn01,
                        sys.dd01, sys.Kd01, sys.Kd11, sys.gofts, sys.fptn, 
                        sys.SS1, sys.Ueet2, it2+1, t1s, t3s, rwa, rmin, resp)


def nr3_r1fs_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp):
    """ Calculates R1f* transfer response function
    
    """

    nr3td_fic.nr3_r1fs_trans_fic(lab.orient_aver, sys.Ns, sys.om01, sys.om12, 
                           sys.nn01, sys.dd01, sys.nn12, sys.dd12, sys.Kd01,
                           sys.Kd11, sys.Kd12, sys.gofts, sys.fptn, 
                           sys.SS1, sys.SS2, sys.Ueet2, it2+1, t1s, t3s, rwa,
                           rmin, resp)


def nr3_r2fs_trans(lab, sys, it2, t1s, t3s, rwa, rmin, resp):
    """ Calculates R1f* transfer response function
    
    """

    nr3td_fic.nr3_r2fs_trans_fic(lab.orient_aver, sys.Ns, sys.om01, sys.om12, 
                           sys.nn01, sys.dd01, sys.nn12, sys.dd12, sys.Kd01,
                           sys.Kd11, sys.Kd12, sys.gofts, sys.fptn, 
                           sys.SS1, sys.SS2, sys.Ueet2, it2+1, t1s, t3s, rwa,
                           rmin, resp)


