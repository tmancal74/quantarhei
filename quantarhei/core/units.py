# -*- coding: utf-8 -*-

import scipy.constants as const

# frequency/energy
conversion_facs_frequency = {
    "2pi/fs" : 1.0, 
    "1/cm"   : 2.0*const.pi*const.c*1.0e-13, 
    "Thz"    : 100
    } 

conversion_facs_position = {}

conversion_facs_time = {}

#
# Some useful factors (deprecated)
#
cm2int = conversion_facs_frequency["1/cm"]
int2cm = 1.0/cm2int


#
# Important constants in internal units
#

# Boltzmann constant
kB_cmK = 0.69503476      # deprecated
kB_intK = kB_cmK*cm2int  # deprecated

# Planck constant
hbar_int = 1.0

# Speed of light
c_int = 1.0

