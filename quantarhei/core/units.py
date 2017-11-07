# -*- coding: utf-8 -*-

import scipy.constants as const

# frequency
conversion_facs_frequency = {
    "int"    : 1.0,
    "1/fs" : 1.0, 
    "1/cm"   : 2.0*const.pi*const.c*1.0e-13, 
    "THz"    : 2.0*const.pi*1.0e-03,
    "Hz"     : 2.0*const.pi,
    "SI"     : 2.0*const.pi,
    "nm"     : 1.0/(1.0e7*2.0*const.pi*const.c*1.0e-13),
    "a.u."   : 27.21138602*1.0e-15*const.e/const.hbar,
    "Ha"     : 27.21138602*1.0e-15*const.e/const.hbar      
    } 

# energy
conversion_facs_energy = {
    "int"    : 1.0,
    "1/fs" : 1.0, 
    "1/cm"   : 2.0*const.pi*const.c*1.0e-13, 
    "THz"    : 2.0*const.pi*1.0e-03,
    "eV"     : 1.0e-15*const.e/const.hbar,
    "meV"    : 1.0e-18*const.e/const.hbar,
    "J"      : 1.0e-15/const.hbar,
    "SI"     : 1.0e-15/const.hbar,
    "nm"     : 1.0/(1.0e7*2.0*const.pi*const.c*1.0e-13),
    "a.u."   : 27.21138602*1.0e-15*const.e/const.hbar,
    "Ha"     : 27.21138602*1.0e-15*const.e/const.hbar    
    } 


#conversion_facs_position = {
#    "int" : 1.0,
#    "A"   : 1.0,
#    "m"   : 1.0e10,
#    "SI"  : 1.0e10
#}


conversion_facs_time = {
    "fs" : 1.0,
    "ps" : 1000.0,
    "ns" : 1.0e6,
    "s"  : 1.0e15,
    "SI" : 1.0e15
}


conversion_facs_temperature = {
    "K" : 1.0,
    "C" : 1.0,
    "F" : 0.0
}

conversion_offs_temperature = {
    "K" : 0.0,
    "C" : 273.15,
    "F" : 0.0
}

conversion_facs_edipole = {
    "int": 1.0,
    "au" : 1.0,
    "D"  : 1.0/0.20819434,
    "Cm" : 1.0e-21/const.c,
    "SI" : 1.0e-21/const.c
}

conversion_facs_length = {
    "int": 1.0,
    "A" : 1.0,
    "Bohr"  : 0.52917721067,
    "a.u." : 0.52917721067,
    "nm" : 10.0,
    "m"  : 1.0e10,
    "SI" : 1.0e10
}

#
# Some useful factors (deprecated)
#
cm2int = conversion_facs_energy["1/cm"]
int2cm = 1.0/cm2int
J2int = conversion_facs_energy["J"]

#
# Important constants in internal units
#

# Boltzmann constant
kB_SI = const.k
kB_int = const.k*conversion_facs_energy["J"]

kB_cmK = 0.69503476      # deprecated
kB_intK = kB_cmK*cm2int  # deprecated

# Planck constant
hbar_int = 1.0

# Speed of light (A/fs)
c_SI = const.c
c_int = const.c*1.0e-5

# Permeability of vacuum

eps0_int = 1.0e19/(4.0*const.pi*J2int)



#
# Conversion function
#


def convert(val, in_units, to=None):
    """Converts value in certain units into other units
    
    """
    from .managers import energy_units
    from .managers import Manager
    m = Manager()
    with energy_units(in_units):
        e = m.convert_energy_2_internal_u(val)
    
    if to is None:
        ne = m.convert_energy_2_current_u(e)
    else:
        with energy_units(to):
            ne = m.convert_energy_2_current_u(e)
        
    return ne

def in_current_units(val, in_units):
    """Converts value in certain units into the current units
    
    """
    from .managers import energy_units
    from .managers import Manager

    m = Manager()
    with energy_units(in_units):
        e = m.convert_energy_2_internal_u(val)
    
    return m.convert_energy_2_current_u(e)
         


