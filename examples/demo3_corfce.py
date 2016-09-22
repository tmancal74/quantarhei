# -*- coding: utf-8 -*-
"""
*******************************************************************************

    Demonstration of the correlation functions and spectral densities





*******************************************************************************
"""
from quantarhei import CorrelationFunction
from quantarhei import TimeAxis
from quantarhei import energy_units

from quantarhei import DFunction
from quantarhei.core.units import kB_int 

import numpy

"""

"""
print("\n")
print("*****************************************************************")
print("*                                                               *")
print("*    Demo of the quantarhei CorrelationFunction                 *")
print("*         and Spectral Density classes                          *")
print("*                                                               *")
print("*****************************************************************")

# time interval om which functions will be defined
ta = TimeAxis(0.0,1000,1.0)

# parameters of the correlation function
params = {"ftype":    "OverdampedBrownian",
          "reorg":    20.0,
          "cortime":  100.0,
          "T":        100.0,
          "matsubara":20}
  
# here we create the correlation function assuming that the energy parameters
# are given in 1/cm        
with energy_units("1/cm"):
    cf = CorrelationFunction(ta,params)

#FIXME It should be normal Fourier transform (but it gives a strange result)!!!
# plotting the correlation function
cf.plot(ylabel=r'$C(t)$ [rad$^2\cdot$fs$^{-2}$]',real_only=False)

# Fourier transform of the correlation function    
cF = cf.get_Fourier_transform()

# Plotting the Fourier transform
cF.plot(ylabel=r'$\tilde{C}(\omega)$ [rad$^2\cdot$fs$^{-1}$]',
        real_only=False,
        axis=[-0.4,0.4,-0.005,0.025],show=False)
        
wa = cF.axis

k = 1999
vals = numpy.zeros(2000,dtype=numpy.complex)
l = 1
while k > 0:
    vals[l] = cF.data[k]*numpy.exp(-wa.data[k]/(kB_int*100.0))
    l += 1
    k -= 1
    
tf = DFunction(wa,vals)

tf.plot(ylabel=r'$\tilde{C}(\omega)$ [rad$^2\cdot$fs$^{-1}$]'
    ,axis=[-0.1,0.1,-0.005,0.025],real_only=False)
    

        