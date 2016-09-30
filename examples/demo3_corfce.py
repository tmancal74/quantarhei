# -*- coding: utf-8 -*-
"""
*******************************************************************************

    Demonstration of the correlation functions and spectral densities





*******************************************************************************
"""
from quantarhei import CorrelationFunction
from quantarhei import SpectralDensity
from quantarhei import TimeAxis, FrequencyAxis
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

temperature = 300.0
# parameters of the correlation function
params = {"ftype":    "OverdampedBrownian",
          "reorg":    30.0,
          "cortime":  200.0,
          "T":        temperature,
          "matsubara":20}
  
# here we create the correlation function assuming that the energy parameters
# are given in 1/cm        
with energy_units("1/cm"):
    cf = CorrelationFunction(ta,params)
    #cf.save("ob_20cm_100fs_100K_m20",ext="dat")

# plotting the correlation function
cf.plot(ylabel=r'$C(t)$ [rad$^2\cdot$fs$^{-2}$]',real_only=False)

# Fourier transform of the correlation function    
cF = cf.get_Fourier_transform()

# This is the same thing, just has the units management 
cF1 = cf.get_FTCorrelationFunction()
cF1o = cf.get_OddFTCorrelationFunction()
cF1e = cf.get_EvenFTCorrelationFunction()

# Plotting the Fourier transform
cF.plot(ylabel=r'$\tilde{C}(\omega)$ [rad$^2\cdot$fs$^{-1}$]',
        real_only=False,
        axis=[-0.4,0.4,-0.005,0.025],show=False)

        
wa = cF.axis

# Check the symmetry of the correlation function
k = 1999
vals = numpy.zeros(2000,dtype=numpy.complex)
l = 1
while k > 0:
    vals[l] = cF.data[k]*numpy.exp(-wa.data[k]/(kB_int*temperature))
    l += 1
    k -= 1
    
tf = DFunction(wa,vals)

tf.plot(ylabel=r'$\tilde{C}(\omega)$ [rad$\cdot$fs$^{-1}$]'
    ,axis=[-0.1,0.1,-0.005,numpy.max(numpy.real(cF.data))*1.1],real_only=False)
    

 # Get spectral density       
sd = cf.get_SpectralDensity()

with energy_units("1/cm"):
    sd1 = SpectralDensity(ta,params)
with energy_units("int"):
    sd1.save("ob_sd_30cm_200fs_300K_m20",ext="dat")

# Plot spectral density
with energy_units("1/cm"):
    #cF1.plot(show=False)
    #cF1e.plot(show=False)
    #cF1o.plot(show=False)
    sd1.plot(show=False)
    sd.plot(axis=[-1000,1000,-0.005,0.005]) 
           
with energy_units("eV"):
    tm = TimeAxis(0.0,100,1.0)
    wm = tm.get_FrequencyAxis()

with energy_units("eV"):
    print(wm.data[1]-wm.data[0])
    print(wm.step)

    