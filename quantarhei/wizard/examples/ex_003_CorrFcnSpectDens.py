# -*- coding: utf-8 -*-
"""
*******************************************************************************

    Demonstration of the correlation functions and spectral densities





*******************************************************************************
"""
from quantarhei import CorrelationFunction
from quantarhei import SpectralDensity
from quantarhei import TimeAxis
from quantarhei import energy_units

from quantarhei import DFunction
from quantarhei.core.units import kB_int 

import numpy

_show_plots_ = True

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
          "reorg":    20.0,
          "cortime":  100.0,
          "T":        temperature,
          "matsubara":20}
  
# here we create the correlation function assuming that the energy parameters
# are given in 1/cm        
with energy_units("1/cm"):
    cf = CorrelationFunction(ta,params)
    #cf.save_data("ob_20cm_100fs_300K_m20",ext="dat")

    sd = SpectralDensity(ta,params)
#sd.save_data("ob_sd_20cm_100fs_300K_m20", ext="dat")
    

if _show_plots_:
    # plotting the correlation function
    cf.plot(ylabel=r'$C(t)$ [rad$^2\cdot$fs$^{-2}$]',real_only=False)

# Fourier transform of the correlation function    
cF = cf.get_Fourier_transform()

# This is the same thing, just has the units management 
cF1 = cf.get_FTCorrelationFunction()
cF1o = cf.get_OddFTCorrelationFunction()
cF1e = cf.get_EvenFTCorrelationFunction()

if _show_plots_:
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

if _show_plots_:
    tf.plot(ylabel=r'$\tilde{C}(\omega)$ [rad$\cdot$fs$^{-1}$]',
            axis=[-0.1,0.1,-0.005,numpy.max(numpy.real(cF.data))*1.1],
            real_only=False)
    

 # Get spectral density       
sd = cf.get_SpectralDensity()

with energy_units("1/cm"):
    sd1 = SpectralDensity(ta,params)
    cf1 = CorrelationFunction(ta, params)
    cf1n = sd1.get_CorrelationFunction()
    sd1n = cf1.get_SpectralDensity()

if _show_plots_:
    cf1.plot(real_only=False,show=False)
    cf1n.plot()

lamb_cfm = cf1.measure_reorganization_energy()
lamb_def = cf1.get_reorganization_energy()
lamb_sdm = sd1.measure_reorganization_energy()
print(lamb_cfm, lamb_def, lamb_sdm)
print(numpy.allclose(lamb_cfm,lamb_def,rtol=1.0e-3))
print(numpy.allclose(lamb_sdm,lamb_def,rtol=1.0e-2))

if _show_plots_:
    sd1.plot(show=False)
    sd1n.plot(axis=[-0.1,0.1,-0.004,0.004])



    # Plot spectral density
    with energy_units("1/cm"):
        #cF1.plot(show=False)
        #cF1e.plot(show=False)
        cF1o.plot(show=False)
        sd1.plot(show=False)
        sd.plot(axis=[-1000,1000,-0.008,0.008])
    
df = numpy.max(numpy.abs(cF1o.data-sd1.data))
print(df)
mx = numpy.max(numpy.abs(sd1.data))
print(mx)
print(df/mx)    
           
with energy_units("eV"):
    tm = TimeAxis(0.0,100,1.0)
    wm = tm.get_FrequencyAxis()

with energy_units("eV"):
    print(wm.data[1]-wm.data[0])
    print(wm.step)
 
ftc = sd1.get_FTCorrelationFunction()  
if _show_plots_:
    cF1.plot(show=False)
    ftc.plot(axis=[-0.1,0.1,0,0.06])

 


    