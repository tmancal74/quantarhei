# -*- coding: utf-8 -*-

from quantarhei import CorrelationFunction
from quantarhei import TimeAxis
from quantarhei import energy_units

from matplotlib import pyplot as plt
import numpy


ta = TimeAxis(0.0,1000,1.0)

params = {"ftype":    "OverdampedBrownian",
          "reorg":    20.0,
          "cortime":  100.0,
          "T":        300.0,
          "matsubara":20}
          
with energy_units("1/cm"):
    cf = CorrelationFunction(ta,params)

with energy_units("1/cm"):
    plt.plot(ta.data,numpy.real(cf.data),"-k")
    plt.plot(ta.time,numpy.imag(cf.data),"-b")

plt.show()
    
cf.convert_2_spectral_density()

plt.plot(ta.time,numpy.real(cf.data),"-k")
plt.plot(ta.time,numpy.imag(cf.data),"-b")

plt.show()
    