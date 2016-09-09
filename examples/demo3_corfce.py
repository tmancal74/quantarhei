# -*- coding: utf-8 -*-

from quantarhei import CorrelationFunction
from quantarhei import TimeAxis
from quantarhei import energy_units

from matplotlib import pyplot as plt
import numpy


ta = TimeAxis(0.0,100,10.0)

params = {"ftype":    "OverdampedBrownian",
          "reorg":    20.0,
          "cortime":  100.0,
          "T":        100.0,
          "matsubara":20}
          
with energy_units("1/cm"):
    cf = CorrelationFunction(ta,params)


plt.plot(ta.data,numpy.real(cf.data),"-k")
plt.plot(ta.time,numpy.imag(cf.data),"-b")

plt.show()

with open("/home/tomas/c.dat","wt") as f:
    i = 0
    for t in ta.time:
        print(t,numpy.real(cf.data[i]), numpy.imag(cf.data[i]),file=f)
        i += 1
    
cf.convert_2_spectral_density()

plt.plot(ta.time,numpy.real(cf.data),"-k")
plt.plot(ta.time,numpy.imag(cf.data),"-b")

plt.show()
    
        