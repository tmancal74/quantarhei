# -*- coding: utf-8 -*-

import numpy
import matplotlib.pyplot as plt

from quantarhei import TimeAxis
from quantarhei.core.dfunction import DFunction

"""
*******************************************************************************

    Demonstration of the Fourier transforms with DFunctions





*******************************************************************************
"""

"""

    Parameters of Lorenzian functions which serve as tests of our ability to
    correctly Fourier transform numerically. We use
    
    e^(-|gg|*t)
    
    and
    
    e^(-|gg|*t -i*om0*t)
    
"""
gg = 1.0/50.0
om0 = 0.4

"""

    Numerics is done with a step `dt` abd `Ns` steps
    
"""
dt = 5.0
Ns = 1000


"""
    First we try the TimeAxis of the `complete` type. This represents 
    an arbitrary time interval and behaves as generally expected when
    the function (DFunction) defined on the interval is Fourier transformed.

"""
# symmetrically defined TimeAxis
t = TimeAxis(-(Ns//2)*dt,Ns,dt,atype="complete")

""" Function no. 1 """

# Function values
y = numpy.exp(-numpy.abs(t.time)*gg)

# we define DFunction 
f = DFunction(t,y)

# DFunction can return its Fourier transform
F = f.get_Fourier_transform()

# plot of the original function to be transformed
plt.plot(f.axis.time,f.data)
plt.show()

# plot of the numerical Fourier tranform
plt.plot(F.axis.frequency,numpy.real(F.data))

# calculate the function analytically
k = 0
ff = numpy.zeros(F.axis.length)
for x in F.axis.frequency:
    ff[k] = 2.0*gg/(gg**2 + x**2)
    #print(x,ff[k])
    k += 1
    
# plot of the analytical transform (in the same figure)
plt.plot(F.axis.frequency,ff)
plt.show()

""" Function no. 2 """

# Highly oscillating function; otherwise the same as before
y = numpy.exp(-numpy.abs(t.time)*gg-1j*om0*t.time)

f = DFunction(t,y)
F = f.get_Fourier_transform()

plt.plot(f.axis.time,numpy.real(f.data))
plt.show()

plt.plot(F.axis.frequency,numpy.real(F.data))
k = 0
ff = numpy.zeros(F.axis.length)
for x in F.axis.frequency:
    ff[k] = 2.0*gg/(gg**2 + (x-om0)**2)
    #print(x,ff[k])
    k += 1
    
plt.plot(F.axis.frequency,numpy.real(ff))
plt.show()

"""

    Demonstration of the upper-half type of TimeAxis. 

"""

""" by default atype="upper-half" """
t = TimeAxis(0.0,Ns,dt) 

# Function values
y = numpy.exp(-numpy.abs(t.time)*gg)

# we define DFunction 
f = DFunction(t,y)

# DFunction can return its Fourier transform
F = f.get_Fourier_transform()

# plot of the original function to be transformed
plt.plot(f.axis.time,f.data)
plt.show()

# plot of the numerical Fourier tranform
plt.plot(F.axis.frequency,numpy.real(F.data))

# calculate the function analytically
k = 0
ff = numpy.zeros(F.axis.length)
for x in F.axis.frequency:
    ff[k] = 2.0*gg/(gg**2 + x**2)
    #print(x,ff[k])
    k += 1
    
# plot of the analytical transform (in the same figure)
plt.plot(F.axis.frequency,ff)
plt.show()