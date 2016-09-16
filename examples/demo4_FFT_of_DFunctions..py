# -*- coding: utf-8 -*-
"""
*******************************************************************************

    Demonstration of the Fourier transforms with DFunctions





*******************************************************************************
"""

import numpy
import matplotlib.pyplot as plt

from quantarhei import TimeAxis
from quantarhei.core.dfunction import DFunction

print("\n")
print("*****************************************************************")
print("*                                                               *")
print("*    Fourier transform demo for quantarhei DFunction class      *")
print("*                                                               *")
print("*****************************************************************")
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

print(
"""    

***
        Fourier transform of a function defined
        on a symmetrical time interval          
***

    First we try the TimeAxis of the `complete` type. This represents 
    an arbitrary time interval and behaves as generally expected when
    the function (DFunction) defined on the interval is Fourier transformed.
    
    We Fourier transform a simple even exponential with 
    decay constant \gamma = 1/50 1/fs
""")

# symmetrically defined TimeAxis
t = TimeAxis(-(Ns//2)*dt,Ns,dt,atype="complete")

""" 

    Function no. 1


"""

# Function values
y = numpy.exp(-numpy.abs(t.time)*gg)

# we define DFunction 
f = DFunction(t,y)

# DFunction can return its Fourier transform
F = f.get_Fourier_transform()

# plot of the original function to be transformed
plt.plot(f.axis.time,f.data,"-b")
plt.axis([-1000,1000,0.0,1.1])
plt.title("Even exponential to be Fourier transformed")
plt.text(-800,0.8,r'$f(t) = e^{-\gamma |t|}$',fontdict={'size':20})
font={'size':20}
plt.xlabel('$t$ [fs]',**font)
plt.ylabel('$f(t)$',**font)
plt.show()

print(
"""    
    The result of the Fast Fourier transform is a Lorentzian function. 
    We compare the numerical values (in blue) with analytical result
    which is in green.
""")

# plot of the numerical Fourier tranform
F.plot(show=False,color="-b")
#plt.plot(F.axis.frequency,numpy.real(F.data),"-b")

# calculate the function analytically
k = 0
ff = numpy.zeros(F.axis.length)
for x in F.axis.frequency:
    ff[k] = 2.0*gg/(gg**2 + x**2)
    #print(x,ff[k])
    k += 1
    
# plot of the analytical transform (in the same figure)
fd = DFunction(F.axis,ff)
fd.plot(color="-g",axis=[-0.6,0.6,0.0,120],
        title = "Lorentzian as a result of the Fourier transform",
        text=[-0.5,80,r'$F(\omega) = \frac{2\gamma}{\gamma^2 + \omega^2}}}$'],
        text_font={'size':20})
# alternatively like below
#plt.plot(F.axis.frequency,ff,"-g")
#plt.axis([-0.6,0.6,0.0,120])
#plt.title("Lorentzian as a result of the Fourier transform")
#plt.text(-0.5,80,r'$F(\omega) = \frac{2\gamma}{\gamma^2 + \omega^2}}}$',
#         fontdict={'size':20})
#font={'size':20}
#plt.xlabel('$\omega$ [fs$^{-1}$]',**font)
#plt.ylabel('$F(\omega)$',**font)
#plt.show()

""" 
    
    Function no. 2 
    
    
"""

print(
"""    
    Next, we Fourier transform the same exponential with 
    and addition fast oscillation with frequency \omega_0 = 0.4 
""")
# Highly oscillating function; otherwise the same as before
y = numpy.exp(-numpy.abs(t.time)*gg-1j*om0*t.time)

f = DFunction(t,y)

plt.plot(f.axis.time,numpy.real(f.data))
plt.axis([-1000,1000,-1.1,1.1])
plt.title("Even exponential with fast oscillations")
plt.text(-900,0.8,r'$f(t) = e^{-\gamma |t| -i\omega_0 t}$',
         fontdict={'size':20})
font={'size':20}
plt.xlabel('$t$ [fs]',**font)
plt.ylabel('$f(t)$',**font)
plt.show()

print(
"""    
    Again, we can see an excellent agreement between analytical and
    numerical Fourier transform. Also the frequency of the function
    was correctly identified by the Fourier transform to be 0.4.
""")
F = f.get_Fourier_transform()
F.plot(show=False)
#plt.plot(F.axis.frequency,numpy.real(F.data))

k = 0
ff = numpy.zeros(F.axis.length)
for x in F.axis.frequency:
    ff[k] = 2.0*gg/(gg**2 + (x-om0)**2)
    k += 1
    
plt.plot(F.axis.frequency,numpy.real(ff))
plt.plot(F.axis.frequency,ff,"-g")
plt.axis([-0.6,0.6,0.0,120])
plt.title("Lorentzian as a result of the Fourier transform")
plt.text(-0.5,80,
         r'$F(\omega) = \frac{2\gamma}{\gamma^2 + (\omega-\omega_0)^2}}}$',
         fontdict={'size':20})
font={'size':20}
plt.xlabel('$\omega$ [fs$^{-1}$]',**font)
plt.ylabel('$F(\omega)$',**font)
plt.show()

print(
"""    
    Now we try to inverse Fourier transform our numerical result. 
    On a zoomed in figure, we can see an excellent agreement between
    the inverted function, and the original.
""")
f1 = F.get_inverse_Fourier_transform()

plt.plot(f1.axis.time,numpy.real(f1.data))
plt.plot(f.axis.time,numpy.real(f.data))
plt.axis([-200,200,-1,1.1])
plt.title("Reconstructed even exponential with fast oscillations")
plt.text(-180,0.8,r'$f(t) = e^{-\gamma |t| -i\omega_0 t}$',
         fontdict={'size':20})
font={'size':20}
plt.xlabel('$t$ [fs]',**font)
plt.ylabel('$f(t)$',**font)
plt.show()

print(
"""    
*** 
        Fourier transform of a function defined
        on the upper-half of a time interval    
***
""")

# by default atype="upper-half"
t = TimeAxis(0.0,Ns,dt) 

# Function values
y = numpy.exp(-numpy.abs(t.time)*gg)

# we define DFunction 
f = DFunction(t,y)

# plot of the original function to be transformed
f.plot(title="Positive time part of the even exponential",
       axis=[0,1000,0,1],
       text=[400,0.6,r'$f(t)=e^{-\gamma |t|}$'],
       text_font={"size":20})

# DFunction can return its Fourier transform
F = f.get_Fourier_transform()

# plot of the numerical Fourier tranform
F.plot(show=False,
       title="Lorentzian as a result of the Fourier transform",
       axis=[-0.6,0.6,0,120])

# calculate the function analytically
k = 0
ff = numpy.zeros(F.axis.length)
for x in F.axis.frequency:
    ff[k] = 2.0*gg/(gg**2 + x**2)
    k += 1
    
# plot of the analytical transform (in the same figure)
df = DFunction(F.axis,ff)
df.plot()

# alternatively one can do it without creating a DFunction
#plt.plot(F.axis.frequency,ff)
#plt.show()
