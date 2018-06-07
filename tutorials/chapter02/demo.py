# -*- coding: utf-8 -*-
print(
"""

       Quantarhei Demo for the Chapter 2 of the Quantarhei Theory Text
""")
import numpy
#import matplotlib
#matplotlib.use("wx")
#import matplotlib.pyplot as plt
#plt.ion()

#
# Quantarhei objects needed below
#
from quantarhei import TimeAxis
from quantarhei import Hamiltonian
from quantarhei import StateVector
from quantarhei import StateVectorPropagator
    
#
# Context manager for basis management in Quantarhei
#
from quantarhei import eigenbasis_of

def pause(text=None):
    if text is not None:
        print(text)
    input('Press <ENTER> to continue')
    
pause()
print("""
###############################################################################
#                                                                             #
#    EXAMPLE 2.1                                                              #
#                                                                             #
###############################################################################
""")
#
# Hamiltonian of a 2 level molecule with adiabatic coupling
# Default energy units are radians per femtosecond (rad/fs)
#
print("""
We define a simple 2x2 Hamiltonian and an initial state vector with the
excited state completely populated.
""")

h = [[0.0, 0.4], 
     [0.4, 1.0]] 
H = Hamiltonian(data=h)

#
# Initial state vector to be propagated by the Hamiltonian above
#
psi_0 = StateVector(data=[0.0, 1.0])

#
# Check the Hamiltonian and state vector
#
print("Here are the two objects: ")
print(H)
print(psi_0)


pause("\nWe will proceed to calculate time evolution ...")


print("""
Now we calculate the evolution of the system from the initial state. Read the
code to see how this is done.
""")
#
# Time interval on which things will be calculated 
#
# TimeAxis(start, number of steps, stepsize)
#
time = TimeAxis(0.0, 1000, 0.1)

#
# Using time axis and the Hamiltonian we define a propagator
# for the state vector
#
prop = StateVectorPropagator(time, H)

#
# This is where time evolution of the state vector is calculated
# It is stored in psi_t which is of StateVectorEvolution type
#
psi_t = prop.propagate(psi_0)


print(
"""We plot the time dependence of the state vector elements on the plots below:

First we plot the squares of the absolute values of the state vector elements
in the eigenstate basis of the Hamiltonian. In other words, we plot
the probabilites of finding the system in its eigenstates.

In the same plot, we plot the real parts of the same elements.
""")
with eigenbasis_of(H):
    psi_t.plot(ptype="square", show=False)
    psi_t.plot(ptype="real")
    
    
pause("\nNext we plot the same in the original basis...")


print("""
Then we plot the same thing but in the basis of states in which defined the
values.
""")
psi_t.plot(ptype="square", show=False)
psi_t.plot(ptype="real")


pause("\nLet us add some vibrational states to the studied system ...")

print("""
###############################################################################
#                                                                             #
#    EXAMPLE 2.2                                                              #
#                                                                             #
###############################################################################
""")
#
# The same Hamiltonian using the Molecule class
#

from quantarhei import Molecule

m = Molecule(name="Mol 1", elenergies=[0.0, 1.0])
m.set_adiabatic_coupling(0,1,0.4)

from quantarhei import Mode

vib1 = Mode(frequency=0.01)
m.add_Mode(vib1)

vib1.set_shift(1, 0.5)
vib1.set_nmax(0, 5)
vib1.set_nmax(1, 5)


Hm = m.get_Hamiltonian()
print(Hm)
print(m)

psi_vib = StateVector(10)
psi_vib.data[3] = 1.0

prop_vib = StateVectorPropagator(time, Hm)

psi_vib_t = prop_vib.propagate(psi_vib)

with eigenbasis_of(Hm):
    psi_vib_t.plot(ptype="square")

psi_vib_t.plot(ptype="square")

sm = numpy.zeros(time.length)
for i in range(time.length):
    sm[i] = numpy.sum(numpy.abs(psi_vib_t.data[i,:])**2)

import matplotlib.pyplot as plt

plt.plot(time.data, sm)
plt.axis([0.0,100.0, 0.0, 1.1])
plt.show()


#
# Molecular dimer without vibrations
#
from quantarhei import Aggregate

mol1 = Molecule(name="Mol 1", elenergies=[0.0, 1.0])
mol2 = Molecule(name="Mol 2", elenergies=[0.0, 1.0])

agg = Aggregate(name="Dimer")

agg.add_Molecule(mol1)
agg.add_Molecule(mol2)
agg.set_resonance_coupling(0,1,0.01)

agg.build()

H = agg.get_Hamiltonian()

print(H)

psi = StateVector(3)
psi.data[2] = 1.0

dimer_propagator = StateVectorPropagator(time, H)
psi_t = dimer_propagator.propagate(psi)
psi_t.plot(ptype="square", show=False)

#dat = [[0.0, 0.0, 0.0], [0.0, 1.0, 0.01], [0.0, 0.01, 1.1]]
H._data[2,2] = 1.2 #dat

print(H.data)

dimer_propagator = StateVectorPropagator(time, H)
psi_t = dimer_propagator.propagate(psi)
psi_t.plot(ptype="square")
