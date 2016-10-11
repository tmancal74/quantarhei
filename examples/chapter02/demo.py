# -*- coding: utf-8 -*-

import numpy

from quantarhei import Hamiltonian
from quantarhei import StateVector
from quantarhei import StateVectorPropagator
from quantarhei import TimeAxis

from quantarhei import eigenbasis_of

#
# Hamiltonian of a 2 level molecule with adiabatic coupling
#

H = Hamiltonian(data=[[0.0, 0.4], [0.4, 1.0]])

psi = StateVector(data=[0.0, 1.0])

print(H)
print(psi)

time = TimeAxis(0.0, 1000, 0.1)

prop = StateVectorPropagator(time, H)

psi_t = prop.propagate(psi)

with eigenbasis_of(H):
    psi_t.plot(ptype="square")
psi_t.plot(ptype="square")


#
# The same Hamiltonian using the Molecule class
#

from quantarhei import Molecule

m = Molecule("Mol 1", [0.0, 1.0])
m.set_adiabatic_coupling(0,1,0.4)

from quantarhei import Mode

vib1 = Mode(frequency=0.01)
m.add_Mode(vib1)

vib1.set_shift(1, 0.5)


Hm = m.get_Hamiltonian()
print(Hm)

print(m)

psi_vib = StateVector(4)
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

mol1 = Molecule("Mol 1", [0.0, 1.0])
mol2 = Molecule("Mol 2", [0.0, 1.0])

agg = Aggregate("Dimer")

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
