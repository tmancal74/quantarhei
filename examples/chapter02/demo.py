# -*- coding: utf-8 -*-
from quantarhei import Hamiltonian
from quantarhei import StateVector
#from quantarhei import StateVectorEvolution
from quantarhei import StateVectorPropagator
from quantarhei import TimeAxis

#
# Hamiltonian of a 2 level molecule with adiabatic coupling
#

H = Hamiltonian(data=[[0.0, 0.4], [0.4, 1.0]])

psi = StateVector(data=[0.0, 1.0])

print(H)
print(psi)

psi2 = H.apply(psi)

print(psi2)

time = TimeAxis(0.0, 1000, 0.1)
prop = StateVectorPropagator(time, H)

psi_t = prop.propagate(psi)

print(psi_t)

psi_t.plot()

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

print(psi_vib_t)

psi_vib_t.plot()


"""
  Napsat evoluci systemu s 1. oscilatorem a vazbou
  
  .. potrebujeme napr. trace pres oscilator atd.


"""