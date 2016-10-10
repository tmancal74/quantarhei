# -*- coding: utf-8 -*-
from quantarhei import Hamiltonian
from quantarhei import StateVector
#from quantarhei import StateVectorEvolution
from quantarhei import StateVectorPropagator
from quantarhei import TimeAxis

H = Hamiltonian(data=[[0.0, 0.4], [0.4, 1.0]])

psi = StateVector(data=[0.0, 1.0])

print(H)
print(psi)

psi2 = H.apply(psi)

print(psi2)

time = TimeAxis(0.0, 1000, 0.1)
prop = StateVectorPropagator(time, H)

psi_t = prop.propagate(psi)

psi_t.plot()

print(psi_t)

