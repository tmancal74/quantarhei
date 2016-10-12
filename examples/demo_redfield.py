# -*- coding: utf-8 -*-

from quantarhei import *

time = TimeAxis(0.0, 10000, 0.1)

with energy_units("1/cm"):

    m1 = Molecule("Mol 1", [0.0, 10100.0])
    m2 = Molecule("Mol 2", [0.0, 10050.0])
    m3 = Molecule("Mol 3", [0.0, 10000.0])

    agg = Aggregate("Dimer")
    agg.add_Molecule(m1)
    agg.add_Molecule(m2)
    agg.add_Molecule(m3)
    agg.set_resonance_coupling(0,1,30.0)
    agg.set_resonance_coupling(1,2,30.0)

    params = dict(ftype="OverdampedBrownian", reorg=30, cortime=100, T=100)
    cf = CorrelationFunction(time, params)

m1.set_transition_environment((0,1), cf)
m2.set_transition_environment((0,1), cf)
m3.set_transition_environment((0,1), cf)


agg.build()


sbi = agg.get_SystemBathInteraction()
ham = agg.get_Hamiltonian()
print(ham)

#
# Calculation of various relaxation tensors
#
RRT = qm.RedfieldRelaxationTensor(ham, sbi)
print("\nRates from the full relaxation tensor")
for i in range(1, ham.dim):
    for j in range(1, ham.dim):
        print(i, "<-", j, ":", 1.0/numpy.real(RRT.data[i,i,j,j]))
        
RRM = qm.RedfieldRateMatrix(ham, sbi)
print("\nRates from the rate matrix")
for i in range(1,ham.dim):
    for j in range(1, ham.dim):
        print(i, "<-", j, ":", 1.0/RRM.data[i,j])
           
print("\nRatio")
for i in range(1, ham.dim):
    for j in range(1, ham.dim):
        print(i, "<-", j, ":", RRM.data[i,j]/numpy.real(RRT.data[i,i,j,j]))


TDRRM = qm.TDRedfieldRateMatrix(ham, sbi)
print("\nRates from the rate matrix")
for i in range(1,ham.dim):
    for j in range(1, ham.dim):
        print(i, "<-", j, ":", 1.0/TDRRM.data[time.length-1,i,j])
        
RRT2 = qm.RedfieldRelaxationTensor(ham, sbi)
print("\nRates from the full relaxation tensor")
for i in range(1, ham.dim):
    for j in range(1, ham.dim):
        print(i, "<-", j, ":", 1.0/numpy.real(RRT2.data[i,i,j,j]))

RRT.secularize()

#
# Evolution of reduced density matrix
#

prop = ReducedDensityMatrixPropagator(time, ham, RRT)

rho_i = ReducedDensityMatrix(dim=4)
with eigenbasis_of(ham):
    rho_i.data[3,3] = 1.0
   
print(rho_i)

print(1.0/RRT.data[1,1,2,2])
with eigenbasis_of(ham):
    print(1.0/RRT.data[1,1,2,2])
    rho_t = prop.propagate(rho_i)
    rho_t.plot(coherences=False)
    
#
#  Evolution of populations
#
prop = PopulationPropagator(time, RRM)
pop_t = prop.propagate([0.0, 0.0, 0.0, 1.0])

import matplotlib.pyplot as plt
plt.plot(time.data, pop_t[:,3],'--r')
plt.show()

print(RRM.data[2,3])
print(RRT.data[2,2,3,3])



