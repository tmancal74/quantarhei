# -*- coding: utf-8 -*-
"""

    Propagation with Foerster theory using Aggregate object




"""


import numpy

from quantarhei import *
import quantarhei as qm

print("""
*******************************************************************************
*                                                                             *
*                         Redfield Theory Demo                                *
*                                                                             *                  
*******************************************************************************
""")

time = TimeAxis(0.0, 4000, 1.0)
with energy_units("1/cm"):

    m1 = Molecule("Mol 1", [0.0, 10100.0])
    m2 = Molecule("Mol 2", [0.0, 10050.0])
    m3 = Molecule("Mol 3", [0.0, 10000.0])
    
    m1.position = [0.0, 0.0, 0.0]
    m2.position = [15.0, 0.0, 0.0]
    m3.position = [10.0, 10.0, 0.0]
    m1.set_dipole(0,1,[5.8, 0.0, 0.0])
    m2.set_dipole(0,1,[5.8, 0.0, 0.0])
    m3.set_dipole(0,1,[numpy.sqrt(12.0), 0.0, 0.0])

    agg = Aggregate("Trimer")
    agg.add_Molecule(m1)
    agg.add_Molecule(m2)
    agg.add_Molecule(m3)
    
    agg.set_coupling_by_dipole_dipole(epsr=1.92921)

    params = dict(ftype="OverdampedBrownian", reorg=20, cortime=100, T=100)
    cf = CorrelationFunction(time, params)

m1.set_transition_environment((0,1), cf)
m2.set_transition_environment((0,1), cf)
m3.set_transition_environment((0,1), cf)


agg.build()

H = agg.get_Hamiltonian()
with energy_units("1/cm"):
    print(H)

#
# Aggregate object can return a propagator
#
prop_Redfield = agg.get_ReducedDensityMatrixPropagator(time,
                           relaxation_theory="standard_Redfield")
prop_Foerster = agg.get_ReducedDensityMatrixPropagator(time,
                           relaxation_theory="standard_Foerster") 
with energy_units("1/cm"):
    prop_RF = agg.get_ReducedDensityMatrixPropagator(time,
                           relaxation_theory="combined_RedfieldFoerster",
                           coupling_cutoff=50*qm.core.units.cm2int)
                           
RF = prop_Foerster.RelaxationTensor
RT = prop_Redfield.RelaxationTensor
#RF = prop_RF.RelaxationTensor

print(1.0/RF.data[1,1,2,2], RT.data[1,1,2,2])
print(1.0/RF.data[3,3,1,1])
print(1.0/RF.data[1,1,3,3])
print(1.0/RF.data[2,2,1,1])
print(1.0/RF.data[3,3,2,2])
print(1.0/RF.data[2,2,3,3])
print(1.0/RF.data[1,1,1,1])
print(1.0/RF.data[2,2,2,2])
print(1.0/RF.data[3,3,3,3])        


rho_i1 = ReducedDensityMatrix(dim=4, name="Initial DM")
rho_i1.data[1,1] = 1.0   

with eigenbasis_of(H): 
#if True:   
    rho_t1 = prop_Foerster.propagate(rho_i1,
                                     name="Foerster evolution from aggregate")
                 
#with eigenbasis_of(H):
if True:
    rho_t2 = prop_Redfield.propagate(rho_i1,
                                     name="Redfield evolution from aggregate")
#    rho_t3 = prop_RF.propagate(rho_i1)

#with eigenbasis_of(H):
if True:
    rho_t1.plot(coherences=False)    
    rho_t2.plot(coherences=False)
#rho_t3.plot(coherences=False)
    
