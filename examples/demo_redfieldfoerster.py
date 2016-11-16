# -*- coding: utf-8 -*-
"""

    Propagation with Redfield-Foerster theory using Aggregate object

"""


import numpy

from quantarhei import *
import quantarhei as qm

import matplotlib.pyplot as plt

print("""
*******************************************************************************
*                                                                             *
*                         Refield-Foerster Theory Demo                        *
*                                                                             *                  
*******************************************************************************
""")
Nt = 5000
dt = 1.0
time = TimeAxis(0.0, Nt, dt)
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
    
#    m4 = Molecule("Mol 4", [0.0, 11000.0])
#    m5 = Molecule("Mol 5", [0.0, 11000.0])   
#    m4.position = [15.0, 15.0, 0.0]
#    m5.position = [15.0, 10.0, 0.0]
#    m4.set_dipole(0,1,[5.8, 0.0, 0.0])
#    m5.set_dipole(0,1,[5.8, 0.0, 0.0])
#    agg.add_Molecule(m4)
#    agg.add_Molecule(m5)
    
    agg.set_coupling_by_dipole_dipole()

    params = dict(ftype="OverdampedBrownian", reorg=20, cortime=100, T=300)
    cf = CorrelationFunction(time, params)

m1.set_transition_environment((0,1), cf)
m2.set_transition_environment((0,1), cf)
m3.set_transition_environment((0,1), cf)

#m4.set_transition_environment((0,1), cf)
#m5.set_transition_environment((0,1), cf)

agg.build()

H = agg.get_Hamiltonian()
with energy_units("1/cm"):
    print(H)
   
cutoff = 130.0
#
# Aggregate object can return a propagator
#
with energy_units("1/cm"):
    prop_comb = agg.get_ReducedDensityMatrixPropagator(time,
                           relaxation_theory="combined_RedfieldFoerster",
                           coupling_cutoff=cutoff,
                           secular_relaxation=False)   

#
# Initial density matrix
#
shp = 4
rho_i1 = ReducedDensityMatrix(dim=shp, name="Initial DM")
rho_i1.data[shp-2,shp-2] = 1.0   
   
#
# Propagation of the density matrix
#   
rho_t1 = prop_comb.propagate(rho_i1,
                             name="Combined evolution from aggregate")
    
with energy_units("1/cm"):
    H.remove_cutoff_coupling(cutoff) 
    print(H) 
    
#with eigenbasis_of(H):
if True:
    RR = prop_comb.RelaxationTensor
    for aa in range(RR.dim):
        rsum = 0.0
        for bb in range(RR.dim):
            print(aa,bb,RR.data[aa,aa,bb,bb])
            rsum += RR.data[bb,bb,aa,aa]
        print(aa, " sum ", numpy.real(rsum))        
         
                   
with eigenbasis_of(H):
#if True:
    #rho_eq = agg.get_DensityMatrix(condition_type="")
    Nshow = Nt*dt
    rho_t1.plot(coherences=False, axis=[0,Nshow,0,1.0])