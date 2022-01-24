# -*- coding: utf-8 -*-
"""

    Propagation with Redfield-Foerster theory using Aggregate object

"""


import numpy

from quantarhei import *
import quantarhei as qm

import matplotlib.pyplot as plt

from modelgenerator import ModelGenerator

print("""
*******************************************************************************
*                                                                             *
*                         Refield-Foerster Theory Demo                        *
*                                                                             *                  
*******************************************************************************
""")
Nt = 1000
dt = 1.0
time = TimeAxis(0.0, Nt, dt)


mg = ModelGenerator()
agg = mg.get_Aggregate_with_environment(name="pentamer-1_env", timeaxis=time)

agg.build()

H = agg.get_Hamiltonian()
#with energy_units("1/cm"):
#    print(H)
   
cutoff = 0.0
#
# Aggregate object can return a propagator
#
#Manager().warn_about_basis_changing_objects = True
with energy_units("1/cm"):
    print("Calculating 1")
    prop_comb = agg.get_ReducedDensityMatrixPropagator(time,
                           relaxation_theory="standard_Redfield",
                           coupling_cutoff=cutoff,
                           secular_relaxation=False,
                           time_dependent=True)
    
    print("Calculating 2")    
    prop_comb2 = agg.get_ReducedDensityMatrixPropagator(time,
                           relaxation_theory="combined_RedfieldFoerster",
                           coupling_cutoff=cutoff,
                           secular_relaxation=False,
                           time_dependent=True)
    

#
# Initial density matrix
#
shp = H.dim
rho_i1 = ReducedDensityMatrix(dim=shp, name="Initial DM")
with eigenbasis_of(H):
    rho_i1.data[shp-1,shp-1] = 1.0   
   

# Propagation of the density matrix
#   
    rho_t1 = prop_comb.propagate(rho_i1,
                             name="Combined evolution from aggregate")

    rho_t2 = prop_comb2.propagate(rho_i1,
                             name="Combined evolution from aggregate")
    
#with energy_units("1/cm"):
#    H.remove_cutoff_coupling(cutoff) 
#    print(H) 

with eigenbasis_of(H):
#if True:
    RR = prop_comb.RelaxationTensor
    RR.name = "stR"
    print(RR.data[10,1,1,2,2], id(RR))
#if True:
    RT = prop_comb2.RelaxationTensor
    RT.name = "cRF"
    print(RT.data[10,1,1,2,2], id(RT))

#    for aa in range(RR.dim):
#        rsum = 0.0
#        for bb in range(RR.dim):
#            print(aa,bb,RR.data[aa,aa,bb,bb])
#            rsum += RR.data[bb,bb,aa,aa]
#        print(aa, " sum ", numpy.real(rsum))        
         
                 
#with eigenbasis_of(H):
if True:
    #rho_eq = agg.get_DensityMatrix(condition_type="")
    Nshow = Nt*dt
    rho_t1.plot(coherences=False, axis=[0,Nshow,0,1.0])
    rho_t2.plot(coherences=False, axis=[0,Nshow,0,1.0])