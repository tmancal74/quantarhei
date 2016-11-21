# -*- coding: utf-8 -*-
"""

    Propagation with Foerster theory using Aggregate object




"""


import numpy

from quantarhei import *
import quantarhei as qm

import matplotlib.pyplot as plt

from modelgenerator import ModelGenerator

print("""
*******************************************************************************
*                                                                             *
*                         Foerster Theory Demo                                *
*                                                                             *                  
*******************************************************************************
""")

Nt = 1000
dt = 1.0
time = TimeAxis(0.0, Nt, dt)

mg = ModelGenerator()
agg = mg.get_Aggregate_with_environment(name="pentamer-1_env",
                                        timeaxis=time)

agg.build()



H = agg.get_Hamiltonian()
#with energy_units("1/cm"):
#    print(H)

#
# Aggregate object can return a propagator
#
prop_Foerster = agg.get_ReducedDensityMatrixPropagator(time,
                           relaxation_theory="standard_Foerster",
                           time_dependent=False)   

#
# Initial density matrix
#
shp = H.dim
rho_i1 = ReducedDensityMatrix(dim=shp, name="Initial DM")
rho_i1.data[shp-1,shp-1] = 1.0   
   
#
# Propagation of the density matrix
#   
with eigenbasis_of(H):
    rho_t1 = prop_Foerster.propagate(rho_i1,
                                 name="Foerster evolution from aggregate")
    
rho_t1.plot(coherences=False, axis=[0,Nt*dt,0,1.0])

#
# Thermal excited state to compare with
#
rho0 = agg.get_DensityMatrix(condition_type="thermal_excited_state",
                             relaxation_theory_limit="strong_coupling",
                             temperature=300)
 
with eigenbasis_of(H):
#if True:       
    pop = numpy.zeros((time.length,shp),dtype=numpy.float64)
    for i in range(1, H.dim):
        pop[:,i] = numpy.real(rho0.data[i,i]) 
        plt.plot(time.data,pop[:,i],'--k')

#    # plot the termal distrubution
#    plt.plot(time.data,pop[:,1],'--r')
#    plt.plot(time.data,pop[:,2],'--b')
#    plt.plot(time.data,pop[:,3],'--g')  

RR = prop_Foerster.RelaxationTensor
with eigenbasis_of(H):
    if isinstance(RR, qm.core.time.TimeDependent):
        print(RR.data[:,1,1,2,2])
