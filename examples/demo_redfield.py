# -*- coding: utf-8 -*-

"""

    Propagation with Redfield theory using Aggregate object

"""
import numpy

from quantarhei import *
import quantarhei as qm

import matplotlib.pyplot as plt

from modelgenerator import ModelGenerator

print("""
*******************************************************************************
*                                                                             *
*                         Redfield Theory Demo                                *
*                                                                             *                  
*******************************************************************************
""")

Nt = 1000
dt = 2.0
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
prop_Redfield = agg.get_ReducedDensityMatrixPropagator(time,
                           relaxation_theory="standard_Redfield",
                           time_dependent=True, secular_relaxation=True)   

#
# Initial density matrix
#
shp = H.dim
rho_i1 = ReducedDensityMatrix(dim=shp, name="Initial DM")
rho_i1.data[shp-1,shp-1] = 1.0   
   
#
# Propagation of the density matrix
#   
rho_t1 = prop_Redfield.propagate(rho_i1,
                                 name="Redfield evolution from aggregate")
rho_t1.plot(coherences=False, axis=[0,Nt*dt,0,1.0], show=False)

#
# Thermal excited state to compare with
#
with eigenbasis_of(H):
    rho0 = agg.get_DensityMatrix(condition_type="thermal_excited_state",
                             relaxation_theory_limit="weak_coupling",
                             temperature=300)
 
#with eigenbasis_of(H):
if True:       
    pop = numpy.zeros((time.length,shp),dtype=numpy.float64)
    for i in range(1, H.dim):
        pop[:,i] = numpy.real(rho0.data[i,i]) 
        plt.plot(time.data,pop[:,i],'--k')

