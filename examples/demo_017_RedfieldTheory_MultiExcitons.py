# -*- coding: utf-8 -*-

"""

    Propagation with Redfield theory using Aggregate object

"""

_show_plots_ = False

import time

import numpy

import quantarhei as qr
import quantarhei.models as models
import matplotlib.pyplot as plt

print("""
*******************************************************************************
*                                                                             *
*           Multi-exciton Redfield Theory Demo                                *
*                                                                             *
*                                                                             *                  
*******************************************************************************
""")

Nt = 1000
dt = 2.0
timea = qr.TimeAxis(0.0, Nt, dt)

mg = models.ModelGenerator()
print("Generating model aggregate:")
agg = mg.get_Aggregate_with_environment(name="pentamer-1_env",
                                        timeaxis=timea)
print("... done")

print("Building aggregate internals:")
t1 = time.time()
agg.build()
t2 = time.time()
print("...done in", t2-t1,"sec")

H = agg.get_Hamiltonian()
#with energy_units("1/cm"):
#    print(H)

#
# Aggregate object can return a propagator
#
print("Setting up Redfield propagator:")
t1 = time.time()
prop_Redfield = agg.get_ReducedDensityMatrixPropagator(timea,
                           relaxation_theory="standard_Redfield",
                           time_dependent=False, secular_relaxation=True)   
t2 = time.time()
print("...done in", t2-t1, "sec")

#
# Initial density matrix
#
print("Initial condition for the density matrix")
shp = H.dim
rho_i1 = qr.ReducedDensityMatrix(dim=shp, name="Initial DM")

#
# Initial condition should be set to 2-exciton band
#
rho_i1.data[shp-1,shp-1] = 1.0   
print("...done")
   
#
# Propagation of the density matrix
#
print("Propagating density matrix:")  
t1 = time.time() 
rho_t1 = prop_Redfield.propagate(rho_i1,
                                 name="Redfield evolution from aggregate")
t2 = time.time()
print("...done in", t2-t1, "sec")

if _show_plots_: 
    
    rho_t1.plot(coherences=False, axis=[0,Nt*dt,0,1.0], show=False)

    #
    # Thermal excited state to compare with
    #
    with qr.eigenbasis_of(H):
        rho0 = agg.get_DensityMatrix(condition_type="thermal_excited_state",
                                     relaxation_theory_limit="weak_coupling",
                                     temperature=300)
      
    pop = numpy.zeros((timea.length,shp),dtype=numpy.float64)
    for i in range(1, H.dim):
        pop[:,i] = numpy.real(rho0.data[i,i]) 
        plt.plot(timea.data,pop[:,i],'--k')

