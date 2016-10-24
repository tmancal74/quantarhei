# -*- coding: utf-8 -*-
import numpy

from quantarhei import *

print("""
*******************************************************************************
*                                                                             *
*                         Redfield Theory Demo                                *
*                                                                             *                  
*******************************************************************************
""")

time = TimeAxis(0.0, 100000, 0.1)

with energy_units("1/cm"):

    m1 = Molecule("Mol 1", [0.0, 10100.0])
    m2 = Molecule("Mol 2", [0.0, 10050.0])
    m3 = Molecule("Mol 3", [0.0, 10000.0])
    
    m1.position = [0.0, 0.0, 0.0]
    m2.position = [10.0, 0.0, 0.0]
    m3.position = [5.0, 5.0, 0.0]
    m1.set_dipole(0,1,[numpy.sqrt(12.0), 0.0, 0.0])
    m2.set_dipole(0,1,[numpy.sqrt(12.0), 0.0, 0.0])
    m3.set_dipole(0,1,[numpy.sqrt(12.0), 0.0, 0.0])

    agg = Aggregate("Dimer")
    agg.add_Molecule(m1)
    agg.add_Molecule(m2)
    agg.add_Molecule(m3)
    
    #agg.set_resonance_coupling(0,1,30.0)
    #agg.set_resonance_coupling(1,2,30.0)
    
    agg.set_coupling_by_dipole_dipole(epsr=1.0)

    params = dict(ftype="OverdampedBrownian", reorg=20, cortime=100, T=100)
    cf = CorrelationFunction(time, params)

m1.set_transition_environment((0,1), cf)
m2.set_transition_environment((0,1), cf)
m3.set_transition_environment((0,1), cf)


agg.build()


#
# Aggregate object can return a propagator
#
prop1 = agg.get_ReducedDensityMatrixPropagator(time,
                           relaxation_theory="standard_Redfield",
                           time_dependent=False,
                           secular_relaxation=True)

prop_nonsec = agg.get_ReducedDensityMatrixPropagator(time,
                           relaxation_theory="standard_Redfield")



sbi = agg.get_SystemBathInteraction()
ham = agg.get_Hamiltonian()
ham.set_name("Hamiltonian")

print(">>> Hamiltonian ")
with energy_units("1/cm"):
    print(ham)

raise Exception()

print("""
*******************************************************************************

                    Calculating relaxation tensor
                  
*******************************************************************************
""")

m = Manager()
m.warn_about_basis_change = False 

sb_reference = BasisReferenceOperator(ham.dim,
                                      name="site basis reference")

#
# Calculation of various relaxation tensors
#

ham.protect_basis()
with eigenbasis_of(ham):
    
    RRT = qm.RedfieldRelaxationTensor(ham, sbi, name="Tensor 1")
    
    print("\nRates from the full relaxation tensor")
    for i in range(1, ham.dim):
        for j in range(1, ham.dim):
            print(i, "<-", j, ":", 1.0/numpy.real(RRT.data[i,i,j,j]))
        
    print("\nCalculating relaxation rates")
    
    RRM = qm.RedfieldRateMatrix(ham, sbi)
    print("\nRates from the rate matrix")
    for i in range(1,ham.dim):
        for j in range(1, ham.dim):
            print(i, "<-", j, ":", 1.0/RRM.data[i,j])

    print("\nComparison of the results: ratio of rates")
    for i in range(1, ham.dim):
        for j in range(1, ham.dim):
            print(i, "<-", j, ":", RRM.data[i,j]/numpy.real(RRT.data[i,i,j,j]))

    TDRRM = qm.TDRedfieldRateMatrix(ham, sbi)
    print("\nRates from the rate matrix")
    for i in range(1,ham.dim):
        for j in range(1, ham.dim):
            print(i, "<-", j, ":", 1.0/TDRRM.data[time.length-1,i,j])

    RRT2 = qm.RedfieldRelaxationTensor(ham, sbi, name="Tensor 2")
    print("\nRates from the full relaxation tensor")
    for i in range(1, ham.dim):
        for j in range(1, ham.dim):
            print(i, "<-", j, ":", 1.0/numpy.real(RRT2.data[i,i,j,j]))

    #RRT2.protect_basis()
    #RRT.protect_basis()

ham.unprotect_basis()
with eigenbasis_of(ham):
    
    #
    # Evolution of reduced density matrix
    #

    prop = ReducedDensityMatrixPropagator(time, ham, RRT2)

    rho_i = ReducedDensityMatrix(dim=4, name="Initial DM")
    rho_i.data[3,3] = 1.0
   
    # FIXME: unprotecting does not work correctly
    #RRT.unprotect_basis()
    
    with eigenbasis_of(sb_reference):
        print(" Rate site basis: ", 1.0/RRT.data[1,1,2,2])

    RRT2.secularize()
    print(" Rate exciton basis: ", 1.0/RRT.data[1,1,2,2])
    rho_t = prop.propagate(rho_i, name="Redfield evolution")

    rho_t.plot(coherences=False)
    
    rho_i1 = ReducedDensityMatrix(dim=4, name="Initial DM")
    rho_i1.data[3,3] = 1.0   
    
    rho_t1 = prop1.propagate(rho_i1, name="Redfield evolution from aggregate")
    rho_t1.plot(coherences=False)
    
    rho_i2 = ReducedDensityMatrix(dim=4, name="Initial DM")
    rho_i2.data[3,3] = 1.0   
    
    rho_t2 = prop_nonsec.propagate(rho_i1,
             name="Redfield evolution from aggregate, non-secular")
    rho_t2.plot(coherences=False)    
 
#rho_t.plot(coherences=False)



#
#  Evolution of populations
#

prop = PopulationPropagator(time, RRM)
pop_t = prop.propagate([0.0, 0.0, 0.0, 1.0])

import matplotlib.pyplot as plt
plt.plot(time.data, pop_t[:,3],'--r')
plt.show()

print(RRM.data[2,3])
with eigenbasis_of(ham):
    print(RRT.data[2,2,3,3])

   

