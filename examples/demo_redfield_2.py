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

time = TimeAxis(0.0, 3000, 1.0)

with energy_units("1/cm"):

    m1 = Molecule("Mol 1", [0.0, 10100.0])
    m2 = Molecule("Mol 2", [0.0, 10050.0])
    m3 = Molecule("Mol 3", [0.0, 10000.0])
    
    m1.position = [0.0, 0.0, 0.0]
    m2.position = [10.0, 0.0, 0.0]
    m3.position = [5.0, 5.0, 0.0]
    m1.set_dipole(0,1,[5.8, 0.0, 0.0])
    m2.set_dipole(0,1,[5.8, 0.0, 0.0])
    m3.set_dipole(0,1,[numpy.sqrt(12.0), 0.0, 0.0])

    agg = Aggregate("Dimer")
    agg.add_Molecule(m1)
    agg.add_Molecule(m2)
    agg.add_Molecule(m3)
    
    #agg.set_resonance_coupling(0,1,30.0)
    #agg.set_resonance_coupling(1,2,30.0)
    
    agg.set_coupling_by_dipole_dipole(epsr=1.92921)

    params = dict(ftype="OverdampedBrownian", reorg=20, cortime=300, T=100)
    cf = CorrelationFunction(time, params)

m1.set_transition_environment((0,1), cf)
m2.set_transition_environment((0,1), cf)
m3.set_transition_environment((0,1), cf)


agg.build()

sbi = agg.get_SystemBathInteraction()
ham = agg.get_Hamiltonian()

ham.protect_basis()
with eigenbasis_of(ham):
    
    TRRT = qm.TDRedfieldRelaxationTensor(ham, sbi, name="Tensor 1",
                                         cutoff_time=1300.0)
 
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
            
ham.unprotect_basis()
with eigenbasis_of(ham):
    
    #
    # Evolution of reduced density matrix
    #
    prop = ReducedDensityMatrixPropagator(time, ham, RRT)
    prop_t = ReducedDensityMatrixPropagator(time, ham, TRRT)

    #TRRT.secularize()
    #RRT.secularize()
    
    rho_i = ReducedDensityMatrix(dim=4, name="Initial DM")
    rho_i.data[3,3] = 1.0
    
    rho_t = prop.propagate(rho_i, name="Redfield evolution")
    rho_tt = prop_t.propagate(rho_i, name="Redfield evolution")

    rho_t.plot(coherences=False)
    rho_tt.plot(coherences=False)
