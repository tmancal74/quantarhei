# -*- coding: utf-8 -*-

_show_plots_ = False

print("""
***********************************************************      
*
*
*          Evolution superoperator demo
*
*
***********************************************************
""")


import quantarhei as qr


with qr.energy_units("1/cm"):
    mol1 = qr.Molecule([0.0, 12000])
    mol2 = qr.Molecule([0.0, 12000])
    mol3 = qr.Molecule([0.0, 12100])
    mol4 = qr.Molecule([0.0, 12100])

agg = qr.Aggregate([mol1, mol2, mol3, mol4])

Ndim = 5

K12 = qr.qm.ProjectionOperator(1, 2, dim=Ndim)
K23 = qr.qm.ProjectionOperator(2, 3, dim=Ndim)
K34 = qr.qm.ProjectionOperator(3, 4, dim=Ndim)
sbi = qr.qm.SystemBathInteraction(sys_operators=[K12, K23, K34], 
                                  rates=(1.0/200, 1.0/100.0, 1.0/150))

agg.set_SystemBathInteraction(sbi)

agg.build()

ham = agg.get_Hamiltonian()
rt = qr.qm.LindbladForm(ham,  sbi, as_operators=False)


#
# Time interval for evolution superoperator
#
time_so = qr.TimeAxis(0.0, 11, 50.0)


#
# calculation of the evolution superoperator
#
eS = qr.qm.EvolutionSuperOperator(time_so, ham, rt)
eS.set_dense_dt(100)

eS.calculate()


#
# Control calculations: via propagation
#
time_tot = qr.TimeAxis(0.0, 1000, 0.5)

rho0 = qr.ReducedDensityMatrix(dim=Ndim)
rho0.data[4,4] = 0.5
rho0.data[3,3] = 0.5
rho0.data[3,4] = 0.5
rho0.data[4,3] = 0.5

prop_tot = qr.ReducedDensityMatrixPropagator(time_tot, ham, rt)

rhot = prop_tot.propagate(rho0)

if _show_plots_:
    rhot.plot(show=False)

#
# Reconstruct density matrix by applying evolution superoperator
#
rhoc = qr.qm.DensityMatrixEvolution(time_so)
rhoc.set_initial_condition(rho0)

for j in range(1, time_so.length):
    
    tj = time_so.data[j]
    sU = eS.at(tj)
    rdm0 = qr.qm.ReducedDensityMatrix(dim=Ndim, data=rho0.data)
    
    rhoc.data[j, :, :] = sU.apply(rdm0).data

if _show_plots_:
    rhoc.plot()
    

print(eS.at(100.0).data.shape)