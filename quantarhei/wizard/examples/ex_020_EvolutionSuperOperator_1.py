# -*- coding: utf-8 -*-
import numpy

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

qr.Manager().gen_conf.legacy_relaxation = True

print("Preparing a model system:")

with qr.energy_units("1/cm"):
    mol1 = qr.Molecule([0.0, 12000])
    mol2 = qr.Molecule([0.0, 12000])
    mol3 = qr.Molecule([0.0, 12100])
    mol4 = qr.Molecule([0.0, 12100])

agg = qr.Aggregate([mol1, mol2, mol3, mol4])
agg.set_resonance_coupling(2,3,qr.convert(100.0,"1/cm","int"))
agg.set_resonance_coupling(1,3,qr.convert(100.0,"1/cm","int"))

qr.save_parcel(agg,"agg.qrp")
agg2 = qr.load_parcel("agg.qrp")
agg2.build()

H = agg2.get_Hamiltonian()


print("...done")

print("Setting up Lindblad form relaxation:")
Ndim = 5
with qr.eigenbasis_of(H):
    K12 = qr.qm.ProjectionOperator(1, 2, dim=Ndim)
    K23 = qr.qm.ProjectionOperator(2, 3, dim=Ndim)
    K34 = qr.qm.ProjectionOperator(3, 4, dim=Ndim)
    sbi = qr.qm.SystemBathInteraction(sys_operators=[K12, K23, K34], 
                                      rates=(1.0/200, 1.0/100.0, 1.0/150))

agg.set_SystemBathInteraction(sbi)

agg.build()

ham = agg.get_Hamiltonian()
rt = qr.qm.LindbladForm(ham,  sbi) # as_operators=False)

print("...done")

#
# Time interval for evolution superoperator
#
time_so = qr.TimeAxis(0.0, 11, 50.0)
time_tot = qr.TimeAxis(0.0, 1000, 0.5)





print("Calculating evolution superoperator:")
#
# calculation of the evolution superoperator
#
eS = qr.qm.EvolutionSuperOperator(time_so, ham, rt)
eS.set_dense_dt(100)

eS.calculate()

print("...done")

print("Reference calculation:")
#
# Reference calculations: via propagation
#

rho0 = qr.ReducedDensityMatrix(dim=Ndim)
rho0.data[4,4] = 0.5
rho0.data[3,3] = 0.5
rho0.data[3,4] = 0.5
rho0.data[4,3] = 0.5

prop_tot = qr.ReducedDensityMatrixPropagator(time_tot, ham, rt)

rhot = prop_tot.propagate(rho0)


if _show_plots_:
    print("Comparing in a plot:")
    with qr.eigenbasis_of(ham):
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
    with qr.eigenbasis_of(ham):
        rhoc.plot()
    print("...done")
    

#print(eS.at(100.0).data.shape)

print("\nThe same with additional pure dephasing:")

rho0 = qr.ReducedDensityMatrix(dim=Ndim)
rho0.data[4,4] = 0.5
rho0.data[3,3] = 0.5
rho0.data[3,4] = 0.5
rho0.data[4,3] = 0.5

print("Setting up pure dephasing:")

td = 30.0

dd = numpy.zeros((5,5), dtype=qr.REAL)
dd[1,2] = 1.0/td
dd[2,1] = 1.0/td
dd[1,3] = 1.0/td
dd[3,1] = 1.0/td
dd[2,3] = 1.0/td
dd[3,2] = 1.0/td
dd[4,3] = 1.0/td
dd[3,4] = 1.0/td
D = None #qr.qm.PureDephasing(drates=dd)
print("...done")

print("Evolution superoperator with pure dephasing:")
eS2 = qr.qm.EvolutionSuperOperator(time_so, ham, rt, pdeph=D, mode="jit")
eS2.set_dense_dt(100)
eS = qr.qm.EvolutionSuperOperator(time_so, ham, rt, pdeph=D)
eS.set_dense_dt(100)

#with qr.eigenbasis_of(ham):
for ii in range(1,time_so.length):
    eS2.calculate_next()
    eS.data[ii,:,:,:,:] = eS2.data

#eS.calculate()


rhoc = qr.qm.DensityMatrixEvolution(time_so)
rhoc.set_initial_condition(rho0)

for j in range(1, time_so.length):
    
    tj = time_so.data[j]
    sU = eS.at(tj)
    rdm0 = qr.qm.ReducedDensityMatrix(dim=Ndim, data=rho0.data)
    
    rhoc.data[j, :, :] = sU.apply(rdm0).data

print("...done")

print("Reference propagation with pure dephasing:")
prop_tot = qr.ReducedDensityMatrixPropagator(time_tot, ham, rt, PDeph=D)
#prop_tot.setDtRefinement(10)

with qr.eigenbasis_of(ham):
    rhot = prop_tot.propagate(rho0)
print("...done")

if _show_plots_:
    print("Comparing in plots:")
    with qr.eigenbasis_of(ham):
    #if True:
        #rhoc.plot(coherences=True, populations=False, show=False)
        #rhot.plot(coherences=True, populations=False)
        rhoc.plot(show=False)
        rhot.plot()

    print("...done")