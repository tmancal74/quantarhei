# -*- coding: utf-8 -*-
_show_plots_ = False
import numpy

import quantarhei as qr

qr.Manager().gen_conf.legacy_relaxation = True

print("Preparing a model system:")

with qr.energy_units("1/cm"):
    mol1 = qr.Molecule([0.0, 12010])
    mol2 = qr.Molecule([0.0, 12000])
    mol3 = qr.Molecule([0.0, 12100])
    mol4 = qr.Molecule([0.0, 12110])

agg = qr.Aggregate([mol1, mol2, mol3, mol4])
agg.set_resonance_coupling(2,3,qr.convert(100.0,"1/cm","int"))
agg.set_resonance_coupling(1,3,qr.convert(100.0,"1/cm","int"))
agg.set_resonance_coupling(1,2,qr.convert(0.0,"1/cm","int"))

#qr.save_parcel(agg,"agg.qrp")
#agg2 = qr.load_parcel("agg.qrp")
agg2 = agg.deepcopy()
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

#
# Reconstruct density matrix by applying evolution superoperator
#
rhoc1 = qr.qm.DensityMatrixEvolution(time_so)
rhoc1.set_initial_condition(rho0)

for j in range(1, time_so.length):
    
    tj = time_so.data[j]
    sU = eS.at(tj)
    rdm0 = qr.qm.ReducedDensityMatrix(dim=Ndim, data=rho0.data)
    
    rhoc1.data[j, :, :] = sU.apply(rdm0).data

rhoc_0_data = numpy.zeros((Ndim, Ndim), dtype=qr.COMPLEX)
rhoc_1_data = numpy.zeros((Ndim, Ndim), dtype=qr.COMPLEX)
rhoc_2_data = numpy.zeros((Ndim, Ndim), dtype=qr.COMPLEX)
rhot_2_data = numpy.zeros((Ndim, Ndim), dtype=qr.COMPLEX)    

print(ham)
if _show_plots_:
    print("Comparing in a plots 1:")
    with qr.eigenbasis_of(ham):
    #if True:
        rhot.plot(populations=False, show=False)
        rhoc1.plot(populations=False)
        rhoc_0_data[:,:] = rhoc1.data[0,:,:]
    print("...done")

with qr.eigenbasis_of(ham):
    rhoc_1_data[:,:] = rhoc1.data[0,:,:]

print("Difference: ", numpy.amax(numpy.abs(rhoc_1_data - rhoc_0_data)))    

D = None

eS = qr.qm.EvolutionSuperOperator(time_so, ham, rt, pdeph=D)
eS.set_dense_dt(100)
with qr.eigenbasis_of(ham):
    eS.calculate()


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

print(ham)
if _show_plots_:
    print("Comparing in plots 2:")
    with qr.eigenbasis_of(ham):
    #if True:
        #rhoc.plot(coherences=True, populations=False, show=False)
        rhot.plot(coherences=True, populations=False, show=False)
        rhoc_2_data[:,:] = rhoc1.data[0,:,:]
        rhot_2_data[:,:] = rhot.data[0,:,:]
        rhoc1.plot(populations=False)
        
        #rhoc.plot(show=False)
        #rhot.plot(show=False)
        #rhoc1.plot()

    print("...done")    

print("Difference: ", numpy.amax(numpy.abs(rhoc_2_data - rhoc_0_data)))   
print("Difference: ", numpy.amax(numpy.abs(rhot_2_data - rhoc_0_data))) 


eS = qr.qm.EvolutionSuperOperator(time_so, ham, rt)
eS.set_dense_dt(100)

eS.calculate()

prop_tot = qr.ReducedDensityMatrixPropagator(time_tot, ham, rt)

rhot = prop_tot.propagate(rho0)

#
# Reconstruct density matrix by applying evolution superoperator
#
rhoc1 = qr.qm.DensityMatrixEvolution(time_so)
rhoc1.set_initial_condition(rho0)

for j in range(1, time_so.length):
    
    tj = time_so.data[j]
    sU = eS.at(tj)
    rdm0 = qr.qm.ReducedDensityMatrix(dim=Ndim, data=rho0.data)
    
    rhoc1.data[j, :, :] = sU.apply(rdm0).data
    
if _show_plots_:
    print("Comparing in plots 2:")
    with qr.eigenbasis_of(H):
    #if True:
        #rhoc.plot(coherences=True, populations=False, show=False)
        rhot.plot(coherences=True, populations=False, show=False)
        rhoc_2_data[:,:] = rhoc1.data[0,:,:]
        rhot_2_data[:,:] = rhot.data[0,:,:]
        rhoc1.plot(populations=False)

print(H)
print(ham)
print("Copy")
diff = H.data - ham.data
ham.data = H.data[:,:]
print(numpy.max(numpy.abs(H.data-ham.data)))
print(numpy.linalg.eig(H.data))
print(numpy.linalg.eig(ham.data))
print(numpy.linalg.eig(ham.data - diff))
print(diff)
print(numpy.max(diff))