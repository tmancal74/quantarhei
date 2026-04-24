# -*- coding: utf-8 -*-

from quantarhei import Molecule, Aggregate, CorrelationFunction
from quantarhei import TimeAxis
from quantarhei import energy_units, eigenbasis_of
  
#
# Define and build the system
#  
en = [0.0, 1.0]
m1 = Molecule("Mol1",en)
m2 = Molecule("Mol2",en)
m3 = Molecule("Mol3",en)
m1.set_dipole(0,1,[1.0, 0.0, 0.0])

time = TimeAxis(0.0, 1000, 1.0)
bath_params = dict(ftype="OverdampedBrownian", T=300, cortime=100, reorg=30.0)
with energy_units("1/cm"):
    cf = CorrelationFunction(time, bath_params)
m1.set_transition_environment((0,1),cf)
m2.set_transition_environment((0,1),cf)
m3.set_transition_environment((0,1),cf)

ag = Aggregate("Homodimer")
ag.add_Molecule(m1)
ag.add_Molecule(m2)
ag.add_Molecule(m3)

ag.set_resonance_coupling(0,1,0.1)

mult = 2

ag.build(mult=mult,sbi_for_higher_ex=False)



#
# Look at its various components
#
H = ag.get_Hamiltonian()
print("Shape of the operators")
print(H.data.shape)
with energy_units("1/cm"):
    print(H)
    
D = ag.get_TransitionDipoleMoment()


with eigenbasis_of(H):
    with energy_units("1/cm"):
        print(H)
    print("\nTransition dipole moments")
    print(D.data[0,1,:])
    print(D.data[0,2,:])
    print(D.data[0,3,:])
    

SS = H.diagonalize()
H.undiagonalize()

print("\nHamiltonian in internal units")
print(H)

print("\nTransformation matrix")
print(SS)

print("\nAll states")
k = 0
for st in ag.allstates(mult=mult, save_indices=True):
    print(st[0], ag.elsigs[k])
    k += 1

sbi = ag.get_SystemBathInteraction()
cfm = sbi.CC
print("\nReorganization energies")
for ic in range(cfm.nof):
    with energy_units("1/cm"):
        print(ic, cfm.cfuncs[ic+1].get_reorganization_energy())
    
print("\nSBI:")
print(sbi)

print(ag)

print(ag.which_band)
print(ag.vibindices)