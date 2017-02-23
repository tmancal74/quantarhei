# -*- coding: utf-8 -*-

from quantarhei import Molecule, Aggregate
from quantarhei import energy_units, eigenbasis_of
    
en = [0.0, 1.0]
m1 = Molecule("Mol1",en)
m2 = Molecule("Mol2",en)
m3 = Molecule("Mol3",en)
m1.set_dipole(0,1,[1.0, 0.0, 0.0])

ag = Aggregate("Homodimer",maxband=1)
ag.add_Molecule(m1)
ag.add_Molecule(m2)
ag.add_Molecule(m3)


ag.set_resonance_coupling(0,1,0.1)

ag.build(mult=2)

H = ag.get_Hamiltonian()

with energy_units("1/cm"):
    print(H.data.shape)
    print(H)
    
D = ag.get_TransitionDipoleMoment()


with eigenbasis_of(H):
    with energy_units("1/cm"):
        print(H)
    print(D.data[0,1,:])
    print(D.data[0,2,:])
    print(D.data[0,3,:])
    

SS = H.diagonalize()
H.undiagonalize()

print(H)

print(SS)