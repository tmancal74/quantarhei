# -*- coding: utf-8 -*-
from quantarhei import Molecule, Aggregate
from quantarhei import energy_units
    
en = [0.0, 12500.0]

with energy_units("1/cm"):
    m1 = Molecule("Mol1",en)
    m2 = Molecule("Mol2",en)

ag = Aggregate("Homodimer",maxband=1)
ag.add_Molecule(m1)
ag.add_Molecule(m2)

with energy_units("1/cm"):
    ag.set_resonance_coupling(0,1,100)

ag.build()

H = ag.get_Hamiltonian()

print(H)
with energy_units("THz"):
    print(H)
with energy_units("1/cm"):
    print(H)
    