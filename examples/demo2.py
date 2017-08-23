# -*- coding: utf-8 -*-

from quantarhei import Molecule, Aggregate
from quantarhei import energy_units
    
en = [0.0, 1.0]
m1 = Molecule("Mol1",en)
m2 = Molecule("Mol2",en)

ag = Aggregate("Homodimer",maxband=1)
ag.add_Molecule(m1)
ag.add_Molecule(m2)


ag.set_resonance_coupling(0,1,0.1)

ag.build()

H = ag.get_Hamiltonian()

with energy_units("1/cm"):
    print(H)