# -*- coding: utf-8 -*-

#<remove>
_show_plots_ = False
#</remove>

from quantarhei import Molecule, Aggregate
#from quantarhei import energy_units
    
en = [0.0, 1.0]
m1 = Molecule(name="Mol1",elenergies=en)
m2 = Molecule(name="Mol2",elenergies=en)

ag = Aggregate(name="Homodimer")
ag.add_Molecule(m1)
ag.add_Molecule(m2)


ag.set_resonance_coupling(0,1,0.1)

ag.build(mult=1)

H = ag.get_Hamiltonian()

#with energy_units("1/cm"):
#    print(H)