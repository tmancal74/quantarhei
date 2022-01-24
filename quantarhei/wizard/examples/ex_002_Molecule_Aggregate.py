# -*- coding: utf-8 -*-

#<remove>
_show_plots_ = False
#</remove>

import quantarhei as qr

    
en = [0.0, 1.0]
m1 = qr.Molecule(name="Mol1",elenergies=en)
m2 = qr.Molecule(name="Mol2",elenergies=en)

ag = qr.Aggregate(name="Homodimer")
ag.add_Molecule(m1)
ag.add_Molecule(m2)


ag.set_resonance_coupling(0,1,0.1)

ag.build(mult=1)

H = ag.get_Hamiltonian()

#with qr.energy_units("1/cm"):
#    print(H)

#
# Here we test generation of states with 3 level molecules
#

en = [0.0, 10100.0] #, 20200.0]
with qr.energy_units("1/cm"):
    m1 = qr.Molecule(name="Mol1",elenergies=en)
    m2 = qr.Molecule(name="Mol2",elenergies=en)
    m3 = qr.Molecule(name="Mol3",elenergies=en)

ag = qr.Aggregate(name="Trimer-3-lev")
ag.add_Molecule(m1)
ag.add_Molecule(m2)
ag.add_Molecule(m3)

ii = 0
for sig in ag.elsignatures(mult=4):
    print(ii, sig)
    ii += 1
    
ag.build(mult=2)

H = ag.get_Hamiltonian()

with qr.energy_units("1/cm"):
    print(H)
    print(H.dim)