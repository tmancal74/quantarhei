# -*- coding: utf-8 -*-

#<remove>
_show_plots_ = False
#</remove>

import quantarhei as qr

en = [0.0, 1.0]

M = qr.Molecule(elenergies=en)

H = M.get_Hamiltonian()

print(H)

print("version = ", qr.Manager().version)
