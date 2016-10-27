# -*- coding: utf-8 -*-

import quantarhei as qr


en = [0.0, 1.0]

M = qr.Molecule("My first two-level molecule",en)

H = M.get_Hamiltonian()

print(H)


