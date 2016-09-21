# -*- coding: utf-8 -*-

import quantarhei as qr
import matplotlib.pyplot as plt
import numpy

ta = qr.TimeAxis(0.0, 1000, 1.0)

"""

    Absorption of a monomeric two-level molecule


"""
cfce_params1 = dict(ftype="OverdampedBrownian",
                   reorg=20.0,
                   cortime=100.0,
                   T=100,matsubara=20)

en = 10000.0

with qr.energy_units("1/cm"):
    m = qr.Molecule("Molecule",[0.0,en])
    with qr.energy_units("1/cm"):
        cfce1 = qr.CorrelationFunction(ta,cfce_params1)
    
m.set_egcf((0,1),cfce1)   
m.set_dipole(0,1,[0.0, 1.0, 0.0])

a1 = qr.AbsSpect(ta,m) 

with qr.energy_units("1/cm"):    
    a1.calculate(rwa=en)


with qr.frequency_units("1/cm"):
    a1.plot()
    
save_load = False

if save_load:

    with qr.frequency_units("1/cm"):
        a1.save("some_file",ext="npy")
    
    with qr.frequency_units("1/cm"):
    
        f = qr.DFunction()
    
        f.load("some_file",ext="npy",axis="frequency")
    
        f.plot()
    
"""

    Absorption of a simple trimeric aggregate of two-level molecules


"""