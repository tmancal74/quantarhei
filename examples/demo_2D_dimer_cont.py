# -*- coding: utf-8 -*-
from quantarhei.spectroscopy.twod2 import TwoDSpectrumContainer
from quantarhei import energy_units

rwa_cm = 12000
w1_min = rwa_cm - 700.0
w1_max = rwa_cm + 700.0
w3_min = rwa_cm - 700.0
w3_max = rwa_cm + 700.0

window_2D = [w1_min, w1_max, w3_min, w3_max]

newtw = TwoDSpectrumContainer()
newtw.load("allspectra.hdf5")

sp = newtw.get_spectrum(40)
with energy_units("1/cm"):
    sp.plot(axis=window_2D)
