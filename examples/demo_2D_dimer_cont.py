# -*- coding: utf-8 -*-
from quantarhei.spectroscopy.twod2 import TwoDSpectrumContainer
from quantarhei import energy_units

rwa_cm = 12000
w1_min = rwa_cm - 700.0
w1_max = rwa_cm + 700.0
w3_min = rwa_cm - 700.0
w3_max = rwa_cm + 700.0

window_2D = [w1_min, w1_max, w3_min, w3_max]
off = 0
window_trim = [w1_min-off, w1_max+off, w3_min-off, w3_max+off]

newtw = TwoDSpectrumContainer()
newtw.load("allspectra.hdf5")

sp = newtw.get_spectrum(40)

sp.save("sp_orig.hdf5")

with energy_units("1/cm"):
    sp.plot(window=window_2D)

    sp.trim_to(window=window_trim)
    sp.save("sp_trimmed.hdf5")

    sp.plot(window=window_2D)

    newtw.trimall_to(window=window_trim)
    newtw.save("allspectra_trimmed.hdf5")


testtw = TwoDSpectrumContainer()
testtw.load("allspectra_trimmed.hdf5")

#with energy_units("1/cm"):
print("Plotting all trimmed spectra:\n")
for s in testtw.get_spectra():
    s.plot()