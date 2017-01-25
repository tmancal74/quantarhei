# -*- coding: utf-8 -*-

import numpy
from quantarhei.models.spectdens import SpectralDensityDB
from quantarhei import TimeAxis, energy_units

db = SpectralDensityDB()
axis = TimeAxis(0.0, 10000, 1.0)
sd = db.get_SpectralDensity(axis, "Wendling_JPCB_104_2000_5825")

with energy_units("1/cm"):
    sd.plot(axis=[0, 1000, 0.0, numpy.max(sd.data)])
    