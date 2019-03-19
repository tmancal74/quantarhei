# -*- coding: utf-8 -*-

import numpy
from quantarhei.models.spectdens import SpectralDensityDB 
#from quantarhei.models.spectdens import CorrelationFunctionDB
from quantarhei import SpectralDensity
from quantarhei import TimeAxis, energy_units

_show_plots_ = True

axis = TimeAxis(0.0, 10000, 1.0)

db = SpectralDensityDB(verbose=True)

sdw = db.get_SpectralDensity(axis, "Wendling_JPCB_104_2000_5825")
#sdw = db.get_SpectralDensity(axis, "Renger_JCP_2002")

params = dict(ftype="OverdampedBrownian", reorg=300.0, cortime=100.0)
with energy_units("1/cm"):
    sdob = SpectralDensity(axis, params)

ax = sdw.axis
sdob.axis = ax

sd_tot = sdob + sdw

with energy_units("1/cm"):
    if _show_plots_:
        sd_tot.plot(axis=[0, 1300, 0.0, numpy.max(sd_tot.data)])

    # here we get correlation function at a given temperature in K
    cf = sd_tot.get_CorrelationFunction(300.0)
    
print(db.get_status_string())

if _show_plots_:
    cf.plot()