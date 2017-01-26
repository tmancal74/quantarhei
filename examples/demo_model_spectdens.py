# -*- coding: utf-8 -*-

import numpy
from quantarhei.models.spectdens import SpectralDensityDB 
from quantarhei.models.spectdens import CorrelationFunctionDB
from quantarhei import SpectralDensity
from quantarhei import TimeAxis, energy_units

axis = TimeAxis(0.0, 10000, 1.0)

db = CorrelationFunctionDB()
#sdw = db.get_SpectralDensity(axis, "Wendling_JPCB_104_2000_5825")
sdw = db.get_SpectralDensity(axis, "Renger_JCP_2002")

params = dict(ftype="OverdampedBrownian", reorg=300.0, cortime=100.0)
with energy_units("1/cm"):
    sdob = SpectralDensity(axis, params)

ax = sdw.axis
sdob.axis = ax

sd_tot = sdob + sdw

with energy_units("1/cm"):
    sd_tot.plot(axis=[0, 1300, 0.0, numpy.max(sd_tot.data)])
    
    #fc = sd_tot.get_CorrelationFunction(temperature=300)
    #fc.plot()
    
    #fc.save("corfce.dat")