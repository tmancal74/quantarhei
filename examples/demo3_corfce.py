# -*- coding: utf-8 -*-

from quantarhei import CorrelationFunction
from quantarhei import TimeAxis
from quantarhei import energy_units


ta = TimeAxis(0.0,100,10.0)

params = {"ftype":    "OverdampedBrownian",
          "reorg":    20.0,
          "cortime":  100.0,
          "T":        100.0,
          "matsubara":20}
          
with energy_units("1/cm"):
    cf = CorrelationFunction(ta,params)

cf.plot(ylabel=r'$C(t)$')

    
        