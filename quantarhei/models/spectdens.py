# -*- coding: utf-8 -*-
import numpy
from ..qm.corfunctions import SpectralDensity, CorrelationFunction
from ..core.managers import energy_units

class SpectralDensityDB:
    
    def __init__(self):
        pass
    
    
    
    def get_SpectralDensity(self, axis, ident=None):
        """
        
        """
        
        if ident is not None:
            
            if ident == "Wendling_JPCB_104_2000_5825":
                
                data = []
                omegas = [36,   70,   117, 173, 185, 195, 237, 260, 284, 327 ]
                fcs    = [0.01, 0.01, 0.0055, 0.008, 0.008, 0.011, 0.005, 0.0025, 0.005, 0.0015]
                for i in range(len(omegas)):
                    data.append((omegas[i],fcs[i]))
                    
                params = dict(ftype="UnderdampedBrownian",
                              gamma=1.0/3000.0)
                
                k = 0
                for d in data:
                    params["freq"] = d[0]
                    params["reorg"] = d[0]*d[1]
                    with energy_units("1/cm"):
                        sd = SpectralDensity(axis, params)
                    if k == 0:
                        ret = sd
                        ax = sd.axis
                    else:
                        sd.axis = ax
                        ret += sd
                    k += 1
                    
            else:
                
                Exception("Unknown spectral density")
            
        else:
            Exception("Spectral density identificator not specified")

        return ret

    
