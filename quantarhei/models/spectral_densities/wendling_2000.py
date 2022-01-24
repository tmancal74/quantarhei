# -*- coding: utf-8 -*-

from quantarhei.models.spectdens import SpectralDensityDatabaseEntry 
from quantarhei import SpectralDensity

from quantarhei import energy_units


class wendling_2000a(SpectralDensityDatabaseEntry):
    """Spectral density of bacteriochlorophyll in FMO complex 
    
    Spectral density derived in 
    
    Markus Wendling, Tonu Pullerits, Milosz A. Przyjalgowski,
    Simone I. E. Vulto, Thijs J. Aartsma, Rienk van Grondelle,
    and Herbert van Amerongen, "Electron-Vibrational Coupling 
    in the Fenna-Matthews-Olson Complex of Prosthecochloris aestuarii
    Determined by Temperature-Dependent Absorption and Fluorescence
    Line-Narrowing Measurements"
    J. Phys. Chem. B, 2000, 104, 5825-5831
    
    
    """
    def __init__(self):
        
        self.identificator = "Wendling_JPCB_104_2000_5825"
        self.alt_ident = ["Wendling", "Wendling2000"]
        
        
    def get_SpectralDensity(self, axis):

        #
        # Data from the paper
        #
        omegas = [36,   70,   117, 173, 185, 195, 237, 260, 284, 327,
                  365, 381, 479, 541, 565, 580, 635, 714, 723, 730,
                  747, 759, 768, 777, 819, 859, 896, 1158, 1176, 1216]
        fcfs   = [0.01, 0.01, 0.0055, 0.008, 0.008, 0.011, 0.005, 
                  0.0025, 0.005, 0.0015,
                  0.002, 0.002, 0.001, 0.001, 0.002, 0.001, 0.003,
                  0.002, 0.003, 0.001, 0.002, 0.002, 0.004, 0.0015,
                  0.002, 0.0025, 0.002, 0.004, 0.003, 0.002]
        
        #
        # Guessed dephasing times of the oscillators
        #
        gammas = [1.0/3000.0]*len(omegas)
        
        data = []
        for i in range(len(omegas)):
            data.append((omegas[i],fcfs[i],gammas[i]))
            
        #
        # All contributions consist of underdamped harmonic oscillators
        #
        params = dict(ftype="UnderdampedBrownian")
        
        #
        # Create and sum-up spectral densities
        #
        k = 0
        for d in data:
            
            params["freq"] = d[0]
            params["reorg"] = d[0]*d[1]
            params["gamma"] = d[2] #1.0/3000.0
            
            with energy_units("1/cm"):
                sd = SpectralDensity(axis, params)
                
            if k == 0:
                ret = sd
                ax = sd.axis
            else:
                sd.axis = ax
                ret = ret + sd
            k +=1

        return ret


