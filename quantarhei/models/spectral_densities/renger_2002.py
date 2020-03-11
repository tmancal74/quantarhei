# -*- coding: utf-8 -*-

from quantarhei.models.spectdens import SpectralDensityDatabaseEntry 
from quantarhei import SpectralDensity
from quantarhei import energy_units


class renger_2002a(SpectralDensityDatabaseEntry):
    """
    Spectral density of bacteriochlorophyll

    Spectral density derived in 
    
    T. Renger and R. A. Marcus, "On the relation of protein dynamics and 
    exciton relaxation in pigment–protein complexes: An estimation of 
    the spectral density and a theory for the calculation of optical 
    spectra" Journal of Chemical Physics, 2002, 116, 9997-10017
    
    """
    def __init__(self):
        
        self.identificator = "Renger_JCP_2002"
        self.alt_ident = ["Renger", "Renger2002"]
        
        
    def get_SpectralDensity(self, axis):
    
        # Calling the spec dens calculated from b777 in rengers paper
        # alternative_form gives the polynomial version of the same spec dens
        # Jang, Newton, Silbey, J Chem Phys. 2007.
        params = {"ftype": "B777",
                  "reorg": 102.0,
                  "alternative_form": True,
                  "T":300,
                  "cortime":100}

        # can vary the paramaters, values from the paper are set by default
        params.update({"freq1": 0.56,
                       "freq2": 1.9,
                        's1':   0.8,
                        's2':   0.5,
                        'om1':  170,
                        'om2':  34,
                        'om3':  69,})

        with energy_units("1/cm"):
            sd = SpectralDensity(axis, params)

        return sd
