# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy
import time

from quantarhei import Molecule
from quantarhei import Aggregate
from quantarhei import energy_units
from quantarhei import TimeAxis
from quantarhei import CorrelationFunction

from quantarhei.spectroscopy.twod2 import TwoDSpectrumCalculator
from aceto.lab_settings import lab_settings

from quantarhei import Manager

print(Manager().version)

t_fstart = time.time()

# Time axis for t1 and t3 times
Nr = 1000
ta = TimeAxis(0.0, Nr, 2.0)

###############################################################################
#
# Define problem
#
###############################################################################

#
# define molecules
#
with energy_units("1/cm"):
    mol1 = Molecule(elenergies=[0.0, 12300.0])
    mol2 = Molecule(elenergies=[0.0, 12000.0])
mol1.position = [0.0, 0.0, 0.0]
## dimer 1
#mol2.position = [0.0, 6.0, 0.0]
# dimer 2
mol2.position = [0.0, 0.0, 6.0]
mol1.set_dipole(0,1,[4.0, 2.0, 0.0])
mol2.set_dipole(0,1,[1.0, 2.0, 0.0])

#
# Setting up laboratory
#
lab = lab_settings(lab_settings.FOUR_WAVE_MIXING)
a_0 = numpy.array([1.0, 0.0, 0.0], dtype=numpy.float64)
lab.set_laser_polarizations(a_0,a_0,a_0,a_0)

#
# System-bath interaction
#      
tsbi = TimeAxis(0.0, 3*Nr, 2.0)
params = dict(ftype="OverdampedBrownian", T=300, reorg=70.0, cortime=100.0)
with energy_units('1/cm'):
    cf = CorrelationFunction(tsbi, params)
    
#
# Set system-bath interaction
#
mol1.set_transition_environment((0,1),cf)
mol2.set_transition_environment((0,1),cf)

#
# Creating aggregate
#      
agg = Aggregate("Dimer", molecules=[mol1, mol2])
agg.set_coupling_by_dipole_dipole()
with energy_units("1/cm"):
    print(agg.get_resonance_coupling(0,1))
agg.build(mult=2)

with energy_units("1/cm"):
    rwa_cm = agg.get_RWA_suggestion()
rwa = agg.get_RWA_suggestion()

#
# Calculate 2D spectra
#

# TimeAxis for t2 waiting time

t2s = TimeAxis(0.0, 2, 100.0)

tcalc = TwoDSpectrumCalculator(t1axis=ta, t2axis=t2s, t3axis=ta,
                               system=agg)

t_start = time.time()
print("Bootstrapping")
tcalc.bootstrap(rwa, verbose=True, lab=lab)
print("...done")
print("Calculating")
twods = tcalc.calculate()
print("...done")
t_end = time.time()

#
# Show 2D spectra (save them)
#    

w1_min = rwa_cm - 700.0
w1_max = rwa_cm + 700.0
w3_min = rwa_cm - 700.0
w3_max = rwa_cm + 700.0

window_2D = [w1_min, w1_max, w3_min, w3_max]
print(window_2D)

with energy_units("1/cm"):
    
    k = 0
    for tt2 in t2s.data:
            
        #
        # Plotting with given units on axes
        #
        twods[k].plot(axis=[w1_min, w1_max, w3_min, w3_max]) #,vmax=1.0, cbmax=cbmax)

        figname = "fig"+str(round(tt2))+".png"
        print("saving file: ", figname, " with 2D spectrum at ", tt2, "fs")
        plt.savefig(figname)
        plt.show()
        
        k += 1

t_fend = time.time()
print("Finished at ", t_fend-t_fstart, " secs")
print("Calculation took ", t_end-t_start, " sec")