


# -*- coding: utf-8 -*-
"""
Calculation of 2D spectra 



"""
import matplotlib.pyplot as plt
import numpy
import time

from quantarhei import Molecule
from quantarhei import Aggregate
from quantarhei import energy_units
from quantarhei import TimeAxis
from quantarhei import CorrelationFunction

from quantarhei.spectroscopy.twod2 import TwoDSpectrumCalculator
from quantarhei.spectroscopy.twod2 import TwoDSpectrumContainer
from quantarhei.spectroscopy.twod2 import TwoDSpectrum

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
# dimer 1
mol2.position = [0.0, 6.0, 0.0]
## dimer 2
#mol2.position = [0.0, 0.0, 6.0]
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
    
params2 = dict(ftype="UnderdampedBrownian", T=300, reorg=10.0, freq=150.0,
               gamma=1.0/10000.0)
with energy_units('1/cm'):
    cf2 = CorrelationFunction(tsbi, params2)
    
cf.add_to_data(cf2)

#
# Set system-bath interaction
#
mol1.set_transition_environment((0,1),cf)
mol2.set_transition_environment((0,1),cf)

#
# Creating aggregate
#      
agg = Aggregate(name="Dimer", molecules=[mol1, mol2])
agg.set_coupling_by_dipole_dipole()
with energy_units("1/cm"):
    print(agg.get_resonance_coupling(0,1))
agg.build(mult=2)

with energy_units("1/cm"):
    rwa_cm = agg.get_RWA_suggestion()
rwa = agg.get_RWA_suggestion()

#
# Prepare for calculation of 2D spectra
#

# TimeAxis for t2 waiting time

t2s = TimeAxis(0.0, 5, 20.0)

#
# Set up calculator
#
tcalc = TwoDSpectrumCalculator(t1axis=ta, t2axis=t2s, t3axis=ta,
                               system=agg)

tcalc.bootstrap(rwa, verbose=True, lab=lab)


#
# Calculate 2D spectra, display and save them
#    

w1_min = rwa_cm - 700.0
w1_max = rwa_cm + 700.0
w3_min = rwa_cm - 700.0
w3_max = rwa_cm + 700.0

window_2D = [w1_min, w1_max, w3_min, w3_max]
t_start = time.time()

#
# Plotting with given units on axes
#
with energy_units("1/cm"):
    
    k = 0
    twods = TwoDSpectrumContainer(t2s)
    
    for tt2 in t2s.data:
        
        # calculate spectra iteratively
        spect = tcalc.calculate_next()
        # save memory by trimming the spectrum
        spect.trim_to(window=window_2D)
        # save it to a container
        twods.set_spectrum(spect)
        # plot it
        #spect.plot(window=window_2D) #,vmax=1.0, cbmax=cbmax)
        # save figure
        #figname = "fig"+str(round(tt2))+".png"
        #print("saving file: ", figname, " with 2D spectrum at ", tt2, "fs")
        #spect.savefig(figname)
        #spect.show()
        
        k += 1
        
t_end = time.time()
t_fend = time.time()
print("Finished at ", t_fend-t_fstart, " secs")
print("Calculation took ", t_end-t_start, " sec")

#pp1 = numpy.zeros(t2s.length)
#pp2 = numpy.zeros(t2s.length)
#k = 0
with energy_units("1/cm"):
        
    pp1 = twods.get_point_evolution(12250,12250,t2s)
    pp2 = twods.get_point_evolution(11900,12250,t2s)    
    pp3 = twods.get_point_evolution(11900,11900,t2s)
    pp4 = twods.get_point_evolution(12000,12500,t2s)
    
    
plt.plot(t2s.data,pp1)
plt.plot(t2s.data,pp2)
plt.plot(t2s.data,pp3)
plt.plot(t2s.data,pp4)
plt.savefig("points.png")

#sp = twods.get_spectrum(t2s.data[-1])
#with energy_units("1/cm"):
#    sp.plot(window=window_2D)
#    
#sp.save("spectrum.hdf5")
#
#rsp = TwoDSpectrum()
#rsp.load("spectrum.hdf5")
#
#with energy_units("1/cm"):
#    rsp.plot(window=window_2D) 

with energy_units("1/cm"):
    twods.trimall_to(window=window_2D)
twods.save("allspectra.hdf5")

#newtw = TwoDSpectrumContainer()
#newtw.load("allspectra.hdf5")
#
#sp = newtw.get_spectrum(t2s.data[-2])
#with energy_units("1/cm"):
#    sp.plot(window=window_2D)


        
    
