


# -*- coding: utf-8 -*-
"""
Calculation of 2D spectra 



"""
import matplotlib.pyplot as plt
import numpy
import time

#from quantarhei import Molecule
#from quantarhei import Aggregate
#from quantarhei import energy_units
#from quantarhei import TimeAxis
#from quantarhei import CorrelationFunction

import quantarhei as qr

#from quantarhei.spectroscopy.twod2 import TwoDSpectrumCalculator
#from quantarhei.spectroscopy.twod2 import TwoDSpectrumContainer
#from quantarhei.spectroscopy.twod2 import TwoDSpectrum

from aceto.lab_settings import lab_settings

from quantarhei import Manager

print(Manager().version)

t_fstart = time.time()

# Time axis for t1 and t3 times
Nr = 1000
ta = qr.TimeAxis(0.0, Nr, 1.0)

###############################################################################
#
# Define problem
#
###############################################################################

#
# define molecules
#
Nmol = 5
Emol = 12500.0
mols = []
with qr.energy_units("1/cm"):
    for ii in range(Nmol):
        mol = qr.Molecule(elenergies=[0.0, Emol+20.0*numpy.random.randn()]) #
        mols.append(mol)
        
#    mol1 = Molecule(elenergies=[0.0, 12550.0])
#    mol2 = Molecule(elenergies=[0.0, 12150.0])
#    mol3 = Molecule(elenergies=[0.0, 12350.0])
#mol1.position = [0.0, 0.0, 0.0]
## dimer 1
#mol2.position = [0.0, 4.0, 0.0]
### dimer 2
##mol2.position = [0.0, 0.0, 6.0]
#mol1.set_dipole(0,1,[2.0, 0.0, 0.0])
#mol2.set_dipole(0,1,[2.0*numpy.sqrt(1.0/2.0), 2.0*numpy.sqrt(1.0/2.0), 0.0])
#mol3.set_dipole(0,1,[0.0, 2.0, 0.0])

for ii in range(Nmol):
    phi = (2.0*numpy.pi/Nmol)*ii
    dip = [numpy.cos(phi), numpy.sin(phi),0.0]
    mols[ii].set_dipole(0,1,dip)

#
# Setting up laboratory
#
lab = lab_settings(lab_settings.FOUR_WAVE_MIXING)
a_0 = numpy.array([1.0, 0.0, 0.0], dtype=numpy.float64)
lab.set_laser_polarizations(a_0,a_0,a_0,a_0)

#
# System-bath interaction
#
   
tsbi = qr.TimeAxis(0.0, 3*Nr, 1.0)
params = dict(ftype="OverdampedBrownian", T=300, reorg=100.0, cortime=50.0)
with qr.energy_units('1/cm'):
    cf = qr.CorrelationFunction(tsbi, params)
    
params2 = dict(ftype="UnderdampedBrownian", T=300, reorg=10.0, freq=150.0,
               gamma=1.0/10000.0)
with qr.energy_units('1/cm'):
    cf2 = qr.CorrelationFunction(tsbi, params2)
    
#cf.add_to_data(cf2)

##
## Set system-bath interaction
##
#mol1.set_transition_environment((0,1),cf)
#mol2.set_transition_environment((0,1),cf)
#mol3.set_transition_environment((0,1),cf)

for ii in range(Nmol):
    mols[ii].set_transition_environment((0,1), cf)

#
# Creating aggregate
#      
#agg = Aggregate(name="Dimer", molecules=[mol1, mol2, mol3])
agg = qr.Aggregate(molecules=mols)
#agg.set_coupling_by_dipole_dipole()
with qr.energy_units("1/cm"):
    #agg.set_resonance_coupling(0,1,-100.0)
    #agg.set_resonance_coupling(1,2,-100.0)
    if Nmol > 1:
        for ii in range(Nmol-1):
            agg.set_resonance_coupling(ii,ii+1,-100.0)
        agg.set_resonance_coupling(0,Nmol-1,-100.0)

        print(agg.get_resonance_coupling(0,1))
agg.build(mult=2)

print(agg.get_Hamiltonian())

with qr.energy_units("1/cm"):
    rwa_cm = agg.get_RWA_suggestion()
rwa = agg.get_RWA_suggestion()

#agg.diagonalize()

#
# Prepare for calculation of 2D spectra
#

# TimeAxis for t2 waiting time

t2s = qr.TimeAxis(0.0, 10, 20.0)

#
# Set up calculator
#
tcalc = qr.TwoDResponseCalculator(t1axis=ta, t2axis=t2s, t3axis=ta,
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
with qr.energy_units("1/cm"):
    
    k = 0
    twods = qr.TwoDResponseContainer(t2s)
    
    for tt2 in t2s.data:
        
        # calculate spectra iteratively
        spect = tcalc.calculate_next()
        # save memory by trimming the spectrum
        spect.trim_to(window=window_2D)
        # save it to a container
        twods.set_spectrum(spect)
        # plot it
        spect2 = spect.get_TwoDSpectrum()
        spect2.plot(window=window_2D, show=True)
        #,vmax=1.0, cbmax=cbmax)
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
#with energy_units("1/cm"):
#        
#    pp1 = twods.get_point_evolution(12250,12250,t2s)
#    pp2 = twods.get_point_evolution(11900,12250,t2s)    
#    pp3 = twods.get_point_evolution(11900,11900,t2s)
#    pp4 = twods.get_point_evolution(12000,12500,t2s)
#    
#    
#plt.plot(t2s.data,pp1)
#plt.plot(t2s.data,pp2)
#plt.plot(t2s.data,pp3)
#plt.plot(t2s.data,pp4)
#plt.savefig("points.png")
#
##sp = twods.get_spectrum(t2s.data[-1])
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

with qr.energy_units("1/cm"):
    twods.trimall_to(window=window_2D)
twods.save("allspectra.qrp")

#newtw = TwoDSpectrumContainer()
newtw = qr.load_parcel("allspectra.qrp")
#twods = newtw.get_TwoDSpectrumContainer()
#
#sp = twods.get_spectrum(0.0)
#with energy_units("1/cm"):
#    sp.plot(window=window_2D)


        
    
