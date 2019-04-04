# -*- coding: utf-8 -*-
import time

import quantarhei as qr
from quantarhei import LabSetup
from quantarhei.utils.vectors import X #, Y, Z

###############################################################################
#
#
#  PARAMETERS
#
#
###############################################################################

omega = 500.0           # mode frequency in 1/cm
hr_factor = 0.02        # Huag-Rhys factor
Ng = 4                  # number of vibrational states in the ground state
Ne = 4                  # number of vibrational states in the excited state

E1 = 12500.0            # electronic transition energy
width = 100             # spectral width of the electronic transition in 1/cm

Nt2 = 100               # number of steps in t2
dt2 = 5.0               # time step in t2

plot_window = [11500, 13500, 11500, 13500]

fft_of = "total"

###############################################################################
#
#
#  MODEL DEFINITION
#
#
###############################################################################

#
# Create a two-level molecule with one intra-molecular harmonic mode 
#
with qr.energy_units("1/cm"):
    mol = qr.Molecule(elenergies=[0.0, E1])
    mol.set_dipole(0, 1, [1.0, 0.0, 1.0])
    
    mod = qr.Mode(frequency=omega)
    
    mol.add_Mode(mod)
    
    mod.set_nmax(0, Ng)
    mod.set_nmax(1, Ne)
    
    mod.set_HR(1, hr_factor)
    
    mol.set_transition_width((0,1), width)
    
#
# Create an aggregate
#
agg = qr.Aggregate(molecules=[mol])

agg.build()

#
# Time axes and the calculator
#
t1axis = qr.TimeAxis(0.0, 1000, 10.0)
t3axis = qr.TimeAxis(0.0, 1000, 10.0)
t2axis = qr.TimeAxis(0.0, Nt2, dt2)

msc = qr.MockTwoDSpectrumCalculator(t1axis, t2axis, t3axis)
msc.bootstrap(rwa=qr.convert(E1,"1/cm","int"), 
              all_positive=False, shape="Gaussian")

#
# Laboratory setup
#

lab = LabSetup()
lab.set_polarizations(pulse_polarizations=[X,X,X], detection_polarization=X)

#
# Hamiltonian is required
#
H = agg.get_Hamiltonian()

#
# Evolution superoperator
#
eUt = qr.EvolutionSuperOperator(t2axis, H)
eUt.set_dense_dt(10)
eUt.calculate(show_progress=False)

#
# We calulate all 2D spectra here
#
print("Calculating 2D spectra ...")
t1 = time.time()
cont = msc.calculate_all_system(agg, H, eUt, lab, show_progress=True)
t2 = time.time()
print("... done in", t2-t1, "sec")

cont.save("container.qrp")

###############################################################################
#
#
#     OUTPUT
#
#
###############################################################################

#
# Example spectrum to plot and a movie of the time evolution
#

    
twd = cont.get_spectrum(tag=30.0)
with qr.energy_units("1/cm"):
    twd.plot()
    
    print("Creating a movie ...")
    cont.make_movie("movie.mp4")
    print("... done")

#
# Trim the spectra to a smaller region
#
with qr.energy_units("1/cm"): 
    cont.trimall_to(window=plot_window)
    
#
# Window function for subsequenty FFT
#
import quantarhei.functions as func
window = func.Tukey(t2axis, r=0.3, sym=False)

#
# FFT of the spectra
#
print("\nCalculating FFT of the 2D maps")
fcont = cont.fft(window=window, dtype=fft_of, offset=0.0)

#
# Have a look which frequencies we actually have
#
Ndat = len(fcont.axis.data)
print("\nNumber of frequency points:", Ndat)
print("In 1/cm they are:")
with qr.energy_units("1/cm"):
    for k_i in range(Ndat):
        print(k_i, fcont.axis.data[k_i])

#
# Which spectrum we want to see
#
with qr.frequency_units("1/cm"):
    sp1, show_Npoint1 = fcont.get_nearest(500.0)
    sp2, show_Npoint2 = fcont.get_nearest(-500.0)

#
# Plotting the corresponding map
#
units = "1/cm"
with qr.energy_units(units):
    print("\nPlotting spectrum at frequency:", 
          fcont.axis.data[show_Npoint1], units)
    sp1.plot(window=plot_window, Npos_contours=20, 
              stype="total", spart="abs")
    fftfile = "twod_fft_map_1.png"
    sp1.savefig(fftfile)
    print("... saved into: ", fftfile)
    print("\nPlotting spectrum at frequency:", 
          fcont.axis.data[show_Npoint2], units)
    sp2.plot(window=plot_window, Npos_contours=20, 
              stype="total", spart="abs")
    fftfile = "twod_fft_map_2.png"
    sp2.savefig(fftfile)
    print("... saved into: ", fftfile)
    
###############################################################################