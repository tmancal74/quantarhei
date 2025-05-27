# -*- coding: utf-8 -*-
import time

import numpy

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
hr_factor = 0.01        # Huag-Rhys factor
Ng = 2                  # number of vibrational states in the ground state
Ne = 2                  # number of vibrational states in the excited state

E1 = 12500.0            # electronic transition energy
width = 100             # spectral width of the electronic transition in 1/cm

Nt2 = 100               # number of steps in t2
dt2 = 10.0               # time step in t2

plot_window = [11500, 13500, 11500, 13500]

dOmega = 50.0

fft_of = qr.signal_REPH  # qr.signal_NONR, qr.signal_TOTL

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

# FIXME: TwoDResponseCalculator
msc = qr.MockTwoDResponseCalculator(t1axis, t2axis, t3axis)
msc.bootstrap(rwa=qr.convert(E1,"1/cm","int"), shape="Gaussian")

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
# FIXME: relaxation in the oscillator
eUt = qr.EvolutionSuperOperator(t2axis, H)
eUt.set_dense_dt(10)
eUt.calculate(show_progress=False)

#
# We calulate all 2D spectra here
#
print("Calculating 2D spectra ...")
t1 = time.time()
    
om1 = qr.convert(omega-dOmega,"1/cm","int")
om2 = qr.convert(omega+dOmega,"1/cm","int")

print(" - positive frequency")
cont1 = msc.calculate_all_system(agg, eUt, lab, 
                                 selection=[["omega2",[om1, om2]],
                                            ["order"]])

print(" - negative frequency")
cont2 = msc.calculate_all_system(agg, eUt, lab, 
                                 selection=[["omega2",[-om2, -om1]],
                                            ["order"]])

t2 = time.time()
print("... done in", t2-t1, "sec")

#cont.save("container.qrp")

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

    
#res = cont1.get_response(tag=0.0)
#specttwd = res.get_TwoDSpectrum(dtype=qr.signal_TOTL)
#with qr.energy_units("1/cm"):
#    res.plot(stype=qr.signal_TOTL)
#    specttwd.plot()
    
#    print("Creating a movie ...")
#    cont.make_movie("movie.mp4")
#    print("... done")

#
# The containers will now keep only the 2D spectrum we want
#
cont1 = cont1.get_TwoDSpectrumContainer(stype=fft_of)
cont2 = cont2.get_TwoDSpectrumContainer(stype=fft_of)
#
# Trim the spectra to a smaller region
#
with qr.energy_units("1/cm"): 
    cont1.trimall_to(window=plot_window)
    cont2.trimall_to(window=plot_window)
    
#
#
# Global fit
#    cont_residue contains the oscillatory residue
#
#params_guess = []
#params_out, cont_residue = cont.global_fit_exponential(params_guess)

    
#
# Window function for subsequenty FFT
#
import quantarhei.functions as func
window = func.Tukey(t2axis, r=0.3, sym=False)

#
# FFT of the spectra
#
print("\nCalculating FFT of the 2D maps")
fcont1 = cont1.fft(window=window, offset=0.0)
fcont2 = cont2.fft(window=window, offset=0.0)


#
# Have a look which frequencies we actually have
#
#Ndat = len(fcont1.axis.data)
#print("\nNumber of frequency points:", Ndat)
#print("In 1/cm they are:")
#with qr.energy_units("1/cm"):
#    for k_i in range(Ndat):
#        print(k_i, fcont1.axis.data[k_i])

#
# Which spectrum we want to see
#
with qr.frequency_units("1/cm"):
    sp21, show_Npoint2 = fcont1.get_nearest(500.0)
    sp11, show_Npoint1 = fcont2.get_nearest(-500.0)


#
# Plotting the corresponding map
#
units = "1/cm"
with qr.energy_units(units):
    
    print("\nPlotting spectrum at frequency:", 
          fcont1.axis.data[show_Npoint1], units)
    sp11.plot(window=plot_window, Npos_contours=20, spart=qr.part_ABS)
    fftfile = "twod_fft_map_1.png"
    sp11.savefig(fftfile)
    print("... saved into: ", fftfile)

    evol1 = fcont2.get_point_evolution(13000, 12500,
                                       fcont1.axis)
    evol1.data = numpy.abs(evol1.data)
    evol1.plot()    
    
    print("\nPlotting spectrum at frequency:", 
          fcont1.axis.data[show_Npoint2], units)
    sp21.plot(window=plot_window, Npos_contours=20, spart=qr.part_ABS)
    fftfile = "twod_fft_map_2.png"
    sp21.savefig(fftfile)
    print("... saved into: ", fftfile)

    evol2 = fcont1.get_point_evolution(12500,13000,
                                       fcont2.axis)
    evol2.data = numpy.abs(evol2.data)
    evol2.plot()

###############################################################################