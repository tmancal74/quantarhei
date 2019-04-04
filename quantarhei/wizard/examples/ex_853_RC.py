# -*- coding: utf-8 -*-

import quantarhei as qr

use_vib = True
make_movie = True

# dimer of molecules
with qr.energy_units("1/cm"):
    mol1 = qr.Molecule([0.0, 10000.0])
    mol2 = qr.Molecule([0.0, 10500.0])
    
    mod2 = qr.Mode(500.0)
    


mol1.set_transition_width((0,1), qr.convert(200.0, "1/cm", "int"))
mol1.set_dipole(0,1,[1.0, 0.0, 0.0])

mol2.set_transition_width((0,1), qr.convert(200.0, "1/cm", "int"))
mol2.set_dipole(0,1,[1.0, 1.0, 0.0])

agg = qr.Aggregate([mol1, mol2])

with qr.energy_units("1/cm"):
    agg.set_resonance_coupling(0,1,100.0)

agg.save("agg.qrp")

if use_vib:
    mol2.add_Mode(mod2)
    mod2.set_nmax(0, 3)
    mod2.set_nmax(1, 3)
    mod2.set_HR(1, 0.01)

agg_el = qr.load_parcel("agg.qrp")

agg.save("agg_vib.qrp")
agg3 = qr.load_parcel("agg_vib.qrp")

agg.build(mult=1)
agg_el.build(mult=1)

HH = agg.get_Hamiltonian()
He = agg_el.get_Hamiltonian()



with qr.energy_units("1/cm"):
    print(HH)

#
# Laboratory setup
#

from quantarhei import LabSetup
from quantarhei.utils.vectors import X #, Y, Z

lab = LabSetup()
lab.set_polarizations(pulse_polarizations=[X,X,X], detection_polarization=X)

time2 = qr.TimeAxis(0.0, 100, 10.0)

cont = qr.TwoDSpectrumContainer(t2axis=time2)
#
# spectra will be indexed by the times in the time axis `time2`
#
cont.use_indexing_type(time2)

#
# We define two-time axes, which will be FFTed and will define the omega_1 and
# omega_3 axes of the 2D spectrum
#
t1_N_steps = 100
t1_time_step = 10.0
t3_N_steps = 100
t3_time_step = 10.0
t1axis = qr.TimeAxis(0.0, t1_N_steps, t1_time_step)
t3axis = qr.TimeAxis(0.0, t3_N_steps, t3_time_step)

#
# This calculator calculated 2D spectra from the effective width defined above
#
msc = qr.MockTwoDSpectrumCalculator(t1axis, time2, t3axis)
msc.bootstrap(rwa=qr.convert(10250.0,"1/cm","int"), 
              all_positive=False, shape="Gaussian")

operators = []
rates = []
with qr.eigenbasis_of(He):
    operators.append(qr.qm.ProjectionOperator(1, 2, dim=He.dim))
rates.append(1.0/200.0)
print(HH.dim)
print(operators[0])
#
# System-bath interaction including vibrational states
#
sbi = qr.qm.SystemBathInteraction(sys_operators=operators,
                                  rates=rates)
sbi.set_system(agg)

LF = qr.qm.ElectronicLindbladForm(HH, sbi, as_operators=True)

eUt = qr.qm.EvolutionSuperOperator(time2, HH, relt=LF, mode="all")
eUt.set_dense_dt(100)
eUt.calculate(show_progress=False)

agg3.build(mult=2)
agg3.diagonalize()

for t2 in time2.data:
    print(t2)
    twod = msc.calculate_one_system(t2, agg3, HH, eUt, lab)
    cont.set_spectrum(twod)

if make_movie:
    with qr.energy_units("1/cm"):
        cont.make_movie("mov.mp4")
    
#
# Window function for subsequenty FFT
#
import quantarhei.functions as func
window = func.Tukey(time2, r=0.3, sym=False)

#
# FFT with the window function
#
# Specify REPH, NONR or `total` to get different types of spectra
#
print("\nCalculating FFT of the 2D maps")
fcont = cont.fft(window=window, dtype="REPH", dpart="real", offset=0.0)

show_omega = 500.0

#
# Have a look which frequencies we actually have
#
Ndat = len(fcont.axis.data)
print("\nNumber of frequency points:", Ndat)
print("In 1/cm they are:")
with qr.energy_units("1/cm"):
    for k_i in range(Ndat):
        print(k_i, fcont.axis.data[k_i])

with qr.frequency_units("1/cm"):
    sp, show_Npoint = fcont.get_nearest(show_omega)

units = "1/cm"
with qr.energy_units(units):
    print("\nPlotting spectrum at frequency:", fcont.axis.data[show_Npoint], units)
    #sp.plot(Npos_contours=10, window=[9500,11500, 9500, 11500],
    #          stype="total", spart="abs")
    sp.plot(Npos_contours=10, 
            stype="total", spart="abs")   
    fftfile = "twod_fft_map.png"
    #sp.savefig(os.path.join(pre_out, fftfile))
    print("... saved into: ", fftfile)