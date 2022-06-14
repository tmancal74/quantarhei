# -*- coding: utf-8 -*-
import numpy
import quantarhei as qr

_show_plots_ = True
_save_2D_ = True
_use_disorder_ = True

E0 = 12000.0
with qr.energy_units("1/cm"):
    # two two-level molecules
    m1 = qr.Molecule([0.0, E0])
    
    # transitions will have Gaussian lineshape with a width specified here
    m1.set_transition_width((0,1), 100.0) 

# we create an aggregate from the two molecules
agg = qr.Aggregate(molecules=[m1])

# we set transition dipole moment orientations for the two molecules
m1.set_dipole(0,1,[1.0, 0.8, 0.8])

# time axes of the propagation in t1 and t3 times
t2_axis = qr.TimeAxis(0.0, 100, 10.0)
t1_axis = qr.TimeAxis(0.0, 100, 10.0)
t3_axis = qr.TimeAxis(0.0, 100, 10.0)

from quantarhei.spectroscopy.mocktwodcalculator \
    import MockTwoDResponseCalculator as TwoDResponseCalculator
    
from quantarhei.spectroscopy import X

calc = TwoDResponseCalculator(t1_axis, t2_axis, t3_axis)
with qr.energy_units("1/cm"):
    calc.bootstrap(rwa=12100.0) #qr.convert(12100.0,"1/cm","int"))

agg_1 = agg.deepcopy()
agg_1.build(mult=1)
H = agg_1.get_Hamiltonian()

agg.build(mult=2)
agg.diagonalize()

# laboratory settings
lab = qr.LabSetup()
lab.set_polarizations(pulse_polarizations=(X,X,X), detection_polarization=X)

eUt = qr.EvolutionSuperOperator(t2_axis, H)
eUt.set_dense_dt(10)

eUt.calculate()

tcont = calc.calculate_all_system(agg, eUt, lab)

T2 = 0.0
twod = tcont.get_spectrum(T2)

twod = twod.get_TwoDSpectrum()

if _save_2D_:
    twod.save("twod.qrp")

if _show_plots_:
    plot_window = [11500,13000,11500,13000]
    with qr.energy_units("1/cm"):
        twod.plot(Npos_contours=10, window=plot_window,            
                  stype=qr.signal_TOTL, spart=qr.part_REAL)
        


#
# Apply disorder averaging = integration over E
#
if _use_disorder_:
    twod_dis = twod.deepcopy()
    twod_dis.data[:,:] = 0.0
    twod_aux = twod.deepcopy()
    twod_aux.data[:,:] = 0.0
    
    N_dis = 1000
    dE_dis = 1.0
    with qr.energy_units("1/cm"):
        for i_dis in range(N_dis):
            E = 12250.0 + dE_dis*(-N_dis/2 + i_dis)
            E_shift = E - E0
            
            #print(i_dis, E, E_shift)
            
            twod_aux.data[:,:] = twod.data[:,:]
            twod_aux.shift_energy(E_shift)
        
            twod_dis.data[:,:] += \
            numpy.exp(-((E-12250.0)/100.0)**2)*twod_aux.data[:,:]
            
    if _show_plots_:
        plot_window = [11500,13000,11500,13000]
        with qr.energy_units("1/cm"):
            twod_dis.plot(Npos_contours=10, window=plot_window,            
                      stype=qr.signal_TOTL, spart=qr.part_REAL)
    
            



