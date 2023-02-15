# -*- coding: utf-8 -*-
"""

    Calculation of 2D spectra with lineshape functions


"""
import copy
import quantarhei as qr

from quantarhei import printlog as print
from quantarhei.spectroscopy import X

import matplotlib
matplotlib.rcParams['text.usetex'] = True

qr.assert_version(">","0.0.62")

print()
qr.timeit("Starting dimer simulation ...", show_stamp=True)

_show_plots_ = True 
_movie_ = False 
_save_2D_ = False

###############################################################################
#
# MODEL: Simple dimer of molecules
#
###############################################################################

Nt2 = 50
dt2 = 20
Npad = 0

Nt = Nt2
dt = dt2

t_axis = qr.TimeAxis(0.0, Nt2, dt2)

with qr.energy_units("1/cm"):
    # two two-level molecules
    m1 = qr.Molecule([0.0, 12000.0])
    m2 = qr.Molecule([0.0, 12300.0])

    
    #  correlation functions the environment
    #
    cfce_params1 = dict(ftype="OverdampedBrownian",
                   reorg=40.0,
                   cortime=100.0,
                   T=100, matsubara=20)
    cfce_params2 = dict(ftype="OverdampedBrownian",
                   reorg=40.0,
                   cortime=100.0,
                   T=100, matsubara=20)
    

    cfce1 = qr.CorrelationFunction(t_axis, cfce_params1)
    cfce2 = qr.CorrelationFunction(t_axis, cfce_params2)
    
    m1.set_transition_environment((0,1), cfce1)
    m2.set_transition_environment((0,1), cfce2)
    

# we create an aggregate from the two molecules
agg = qr.Aggregate(molecules=[m1, m2])

# we set transition dipole moment orientations for the two molecules
m1.set_dipole(0,1,[1.0, 0.8, 0.8])
m2.set_dipole(0,1,[0.8, 0.8, 0.0])

# resonance coupling is set by hand
with qr.energy_units("1/cm"):
    agg.set_resonance_coupling(0,1, 100.0)

# we copy the aggregate before it is built. For the calculation of 2D 
# spectrum, we need to build the aggregate so that it contains two-exciton
# states. But those are irrelevant for single exciton excited state dynamics
# so we make two identical aggregates, one with single-excitons only, and
# one with two-excitons. 
agg_2D = copy.copy(agg)


# the aggregate is built with single exciton states only
agg.build(mult=1)

# we print its Hamiltonian to check everything is alright
H = agg.get_Hamiltonian()
with qr.energy_units("1/cm"):
    print(H)


###############################################################################
#
# EXCITED STATE DYNAMICS: Lindblad relaxation between eigenstates
#
###############################################################################

# time span of the excited state evolution (later t2 time of the 2D spectrum)
t2_axis = qr.TimeAxis(0.0, Nt2, dt2)

# Lindblad relaxation operator
with qr.eigenbasis_of(H):
    K = qr.qm.ProjectionOperator(1,2,dim=H.dim)
rates = [1.0/200.0]

sbi = qr.qm.SystemBathInteraction(sys_operators=[K], rates=rates)

L = qr.qm.LindbladForm(H, sbi)


eUt = qr.EvolutionSuperOperator(time=t2_axis, ham=H, relt=L)
eUt.set_dense_dt(10)

eUt.calculate()

#if _show_plots_:
if False:
    with qr.eigenbasis_of(H):
        eUt.plot_element((2,2,2,2), show=False)
        eUt.plot_element((1,1,1,1), show=False)
        eUt.plot_element((1,1,2,2))
        eUt.plot_element((1,2,1,2))
   
    
###############################################################################
#
# 2D SPECTRUM: effective lineshape 2D spectrum
#
###############################################################################

# time axes of the propagation in t1 and t3 times

t1_axis = qr.TimeAxis(0.0, Nt, dt)
t3_axis = qr.TimeAxis(0.0, Nt, dt)

agg_2D.build(mult=2)
agg_2D.diagonalize()

# laboratory settings
lab = qr.LabSetup()
lab.set_polarizations(pulse_polarizations=(X,X,X), detection_polarization=X)

#
#
#
#
#
#from quantarhei.spectroscopy.mocktwodcalculator \
#    import MockTwoDResponseCalculator as TwoDResponseCalculator
#
#calc = TwoDResponseCalculator(t1_axis, t2_axis, t3_axis)
calc = qr.TwoDResponseCalculator(t1_axis, t2_axis, t3_axis, system=agg_2D)
with qr.energy_units("1/cm"):
    calc.bootstrap(rwa=12100.0, pad=Npad) 

#calc.bootstrap(rwa=qr.convert(12100.0,"1/cm","int"), pad=1000)

print("Calculating", Nt2,"spectra")
tcont = calc.calculate()
qr.done_in(True) 

tcont = tcont.get_TwoDSpectrumContainer()

T2 = 3*dt2
twod = tcont.get_spectrum(T2)

### short test of addition of data
twod1 = twod
twod2 = tcont.get_spectrum(0.0)

#twod.add_data(twod2.data)
### end of the test

if _save_2D_:
    twod.save("twod.qrp")

if _show_plots_:
    plot_window = [11500,13000,11500,13000]
    with qr.energy_units("1/cm"):
        twod.plot(Npos_contours=10, #indow=plot_window,            
                  stype=qr.signal_TOTL, spart=qr.part_REAL)
    qr.show_plot()
    
if _movie_:
    def label_func(sp):
        """Function returning label for each frame
        """
        return (r'$t_{2}='+str(sp.get_t2())+r'$ fs', [0.5,0.8])

    plot_window = [11500, 13000, 11500, 13000]
    with qr.energy_units("1/cm"):
        tcont.make_movie("twod.mp4", window=plot_window, label_func=label_func)
 
 
qr.stop() 


###############################################################################
#
# PUMP PROBE SPECRA as projections for 2D spectra
#
###############################################################################

# pump-probe spectrum can be obtained as a projection of a 2D spectrum
pprop = twod.get_PumpProbeSpectrum()

if _show_plots_:
    with qr.energy_units("1/cm"):
        pprop.plot(axis=[11500, 13000, -1200.0, 0.0], vmax=100,
                   title="Pump probe spectrum",
                   show=True)
 
# Pump-probe spectra also have their container
pcont = qr.PumpProbeSpectrumContainer(t2axis=t2_axis)
pcont.set_spectrum(pprop, tag=T2)

T2 = 300.0
twod = tcont.get_spectrum(T2)
pprop = twod.get_PumpProbeSpectrum()
pcont.set_spectrum(pprop, tag=T2)

# Container can be plotted (all curves in one plot)
if _show_plots_:
    with qr.energy_units("1/cm"):
        pcont.plot()
    qr.show_plot()
    
# Pump-probe spectra from 2D spectrum container
pcont2 = tcont.get_PumpProbeSpectrumContainer()

if _show_plots_:
    with qr.energy_units("1/cm"):
        pcont2.plot()

if _movie_:
    with qr.energy_units("1/cm"):   
        pcont2.make_movie("pprob.mp4")




