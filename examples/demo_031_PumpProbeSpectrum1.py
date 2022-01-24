# -*- coding: utf-8 -*-
"""

    Calculation of pump-probe spectra with effective lineshapes 


"""
import copy
import quantarhei as qr
    
from quantarhei.spectroscopy import X

_show_plots_ = True
_movie_ = True

###############################################################################
#
# MODEL: Simple dimer of molecules
#
###############################################################################

with qr.energy_units("1/cm"):
    # two two-level molecules
    m1 = qr.Molecule([0.0, 12000.0])
    m2 = qr.Molecule([0.0, 12300.0])
    m3 = qr.Molecule([0.0, 12600.0])
    
    # transitions will have Gaussian lineshape with a width specified here
    m1.set_transition_width((0,1), 150.0) 
    m2.set_transition_width((0,1), 150.0)
    m3.set_transition_width((0,1), 200.0)

# we create an aggregate from the two molecules
agg = qr.Aggregate(molecules=[m1, m2, m3])

# we set transition dipole moment orientations for the two molecules
m1.set_dipole(0,1,[1.0, 0.8, 0.8])
m2.set_dipole(0,1,[0.8, 0.8, 0.0])
m3.set_dipole(0,1,[0.0, 2.0, 0.0])

# resonance coupling is set by hand
with qr.energy_units("1/cm"):
    agg.set_resonance_coupling(0,1, 100.0)
    agg.set_resonance_coupling(1,2, 50.0)

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
t2_axis = qr.TimeAxis(0.0, 300, 10.0)

# Lindblad relaxation operator
with qr.eigenbasis_of(H):
    K1 = qr.qm.ProjectionOperator(1,2,dim=H.dim)
    K2 = qr.qm.ProjectionOperator(2,3,dim=H.dim)
rates = [1.0/500.0, 1.0/100.0]

sbi = qr.qm.SystemBathInteraction(sys_operators=[K1, K2], rates=rates)

L = qr.qm.LindbladForm(H, sbi)


eUt = qr.EvolutionSuperOperator(time=t2_axis, ham=H, relt=L)
eUt.set_dense_dt(10)

with qr.eigenbasis_of(H):
    eUt.calculate()

if _show_plots_:
    with qr.eigenbasis_of(H):
        eUt.plot_element((2,2,2,2), show=False)
        eUt.plot_element((1,1,1,1), show=False)
        eUt.plot_element((1,1,2,2))
        eUt.plot_element((1,2,1,2))
    
    
###############################################################################
#
# PUMP-PROBE SPECTRUM: effective lineshape pump-probe spectrum
#
###############################################################################

# time axes of the propagation in t1 and t3 times

t1_axis = qr.TimeAxis(0.0, 100, 10.0)
t3_axis = qr.TimeAxis(0.0, 100, 10.0)

agg_2D.build(mult=2)
agg_2D.diagonalize()

# laboratory settings
lab = qr.LabSetup()
lab.set_polarizations(pulse_polarizations=(X,X,X), detection_polarization=X)

#
# Given a molecular system (only Aggregate class so far), we calculate
# 2D spectrum along the whole t2 axis in one shot
#
pcalc = qr.MockPumpProbeSpectrumCalculator(t1_axis, t2_axis, t3_axis)
with qr.energy_units("1/cm"):
    pcalc.bootstrap(rwa=12100.0)
    
pcont3 = pcalc.calculate_all_system(agg_2D, H, eUt, lab)

if _movie_:
    with qr.energy_units("1/cm"):   
        pcont3.make_movie("pprob2.mp4") #, axis=[10500, 13500, -600, 100])